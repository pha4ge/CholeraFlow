/* Screen 1 — Upload & Validate. All real logic, no mock data: classifies
   dropped files, detects paired-end mates from filenames, normalises sample
   names, parses metadata CSVs, and counts reads/contigs client-side. */
(function () {
  var ACCEPTED_EXT = ['.fastq.gz', '.fasta', '.csv'];
  var PAIR_SUFFIX_RE = /_(R)?([12])$/i;
  var SAMPLE_COL_CANDIDATES = ['sample', 'sample_name', 'sample_id', 'name'];

  // ---- module state -------------------------------------------------
  var seqFilesByName = {};        // filename -> File   (.fastq.gz / .fasta)
  var rejectedFilesByName = {};   // filename -> File   (unsupported extension)
  var metadataFile = null;
  var metadataHeaders = [];
  var metadataRecords = [];
  var metadataSampleColumn = null;
  var metadataReadError = null;
  var manualMetadata = { date: '', source: 'Clinical isolate', location: '' };
  var metaTab = 'csv';
  var readCountCache = {}; // filename -> Promise<number|null>

  // ---- file classification ------------------------------------------
  function classifyExtension(filename) {
    var lower = filename.toLowerCase();
    if (lower.endsWith('.fastq.gz')) return 'fastq';
    if (lower.endsWith('.fasta')) return 'fasta';
    if (lower.endsWith('.csv')) return 'csv';
    return null;
  }

  function stripExt(filename, ext) {
    return filename.slice(0, filename.length - ext.length);
  }

  function splitPairSuffix(basename) {
    var m = basename.match(PAIR_SUFFIX_RE);
    if (!m) return { base: basename, mate: null, rPrefix: '' };
    return { base: basename.slice(0, m.index), mate: m[2], rPrefix: m[1] || '' };
  }

  function guessMateName(filename, wantMate) {
    var basename = stripExt(filename, '.fastq.gz');
    var split = splitPairSuffix(basename);
    if (!split.mate) return null;
    return split.base + '_' + split.rPrefix + wantMate + '.fastq.gz';
  }

  function normalizeName(base) {
    return base.trim().replace(/[^A-Za-z0-9._-]+/g, '_').replace(/^_+|_+$/g, '');
  }

  // ---- grouping into samples ------------------------------------------
  function buildSequenceGroups() {
    var groups = {};
    var order = [];
    Object.keys(seqFilesByName).forEach(function (name) {
      var file = seqFilesByName[name];
      var ext = classifyExtension(name);

      if (ext === 'fasta') {
        var fastaBase = stripExt(name, '.fasta');
        var fastaKey = 'fasta:' + fastaBase;
        if (!groups[fastaKey]) {
          groups[fastaKey] = { key: fastaKey, kind: 'fasta', base: fastaBase, files: [] };
          order.push(fastaKey);
        }
        groups[fastaKey].files.push(file);
        return;
      }

      if (ext === 'fastq') {
        var basename = stripExt(name, '.fastq.gz');
        var split = splitPairSuffix(basename);
        var key = 'fastq:' + split.base;
        if (!groups[key]) {
          groups[key] = { key: key, kind: 'fastq', base: split.base, files: [], mate1: null, mate2: null, single: null };
          order.push(key);
        }
        var g = groups[key];
        g.files.push(file);
        if (split.mate === '1') {
          if (g.mate1) g.duplicateMate1 = true; else g.mate1 = file;
        } else if (split.mate === '2') {
          if (g.mate2) g.duplicateMate2 = true; else g.mate2 = file;
        } else {
          if (g.single) g.duplicateSingle = true; else g.single = file;
        }
      }
    });
    return order.map(function (key) { return groups[key]; });
  }

  function evaluateGroup(g) {
    var issues = [];
    var pairStatus;

    if (g.kind === 'fasta') {
      pairStatus = 'fasta';
      if (g.files.length > 1) {
        issues.push('"' + g.base + '" has ' + g.files.length + ' FASTA files with the same sample name (' +
          g.files.map(function (f) { return f.name; }).join(', ') + '). Each sample needs a unique file name — rename the extra file(s) and re-upload.');
      }
    } else {
      var hasM1 = !!g.mate1, hasM2 = !!g.mate2, hasSingle = !!g.single;
      if (hasM1 && hasM2 && !hasSingle) {
        pairStatus = 'paired';
      } else if (hasSingle && !hasM1 && !hasM2) {
        pairStatus = 'single';
      } else if (hasM1 && !hasM2 && !hasSingle) {
        pairStatus = 'unmatched';
        var expected2 = guessMateName(g.mate1.name, '2');
        issues.push('"' + g.mate1.name + '" looks like a paired-end R1 file but no matching R2 file was found. ' +
          'CholeraFlow expected a file named "' + expected2 + '". Upload the missing mate, or rename this file to remove the _R1/_1 suffix if it is actually single-end.');
      } else if (hasM2 && !hasM1 && !hasSingle) {
        pairStatus = 'unmatched';
        var expected1 = guessMateName(g.mate2.name, '1');
        issues.push('"' + g.mate2.name + '" looks like a paired-end R2 file but no matching R1 file was found. ' +
          'CholeraFlow expected a file named "' + expected1 + '". Upload the missing mate, or rename this file to remove the _R2/_2 suffix if it is actually single-end.');
      } else {
        pairStatus = 'conflict';
        issues.push('"' + g.base + '" has conflicting files (' + g.files.map(function (f) { return f.name; }).join(', ') +
          '). CholeraFlow can\'t tell whether this is one paired-end sample or several single-end samples. Rename the files so each sample has either one file, or exactly two files using a consistent _R1/_R2 (or _1/_2) suffix.');
      }
      if (g.duplicateMate1) issues.push('More than one file matches "' + g.base + '" as the R1/_1 mate. Remove the duplicate before continuing.');
      if (g.duplicateMate2) issues.push('More than one file matches "' + g.base + '" as the R2/_2 mate. Remove the duplicate before continuing.');
      if (g.duplicateSingle) issues.push('More than one file resolves to the single-end sample "' + g.base + '". Rename one of them so each sample has a unique name.');
    }

    if (!g.base) {
      issues.push('A file doesn\'t have a usable sample name before its pairing suffix. Rename it to include a sample identifier, e.g. "sample1_R1.fastq.gz".');
    }

    return { pairStatus: pairStatus, issues: issues };
  }

  function buildSamples() {
    var groups = buildSequenceGroups();
    var samples = groups.map(function (g) {
      var evald = evaluateGroup(g);
      return {
        id: g.key,
        rawBase: g.base,
        normalizedName: normalizeName(g.base),
        kind: g.kind,
        group: g,
        pairStatus: evald.pairStatus,
        issues: evald.issues,
      };
    });

    // name collisions across the whole batch (case-insensitive)
    var seen = {};
    samples.forEach(function (s) {
      var key = s.normalizedName.toLowerCase();
      if (!key) return;
      (seen[key] = seen[key] || []).push(s);
    });
    Object.keys(seen).forEach(function (key) {
      var group = seen[key];
      if (group.length > 1) {
        var names = group.map(function (s) { return '"' + (s.group.files[0] ? s.group.files[0].name : s.rawBase) + '"'; }).join(' and ');
        group.forEach(function (s) {
          s.issues.push('Sample name "' + s.normalizedName + '" is used by more than one upload (' + names +
            '). Rename one of the source files so each sample has a unique name.');
        });
      }
    });

    return samples;
  }

  // ---- CSV parsing -----------------------------------------------------
  function parseCsv(text) {
    var rows = [];
    var i = 0, field = '', row = [], inQuotes = false;
    var len = text.length;
    function pushField() { row.push(field); field = ''; }
    function pushRow() { pushField(); rows.push(row); row = []; }
    while (i < len) {
      var c = text[i];
      if (inQuotes) {
        if (c === '"') {
          if (text[i + 1] === '"') { field += '"'; i += 2; continue; }
          inQuotes = false; i++; continue;
        }
        field += c; i++; continue;
      } else {
        if (c === '"') { inQuotes = true; i++; continue; }
        if (c === ',') { pushField(); i++; continue; }
        if (c === '\r') { i++; continue; }
        if (c === '\n') { pushRow(); i++; continue; }
        field += c; i++; continue;
      }
    }
    if (field.length > 0 || row.length > 0) pushRow();
    while (rows.length && rows[rows.length - 1].length === 1 && rows[rows.length - 1][0] === '') rows.pop();
    return rows;
  }

  function rowsToRecords(rows) {
    if (!rows.length) return { headers: [], records: [] };
    var headers = rows[0].map(function (h) { return h.trim(); });
    var records = rows.slice(1).map(function (r) {
      var obj = {};
      headers.forEach(function (h, idx) { obj[h] = (r[idx] !== undefined ? r[idx] : '').trim(); });
      return obj;
    });
    return { headers: headers, records: records };
  }

  function findSampleColumn(headers) {
    var lower = headers.map(function (h) { return h.toLowerCase(); });
    for (var i = 0; i < SAMPLE_COL_CANDIDATES.length; i++) {
      var idx = lower.indexOf(SAMPLE_COL_CANDIDATES[i]);
      if (idx !== -1) return headers[idx];
    }
    return null;
  }

  // ---- read / contig counting ------------------------------------------
  function countFastqReads(file) {
    if (typeof DecompressionStream === 'undefined') return Promise.resolve(null);
    try {
      var stream = file.stream().pipeThrough(new DecompressionStream('gzip'));
      var reader = stream.getReader();
      var decoder = new TextDecoder();
      var newlineCount = 0;
      function pump() {
        return reader.read().then(function (res) {
          if (res.done) return Math.floor(newlineCount / 4);
          var text = decoder.decode(res.value, { stream: true });
          for (var i = 0; i < text.length; i++) if (text.charCodeAt(i) === 10) newlineCount++;
          return pump();
        });
      }
      return pump().catch(function () { return null; });
    } catch (e) {
      return Promise.resolve(null);
    }
  }

  function countFastaContigs(file) {
    try {
      var reader = file.stream().getReader();
      var decoder = new TextDecoder();
      var count = 0, atLineStart = true;
      function pump() {
        return reader.read().then(function (res) {
          if (res.done) return count;
          var text = decoder.decode(res.value, { stream: true });
          for (var i = 0; i < text.length; i++) {
            var ch = text[i];
            if (atLineStart && ch === '>') count++;
            atLineStart = (ch === '\n');
          }
          return pump();
        });
      }
      return pump().catch(function () { return null; });
    } catch (e) {
      return Promise.resolve(null);
    }
  }

  function readCountFor(file) {
    if (!readCountCache[file.name]) {
      var ext = classifyExtension(file.name);
      readCountCache[file.name] = ext === 'fastq' ? countFastqReads(file) : countFastaContigs(file);
    }
    return readCountCache[file.name];
  }

  // ---- removing a mistaken upload (the only way to recover without a
  // full page reload, since filenames are the de-dup key) ----------------
  function removeSampleById(id) {
    var target = buildSamples().filter(function (s) { return s.id === id; })[0];
    if (!target) return;
    target.group.files.forEach(function (f) { delete seqFilesByName[f.name]; delete readCountCache[f.name]; });
  }

  function removeRejected(name) {
    delete rejectedFilesByName[name];
  }

  // ---- DOM refs (resolved on init) -------------------------------------
  var dom = {};

  function q(id) { return document.getElementById(id); }

  function cacheDom() {
    dom.dropzone = q('cf-dropzone');
    dom.fileInput = q('cf-file-input');
    dom.metaTabCsv = q('cf-meta-tab-csv');
    dom.metaTabManual = q('cf-meta-tab-manual');
    dom.metaCsvPanel = q('cf-meta-csv-panel');
    dom.metaManualPanel = q('cf-meta-manual-panel');
    dom.metaDropzone = q('cf-meta-dropzone');
    dom.metaFileInput = q('cf-meta-file-input');
    dom.metaStatus = q('cf-meta-status');
    dom.metaStatusName = q('cf-meta-status-name');
    dom.metaRemoveBtn = q('cf-meta-remove-btn');
    dom.metaError = q('cf-meta-error');
    dom.metaTemplateLink = q('cf-meta-template-link');
    dom.metaDate = q('cf-meta-date');
    dom.metaSource = q('cf-meta-source');
    dom.metaLocation = q('cf-meta-location');
    dom.metaPreviewCard = q('cf-meta-preview-card');
    dom.metaPreviewCount = q('cf-meta-preview-count');
    dom.metaPreviewTableWrap = q('cf-meta-preview-table-wrap');
    dom.metaRemoveBtn2 = q('cf-meta-remove-btn-2');
    dom.validationCard = q('cf-validation-card');
    dom.checksPassedLabel = q('cf-checks-passed-label');
    dom.checksList = q('cf-checks-list');
    dom.sampleCard = q('cf-sample-card');
    dom.sampleCount = q('cf-sample-count');
    dom.sampleRows = q('cf-sample-rows');
    dom.errorBanner = q('cf-error-banner');
    dom.continueHint = q('cf-continue-hint');
    dom.continueBtn = q('cf-continue-btn');
  }

  // ---- rendering --------------------------------------------------------
  function fileCellFor(sample) {
    var g = sample.group;
    if (g.kind === 'fasta') return g.files.map(function (f) { return f.name; }).join(', ');
    if (g.mate1 && g.mate2) return g.mate1.name + ' +1';
    if (g.single) return g.single.name;
    return g.files.map(function (f) { return f.name; }).join(', ');
  }

  function statusBadge(sample) {
    var ok = sample.issues.length === 0;
    return '<span class="cf-pill cf-pill--badge ' + (ok ? 'cf-pill--success' : 'cf-pill--error') + '">' +
      (ok ? 'Accepted' : 'Needs fix') + '</span>' +
      '<button type="button" class="cf-row-remove" data-action="remove-sample" data-key="' + escapeHtml(sample.id) + '" ' +
      'style="display:block;margin:4px 0 0 auto;background:none;border:none;color:var(--cf-text-faint);font-size:11px;cursor:pointer;text-decoration:underline">Remove</button>';
  }

  function renderSampleRows(samples) {
    dom.sampleRows.innerHTML = samples.map(function (sample) {
      var rowClass = sample.issues.length ? ' style="background:var(--cf-error-bg)"' : '';
      return '<tr data-sample-id="' + escapeHtml(sample.id) + '"' + rowClass + '>' +
        '<td class="cf-mono" style="font-weight:600;color:var(--cf-text)">' + escapeHtml(sample.normalizedName || '(unnamed)') + '</td>' +
        '<td class="cf-mono" style="color:var(--cf-text-muted)">' + escapeHtml(fileCellFor(sample)) + '</td>' +
        '<td class="cf-mono js-reads" style="color:#475569"><span class="cf-spinner"></span></td>' +
        '<td class="cf-table__right">' + statusBadge(sample) + '</td>' +
        '</tr>';
    }).join('');

    samples.forEach(function (sample) {
      var g = sample.group;
      var promises;
      if (g.kind === 'fasta') {
        promises = [readCountFor(g.files[0])];
      } else if (g.mate1 && g.mate2) {
        promises = [readCountFor(g.mate1), readCountFor(g.mate2)];
      } else if (g.single) {
        promises = [readCountFor(g.single)];
      } else {
        promises = g.files.map(readCountFor);
      }
      Promise.all(promises).then(function (counts) {
        var cell = dom.sampleRows.querySelector('tr[data-sample-id="' + cssEscape(sample.id) + '"] .js-reads');
        if (!cell) return;
        var valid = counts.filter(function (c) { return typeof c === 'number'; });
        if (!valid.length) { cell.textContent = '—'; return; }
        var total = valid.reduce(function (a, b) { return a + b; }, 0);
        if (g.kind === 'fasta') {
          cell.textContent = total + (total === 1 ? ' contig' : ' contigs');
        } else {
          cell.textContent = formatCount(total);
        }
      });
    });
  }

  function formatCount(n) {
    if (n >= 1e6) return (n / 1e6).toFixed(1) + 'M';
    if (n >= 1e3) return (n / 1e3).toFixed(1) + 'K';
    return String(n);
  }

  function escapeHtml(s) {
    return String(s).replace(/[&<>"']/g, function (c) {
      return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c];
    });
  }

  function cssEscape(s) {
    return String(s).replace(/["\\]/g, '\\$&');
  }

  function checkRowHtml(check) {
    var icon;
    if (check.status === 'pass') icon = '<span class="cf-check-icon cf-check-icon--pass">✓</span>';
    else if (check.status === 'fail') icon = '<span class="cf-check-icon cf-check-icon--fail">✕</span>';
    else if (check.status === 'pending') icon = '<span class="cf-spinner"></span>';
    else icon = '<span class="cf-check-icon" style="background:var(--cf-neutral-bg);color:var(--cf-text-faint)">–</span>';

    var textClass = check.status === 'pass' ? 'cf-check-text--pass' : check.status === 'fail' ? 'cf-check-text--fail' : 'cf-check-text--pending';
    var note = '';
    if (check.status === 'pass' && check.note) note = '<span class="cf-check-note--pass">' + escapeHtml(check.note) + '</span>';
    if (check.status === 'neutral' && check.note) note = '<span style="font-size:11.5px;color:var(--cf-text-faint)">' + escapeHtml(check.note) + '</span>';

    var fixes = '';
    if (check.status === 'fail' && check.fixMessages && check.fixMessages.length) {
      fixes = check.fixMessages.map(function (m) {
        return '<div class="cf-check-fix">⚠ ' + escapeHtml(m) + '</div>';
      }).join('');
    }

    return '<div><div class="cf-check-row">' + icon +
      '<span class="' + textClass + '">' + escapeHtml(check.label) + '</span>' + note + '</div>' + fixes + '</div>';
  }

  function recomputeAndRender() {
    var samples = buildSamples();
    var rejectedFiles = Object.keys(rejectedFilesByName).map(function (n) { return rejectedFilesByName[n]; });
    var hasAnyUpload = samples.length > 0 || rejectedFiles.length > 0 || !!metadataFile;

    dom.validationCard.style.display = hasAnyUpload ? '' : 'none';
    dom.sampleCard.style.display = samples.length > 0 ? '' : 'none';

    if (samples.length > 0) {
      dom.sampleCount.textContent = '· ' + samples.length + (samples.length === 1 ? ' sample' : ' samples');
      renderSampleRows(samples);
    }

    // ---- checks ----
    var formatOk = rejectedFiles.length === 0;
    var formatFixes = rejectedFiles.map(function (f) {
      return '"' + f.name + '" has an unsupported extension. CholeraFlow accepts .fastq.gz for reads, .fasta for assemblies, and .csv for sample metadata. Rename or re-export the file and try again.';
    });
    var formatCheck = {
      label: 'File format valid (FASTQ / FASTA)',
      status: !hasAnyUpload ? 'neutral' : (formatOk ? 'pass' : 'fail'),
      note: formatOk && samples.length ? (samples.length + (samples.length === 1 ? ' file group' : ' file groups')) : null,
      fixMessages: formatFixes,
    };

    var pairIssueSamples = samples.filter(function (s) { return s.pairStatus === 'unmatched' || s.pairStatus === 'conflict'; });
    var pairFixes = [];
    pairIssueSamples.forEach(function (s) { pairFixes = pairFixes.concat(s.issues); });
    var pairedCount = samples.filter(function (s) { return s.pairStatus === 'paired'; }).length;
    var pairCheck = {
      label: 'Paired-end reads matched',
      status: samples.length === 0 ? 'neutral' : (pairFixes.length === 0 ? 'pass' : 'fail'),
      note: pairFixes.length === 0 && pairedCount ? (pairedCount + (pairedCount === 1 ? ' pair' : ' pairs')) : null,
      fixMessages: pairFixes,
    };

    var nameFixes = [];
    samples.forEach(function (s) {
      s.issues.forEach(function (issue) {
        if (issue.indexOf('Sample name') === 0 || issue.indexOf('usable sample name') !== -1) nameFixes.push(issue);
      });
    });
    var nameCheck = {
      label: 'Sample names normalised',
      status: samples.length === 0 ? 'neutral' : (nameFixes.length === 0 ? 'pass' : 'fail'),
      note: nameFixes.length === 0 && samples.length ? (samples.length + ' unique names') : null,
      fixMessages: nameFixes,
    };

    var metaCheck = buildMetadataCheck(samples);

    var checks = [formatCheck, pairCheck, nameCheck, metaCheck];
    var passed = checks.filter(function (c) { return c.status === 'pass'; }).length;

    dom.checksPassedLabel.textContent = passed + ' / ' + checks.length + ' checks passed';
    dom.checksList.innerHTML = checks.map(checkRowHtml).join('');

    // ---- aggregate error banner (specific fix messages, each removable) ----
    var bannerEntries = rejectedFiles.map(function (f) {
      return {
        message: '"' + f.name + '" has an unsupported extension. CholeraFlow accepts .fastq.gz for reads, .fasta for assemblies, and .csv for sample metadata. Rename or re-export the file and try again.',
        action: 'remove-rejected', key: f.name,
      };
    });
    samples.forEach(function (s) {
      s.issues.forEach(function (msg) { bannerEntries.push({ message: msg, action: 'remove-sample', key: s.id }); });
    });

    if (bannerEntries.length) {
      dom.errorBanner.style.display = '';
      dom.errorBanner.innerHTML =
        '<div class="cf-error-banner__title"><span class="cf-check-icon cf-check-icon--fail">✕</span>' +
        (bannerEntries.length === 1 ? '1 issue needs to be fixed before you can continue' : bannerEntries.length + ' issues need to be fixed before you can continue') +
        '</div>' +
        bannerEntries.map(function (entry) {
          return '<div style="display:flex;align-items:flex-start;justify-content:space-between;gap:14px;font-size:12.5px;color:var(--cf-error-text-strong);line-height:1.5">' +
            '<span>' + escapeHtml(entry.message) + '</span>' +
            '<button type="button" data-action="' + entry.action + '" data-key="' + escapeHtml(entry.key) + '" ' +
            'style="flex:none;background:none;border:none;color:var(--cf-error-text-strong);text-decoration:underline;font-size:12px;font-weight:600;cursor:pointer;white-space:nowrap">Remove</button>' +
            '</div>';
        }).join('');
    } else {
      dom.errorBanner.style.display = 'none';
      dom.errorBanner.innerHTML = '';
    }

    var allValid = formatOk && pairFixes.length === 0 && nameFixes.length === 0 && samples.length > 0;
    dom.continueBtn.disabled = !allValid;

    if (allValid) {
      dom.continueHint.textContent = 'All checks must pass before you can continue.';
      dom.continueHint.style.color = 'var(--cf-text-faint)';
    } else if (samples.length === 0 && hasAnyUpload) {
      dom.continueHint.textContent = 'Add at least one FASTQ.GZ or FASTA sequencing file to continue — metadata alone isn\'t enough.';
      dom.continueHint.style.color = 'var(--cf-error-text-strong)';
    } else if (samples.length === 0) {
      dom.continueHint.textContent = 'Drag and drop sequencing files above to get started.';
      dom.continueHint.style.color = 'var(--cf-text-faint)';
    } else {
      dom.continueHint.textContent = 'All checks must pass before you can continue.';
      dom.continueHint.style.color = 'var(--cf-text-faint)';
    }

    persistState(samples, metaCheck);
  }

  function buildMetadataCheck(samples) {
    var manualProvided = manualMetadata.date.trim() !== '' || manualMetadata.location.trim() !== '';
    if (metaTab === 'manual') {
      return {
        label: 'Metadata matched to samples',
        status: manualProvided ? 'pass' : 'neutral',
        note: manualProvided ? 'Manual metadata applied to all samples' : 'Optional — not provided',
        fixMessages: [],
      };
    }

    if (metadataReadError) {
      return { label: 'Metadata matched to samples', status: 'fail', fixMessages: [metadataReadError] };
    }
    if (!metadataFile) {
      return { label: 'Metadata matched to samples', status: 'neutral', note: 'Optional — not provided', fixMessages: [] };
    }
    if (!metadataSampleColumn) {
      return {
        label: 'Metadata matched to samples', status: 'fail',
        fixMessages: ['The metadata CSV doesn\'t have a recognisable sample column. Add a column named "sample" (or "sample_name" / "sample_id") so CholeraFlow can match rows to your uploaded files.'],
      };
    }

    var sampleNames = {};
    samples.forEach(function (s) { if (s.normalizedName) sampleNames[s.normalizedName.toLowerCase()] = false; });

    var unmatchedRows = [];
    metadataRecords.forEach(function (rec) {
      var raw = rec[metadataSampleColumn] || '';
      var norm = normalizeName(raw).toLowerCase();
      if (norm && sampleNames.hasOwnProperty(norm)) sampleNames[norm] = true;
      else unmatchedRows.push(raw || '(blank)');
    });

    var unmatchedSamples = Object.keys(sampleNames).filter(function (k) { return !sampleNames[k]; });

    if (!samples.length) {
      return { label: 'Metadata matched to samples', status: 'neutral', note: 'Upload sequencing files first', fixMessages: [] };
    }

    var fixes = [];
    unmatchedRows.forEach(function (r) { fixes.push('Metadata row "' + r + '" doesn\'t match any uploaded sample. Check that the sample column matches the file-derived sample name exactly.'); });
    unmatchedSamples.forEach(function (s) { fixes.push('Sample "' + s + '" has no matching row in the metadata CSV — it will be analysed without metadata.'); });

    return {
      label: 'Metadata matched to samples',
      status: fixes.length === 0 ? 'pass' : 'fail',
      note: fixes.length === 0 ? (metadataRecords.length + ' / ' + samples.length + ' matched') : null,
      fixMessages: fixes,
    };
  }

  function persistState(samples, metaCheck) {
    var summary = samples.map(function (s) {
      return { name: s.normalizedName, kind: s.kind, pairStatus: s.pairStatus, files: fileCellFor(s), valid: s.issues.length === 0 };
    });
    window.CholeraFlowState.update({
      upload: {
        samples: summary,
        metadataMode: metaTab,
        metadataProvided: metaCheck.status !== 'neutral',
        manualMetadata: manualMetadata,
      },
    });
  }

  // ---- metadata preview --------------------------------------------------
  // The raw table renders full-width (its own card, below the sequencing
  // sample preview) rather than inside the narrow metadata column — a CSV
  // with long unwrappable content (e.g. FTP/HTTP URLs in a real CholeraSeq
  // samplesheet) would otherwise force that 1fr grid column to grow and
  // squeeze the dropzone column next to it.
  function renderMetaPreview(errorMessage) {
    if (errorMessage) {
      dom.metaStatus.style.display = 'none';
      dom.metaPreviewCard.style.display = 'none';
      dom.metaError.style.display = '';
      dom.metaError.innerHTML = '<div class="cf-check-fix">⚠ ' + escapeHtml(errorMessage) + '</div>';
      return;
    }
    dom.metaError.style.display = 'none';
    dom.metaError.innerHTML = '';

    if (!metadataFile) {
      dom.metaStatus.style.display = 'none';
      dom.metaPreviewCard.style.display = 'none';
      return;
    }

    dom.metaStatus.style.display = 'flex';
    dom.metaStatusName.textContent = metadataFile.name + ' · ' + metadataRecords.length + (metadataRecords.length === 1 ? ' row' : ' rows');

    var rowCap = 50;
    var rowsToShow = metadataRecords.slice(0, rowCap);
    dom.metaPreviewCount.textContent = '· ' + metadataFile.name + ' · ' + metadataRecords.length + (metadataRecords.length === 1 ? ' row' : ' rows');
    var html = '<table class="cf-table"><thead><tr>' +
      metadataHeaders.map(function (h) { return '<th>' + escapeHtml(h) + '</th>'; }).join('') + '</tr></thead><tbody>';
    html += rowsToShow.map(function (rec) {
      return '<tr>' + metadataHeaders.map(function (h) { return '<td class="cf-mono" style="white-space:nowrap">' + escapeHtml(rec[h] || '') + '</td>'; }).join('') + '</tr>';
    }).join('');
    html += '</tbody></table>';
    if (metadataRecords.length > rowCap) {
      html += '<div style="padding:10px 18px;font-size:12px;color:var(--cf-text-faint)">+' + (metadataRecords.length - rowCap) + ' more rows not shown</div>';
    }
    dom.metaPreviewTableWrap.innerHTML = html;
    dom.metaPreviewCard.style.display = '';
  }

  function clearMetadata() {
    metadataFile = null;
    metadataHeaders = [];
    metadataRecords = [];
    metadataSampleColumn = null;
    metadataReadError = null;
    renderMetaPreview();
    recomputeAndRender();
  }

  function handleMetadataCsv(file) {
    metadataFile = file;
    metadataReadError = null;
    var reader = new FileReader();
    reader.onload = function () {
      var rows = parseCsv(String(reader.result));
      var parsed = rowsToRecords(rows);
      metadataHeaders = parsed.headers;
      metadataRecords = parsed.records;
      metadataSampleColumn = findSampleColumn(parsed.headers);
      renderMetaPreview();
      recomputeAndRender();
    };
    reader.onerror = function () {
      metadataFile = null;
      metadataHeaders = [];
      metadataRecords = [];
      metadataReadError = 'Could not read "' + file.name + '". Make sure it is a plain-text CSV file and try again.';
      renderMetaPreview(metadataReadError);
      recomputeAndRender();
    };
    reader.readAsText(file);
  }

  // ---- drop / input wiring ------------------------------------------------
  function preventDefaults(e) { e.preventDefault(); e.stopPropagation(); }

  function handleIncomingFiles(files) {
    var csvFile = null;
    files.forEach(function (f) {
      var ext = classifyExtension(f.name);
      if (ext === 'fastq' || ext === 'fasta') {
        seqFilesByName[f.name] = f;
      } else if (ext === 'csv') {
        csvFile = f;
      } else {
        rejectedFilesByName[f.name] = f;
      }
    });
    if (csvFile) {
      setMetaTab('csv');
      handleMetadataCsv(csvFile);
    }
    recomputeAndRender();
  }

  function setMetaTab(tab) {
    metaTab = tab;
    dom.metaTabCsv.classList.toggle('cf-seg__btn--active', tab === 'csv');
    dom.metaTabManual.classList.toggle('cf-seg__btn--active', tab === 'manual');
    dom.metaCsvPanel.style.display = tab === 'csv' ? '' : 'none';
    dom.metaManualPanel.style.display = tab === 'manual' ? '' : 'none';
    recomputeAndRender();
  }

  function downloadTemplate(e) {
    e.preventDefault();
    var csv = 'sample,collection_date,source_type,location\n';
    var blob = new Blob([csv], { type: 'text/csv' });
    var url = URL.createObjectURL(blob);
    var a = document.createElement('a');
    a.href = url;
    a.download = 'choleraflow_metadata_template.csv';
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }

  function wireDropzone(zone, input, onFiles) {
    ['dragenter', 'dragover'].forEach(function (evt) {
      zone.addEventListener(evt, function (e) { preventDefaults(e); zone.classList.add('cf-dropzone--active'); });
    });
    ['dragleave', 'drop'].forEach(function (evt) {
      zone.addEventListener(evt, function (e) { preventDefaults(e); zone.classList.remove('cf-dropzone--active'); });
    });
    zone.addEventListener('drop', function (e) {
      var files = Array.prototype.slice.call((e.dataTransfer && e.dataTransfer.files) || []);
      if (files.length) onFiles(files);
    });
    zone.addEventListener('click', function () { input.click(); });
    zone.addEventListener('keydown', function (e) {
      if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); input.click(); }
    });
    input.addEventListener('change', function (e) {
      var files = Array.prototype.slice.call(e.target.files || []);
      if (files.length) onFiles(files);
      input.value = '';
    });
  }

  function init() {
    cacheDom();

    wireDropzone(dom.dropzone, dom.fileInput, handleIncomingFiles);
    wireDropzone(dom.metaDropzone, dom.metaFileInput, function (files) {
      var csv = files.find(function (f) { return classifyExtension(f.name) === 'csv'; });
      if (csv) handleMetadataCsv(csv);
      else renderMetaPreview('"' + files[0].name + '" is not a .csv file. Drop your sample metadata as a CSV export.');
    });

    dom.metaTabCsv.addEventListener('click', function () { setMetaTab('csv'); });
    dom.metaTabManual.addEventListener('click', function () { setMetaTab('manual'); });
    dom.metaTemplateLink.addEventListener('click', downloadTemplate);
    dom.metaRemoveBtn.addEventListener('click', clearMetadata);
    dom.metaRemoveBtn2.addEventListener('click', clearMetadata);

    dom.metaDate.addEventListener('input', function () { manualMetadata.date = dom.metaDate.value; recomputeAndRender(); });
    dom.metaSource.addEventListener('change', function () { manualMetadata.source = dom.metaSource.value; recomputeAndRender(); });
    dom.metaLocation.addEventListener('input', function () { manualMetadata.location = dom.metaLocation.value; recomputeAndRender(); });

    dom.continueBtn.addEventListener('click', function () {
      if (dom.continueBtn.disabled) return;
      window.CholeraFlowState.update({ maxStepReached: 2 });
      window.location.href = 'screen2_select_mode.html';
    });

    function handleRemoveClick(e) {
      var btn = e.target.closest('[data-action]');
      if (!btn) return;
      var action = btn.getAttribute('data-action');
      var key = btn.getAttribute('data-key');
      if (action === 'remove-rejected') removeRejected(key);
      else if (action === 'remove-sample') removeSampleById(key);
      recomputeAndRender();
    }
    dom.errorBanner.addEventListener('click', handleRemoveClick);
    dom.sampleRows.addEventListener('click', handleRemoveClick);

    recomputeAndRender();
  }

  window.CholeraFlowUpload = { init: init };
})();
