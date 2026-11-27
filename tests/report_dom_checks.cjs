/* Development-only regression checks for the committed synthetic report.
 * Run with Node 22+ and jsdom on NODE_PATH (neither is a report runtime dependency):
 *   NODE_PATH=/path/to/node_modules node tests/report_dom_checks.cjs [report.html]
 * Browser-level checks (layout, Firefox/WebKit) live in tests/browser/.
 */
const fs = require('node:fs');
const path = require('node:path');
const assert = require('node:assert/strict');
const {JSDOM, VirtualConsole} = require('jsdom');

const file = process.argv[2] || path.join(__dirname, '../examples/reports/outbreak_report.html');
const html = fs.readFileSync(file, 'utf8');
const errors = [];
const consoleErrors = [];
const virtualConsole = new VirtualConsole();
virtualConsole.on('jsdomError', error => errors.push(error.message));
virtualConsole.on('error', message => consoleErrors.push(String(message)));
let exported = null;
let exportedName = null;
const dom = new JSDOM(html, {runScripts: 'dangerously', pretendToBeVisual: true, virtualConsole, beforeParse(window) {
  window.Blob = class { constructor(parts, options) { this.parts = parts; this.type = options && options.type; } };
  window.URL.createObjectURL = blob => { exported = blob.parts.map(String).join(''); return 'blob:local-test'; };
  window.URL.revokeObjectURL = () => {};
  window.HTMLAnchorElement.prototype.click = function () { exportedName = this.getAttribute('download'); };
  window.print = () => {};
  window.scrollTo = () => {};
  window.Element.prototype.scrollIntoView = function () {};
  window.matchMedia = () => ({matches: false, addEventListener() {}, removeEventListener() {}});
}});
const {window} = dom;
const doc = window.document;
const CF = window.CF;
const fire = (node, type, init) => node.dispatchEvent(new window.Event(type, Object.assign({bubbles: true}, init || {})));
const key = (node, keyName) => node.dispatchEvent(new window.KeyboardEvent('keydown', {key: keyName, bubbles: true}));
const go = hash => { window.location.hash = hash; fire(window, 'hashchange'); };
const activePages = () => Array.from(doc.querySelectorAll('main > .cf-page.is-active')).map(page => page.id);

/* ── Boot, structure and offline guarantees ─────────────────────────────── */
assert.deepEqual(errors, [], 'no uncaught script errors');
assert.ok(CF && CF.ready, 'runtime booted');
const PAGES = ['overview', 'workflow', 'run', 'samples', 'metadata', 'qc', 'typing', 'virulence', 'amr',
  'classifications', 'phylogeny', 'snp', 'cgmlst', 'network', 'global', 'clusters', 'interpretation',
  'limitations', 'follow-up', 'methods', 'software', 'provenance', 'downloads', 'guide', 'about'];
assert.deepEqual(Array.from(doc.querySelectorAll('main > .cf-page')).map(page => page.id), PAGES);
PAGES.forEach(id => assert.ok(doc.querySelector('.cf-nav-link[href="#' + id + '"]'), 'sidebar link for ' + id));
assert.equal(doc.querySelectorAll('script[src], link[rel=stylesheet], img[src^=http], iframe, object, embed').length, 0);
assert.equal(doc.querySelectorAll('script:not([type])').length, 1, 'one executable inline script');
assert.equal(doc.querySelectorAll('script[type="application/json"]').length, 1, 'one inert data island');
assert.match(doc.querySelector('meta[http-equiv="Content-Security-Policy"]').content, /connect-src 'none'/);
assert.ok(doc.body.textContent.includes('Synthetic demonstration — not validated surveillance results.'));
assert.equal(CF.sampleIds.length, CF.data.sample_count);
assert.ok(doc.documentElement.classList.contains('js'));

/* ── Paged navigation, deep links and legacy anchors ───────────────────── */
assert.deepEqual(activePages(), ['overview']);
go('qc');
assert.deepEqual(activePages(), ['qc']);
assert.equal(doc.querySelector('.cf-nav-link[aria-current]').dataset.page, 'qc');
assert.equal(doc.querySelector('[data-cf-crumb-page]').textContent, 'Quality control');
go('s6');
assert.deepEqual(activePages(), ['phylogeny'], 'legacy #s6 anchor opens the phylogeny page');
go('overview-review');
assert.deepEqual(activePages(), ['overview'], 'element anchors open their page');
key(doc.body, ']');
assert.equal(window.location.hash, '#workflow', '] moves to the next page');

/* ── Search combobox → sample drawer → Escape restores focus ───────────── */
const input = doc.getElementById('cf-search-input');
input.focus();
input.value = 'KENYA';
fire(input, 'input');
const options = Array.from(doc.querySelectorAll('#cf-search-results [role=option]'));
assert.ok(options.length >= 1);
assert.equal(input.getAttribute('aria-expanded'), 'true');
key(input, 'ArrowDown');
assert.ok(input.getAttribute('aria-activedescendant'));
key(input, 'Enter');
const drawer = doc.getElementById('cf-drawer');
assert.equal(drawer.hidden, false, 'drawer opens from search');
assert.equal(doc.querySelector('[data-cf-drawer-title]').textContent, 'VC_KENYA_003');
assert.equal(CF.state.selected, 'VC_KENYA_003');
assert.ok(doc.querySelectorAll('.is-selected[data-cf-sample="VC_KENYA_003"]').length >= 1, 'selection is highlighted in linked views');
key(doc.body, 'Escape');
assert.equal(drawer.hidden, true, 'Escape closes the drawer');

/* ── Global filters are display-only and visible ───────────────────────── */
CF.setFilter('country', ['Kenya']);
assert.equal(doc.getElementById('cf-filterbar').hidden, false);
assert.equal(doc.querySelectorAll('[data-cf-filter-chips] .cf-chip').length, 1);
assert.deepEqual(Array.from(CF.filteredIds()), ['VC_KENYA_003']);
const filterableScope = doc.querySelector('[data-scope="filterable"] [data-cf-scope] span');
if (filterableScope) assert.match(filterableScope.textContent, /Filtered subset · 1 of 3/);
const fixedScope = doc.querySelector('[data-scope="fixed"] [data-cf-scope] span');
if (fixedScope) assert.match(fixedScope.textContent, /Full cohort · not filtered/);
const metricsBefore = JSON.stringify(CF.data.metrics);
doc.querySelector('[data-cf-action="clear-filters"]').click();
assert.equal(CF.filteredIds(), null);
assert.equal(doc.getElementById('cf-filterbar').hidden, true);
assert.equal(JSON.stringify(CF.data.metrics), metricsBefore, 'filters never change reported values');

/* ── Theme, layout and presentation mode ────────────────────────────────── */
doc.querySelector('[data-cf-theme="dark"]').click();
assert.equal(doc.documentElement.dataset.theme, 'dark');
doc.querySelector('[data-cf-theme="light"]').click();
assert.equal(doc.documentElement.dataset.theme, 'light');
doc.querySelector('[data-cf-layout="continuous"]').click();
assert.equal(activePages().length, PAGES.length, 'continuous layout shows every page');
doc.querySelector('[data-cf-layout="paged"]').click();
doc.querySelector('[data-cf-action="presentation"]').click();
assert.ok(doc.documentElement.classList.contains('cf-present'));
key(doc.body, 'Escape');
assert.ok(!doc.documentElement.classList.contains('cf-present'));

/* ── Exports: registry, CSV formula protection, standalone SVG ──────────── */
const exportIds = CF.exports.list().map(item => item.id);
['print', 'samples-csv', 'samples-json', 'newick-cohort', 'matrix-snp', 'matrix-cgmlst'].forEach(id => assert.ok(exportIds.includes(id), 'export ' + id));
CF.exports.run('samples-csv');
const csvRows = exported.replace(/^﻿/, '').trim().split('\r\n');
assert.equal(csvRows.length, CF.sampleIds.length + 1);
assert.match(csvRows[0], /^"sample_id"/);
assert.match(exportedName, /^choleraflow-sample-records-.*\.csv$/);
assert.equal(CF.csv([['=HYPERLINK("untrusted")', '+1', '@x', 'plain']]).replace('﻿', '').trim(),
  '"\'=HYPERLINK(""untrusted"")","\'+1","\'@x","plain"');
CF.exports.run('matrix-snp');
assert.ok(exported.startsWith('\tVC_HAITI_001\tVC_KENYA_003\tVC_MALAWI_002\n'), 'matrix export keeps supplied labels and order');
CF.exports.run('newick-cohort');
assert.equal(exported.trim(), '(VC_HAITI_001:0.001,VC_KENYA_003:0.001,VC_MALAWI_002:0.02);');
const figure = Array.from(doc.querySelectorAll('.cf-figure')).find(node => node.querySelector('svg:not(.cf-icon)'));
assert.ok(figure, 'at least one vector figure');
const svgText = CF.buildFigureSvg(figure);
const parsed = new window.DOMParser().parseFromString(svgText, 'image/svg+xml');
assert.equal(parsed.querySelector('parsererror'), null, 'exported SVG is well-formed XML');
assert.equal(parsed.querySelector('script'), null);
assert.ok(parsed.documentElement.querySelector('style').textContent.length > 0, 'export carries its styles');

/* ── Print: closed details open for printing (except verbatim source text) ─ */
const closed = Array.from(doc.querySelectorAll('main details:not([open])'));
const plain = closed.find(node => !node.querySelector('pre'));
const verbatim = closed.find(node => node.querySelector('pre'));
fire(window, 'beforeprint', {bubbles: false});
if (plain) assert.equal(plain.open, true);
if (verbatim) assert.equal(verbatim.open, false);
fire(window, 'afterprint', {bubbles: false});
if (plain) assert.equal(plain.open, false);

/* Presentation-layer controls must inspect unchanged supplied evidence. */
CF.initAll();
const sampleTable = doc.querySelector('[data-cf-table="sample-table"]');
assert.equal(sampleTable.tBodies[0].rows.length, 3);
const columnToggle = doc.querySelector('[data-sx-controls] input[type="checkbox"]');
const columnWasVisible = !sampleTable.tHead.rows[0].cells[1].hidden;
columnToggle.checked = !columnToggle.checked;
fire(columnToggle, 'change');
assert.equal(sampleTable.tHead.rows[0].cells[1].hidden, columnWasVisible, 'column chooser changes visibility');
columnToggle.checked = !columnToggle.checked;
fire(columnToggle, 'change');
const recordSearch = doc.querySelector('[data-sx-rec-filter]');
recordSearch.value = 'VC_HAITI_001'; fire(recordSearch, 'input');
assert.equal(Array.from(doc.querySelectorAll('.sx-rec')).filter(node=>!node.hidden).length, 1);
recordSearch.value = ''; fire(recordSearch, 'input');

const matrixButton = doc.querySelector('#snp .cf-matrix-cell');
assert.ok(matrixButton, 'numeric matrix cells can be inspected');
matrixButton.click();
assert.equal(CF.state.pair.value, matrixButton.textContent);
assert.match(doc.querySelector('[data-cf-drawer-body]').textContent, /does not establish direct transmission/);
CF.select('VC_HAITI_001');
assert.equal(doc.querySelectorAll('#cf-drawer [role="tab"]').length, 6);
assert.match(doc.querySelector('[data-cf-download-centre]').textContent, /SNP distance matrix/);
const treeLayout = doc.querySelector('[aria-label="Tree layout"]');
treeLayout.value = 'circular';
fire(treeLayout, 'change');
assert.ok(doc.getElementById('cohort-tree').classList.contains('cf-tree-inactive'));
assert.ok(!doc.getElementById('cohort-tree-circular').classList.contains('cf-tree-inactive'));
assert.ok(CF.exports.list().some(item=>item.id === 'workflow-svg'));
assert.ok(CF.exports.list().some(item=>item.id === 'network-edges'));

/* ── Without JavaScript every page and figure is present and readable ──── */
const staticDom = new JSDOM(html);
const staticDoc = staticDom.window.document;
assert.equal(staticDoc.querySelectorAll('main > .cf-page').length, PAGES.length);
assert.equal(staticDoc.querySelectorAll('main > .cf-page[hidden]').length, 0);
assert.ok(staticDoc.querySelectorAll('.cf-figure svg:not(.cf-icon), svg.report-chart, svg.tree-svg, svg.report-network').length >= 10, 'figures are server-rendered');
assert.ok(staticDoc.querySelector('.limitations').textContent.includes('Genomic clustering does not prove transmission'));
staticDom.window.close();
dom.window.close();
process.stdout.write('Report DOM checks passed: ' + PAGES.length + ' pages, navigation, deep links, search, drawer, filters, themes, layouts, exports, SVG, print and no-JavaScript content.\n');
