/* Development-only browser checks for a generated CholeraFlow report.
 *
 * Loads the report from file:// with networking disabled in Chromium, Firefox
 * and/or WebKit (Playwright), walks every page and exercises navigation,
 * search, the details drawer, filters, themes, exports, the mobile layout and
 * print media. Not a runtime dependency of the report or the pipeline.
 *
 *   npm install --prefix /tmp/cf-browser playwright-core@1.56.1
 *   PLAYWRIGHT_BROWSERS_PATH=/tmp/cf-browsers /tmp/cf-browser/node_modules/.bin/playwright-core install firefox webkit
 *   NODE_PATH=/tmp/cf-browser/node_modules PLAYWRIGHT_BROWSERS_PATH=/tmp/cf-browsers \
 *     node tests/browser/report_browser_checks.cjs examples/reports/outbreak_report.html /tmp/cf-shots chromium,firefox,webkit
 *
 * Chromium uses the locally installed Google Chrome (channel "chrome") when
 * available. Pass --shots to save full-page screenshots for visual review.
 */
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const playwright = require('playwright-core');

const [,, reportArg, outArg, browserArg, ...flags] = process.argv;
const report = path.resolve(reportArg || path.join(__dirname, '../../examples/reports/outbreak_report.html'));
const outdir = path.resolve(outArg || path.join(require('node:os').tmpdir(), 'choleraflow-browser-checks'));
const browsers = (browserArg || 'chromium').split(',').filter(Boolean);
const shots = flags.includes('--shots');
fs.mkdirSync(outdir, {recursive: true});

async function launch(name) {
  if (name === 'chromium') {
    try { return await playwright.chromium.launch({channel: 'chrome', headless: true}); }
    catch (error) { return playwright.chromium.launch({headless: true}); }
  }
  return playwright[name].launch({headless: true});
}

async function check(name) {
  const browser = await launch(name);
  const context = await browser.newContext({viewport: {width: 1440, height: 900}, offline: true, reducedMotion: 'reduce', acceptDownloads: true});
  const page = await context.newPage();
  const errors = [];
  const requests = [];
  page.on('pageerror', error => errors.push(String(error && error.message || error)));
  page.on('request', request => { if (!/^(file|data|blob|about):/.test(request.url())) requests.push(request.url()); });
  const started = Date.now();
  await page.goto('file://' + report, {waitUntil: 'load'});
  const loadMs = Date.now() - started;
  assert.equal(await page.evaluate(() => Boolean(window.CF && window.CF.ready)), true, name + ': runtime ready');
  const pages = await page.evaluate(() => Array.from(document.querySelectorAll('main > .cf-page')).map(node => node.id));
  for (const id of pages) {
    await page.evaluate(target => { location.hash = target; }, id);
    await page.waitForTimeout(80);
    const state = await page.evaluate(target => {
      const node = document.getElementById(target);
      return {active: node.classList.contains('is-active'), shown: node.getClientRects().length > 0,
              failed: node.dataset.cfError || null, overflow: document.documentElement.scrollWidth > innerWidth + 1};
    }, id);
    assert.ok(state.active && state.shown, name + ': page ' + id + ' shows');
    assert.equal(state.failed, null, name + ': page ' + id + ' initialised without errors');
    assert.equal(state.overflow, false, name + ': page ' + id + ' does not overflow horizontally');
    if (shots) await page.screenshot({path: path.join(outdir, name + '-' + id + '.png'), fullPage: true});
  }
  const sample = await page.evaluate(() => window.CF.sampleIds[0] || null);
  if (sample) {
    await page.evaluate(() => { location.hash = 'overview'; });
    await page.focus('#cf-search-input');
    await page.keyboard.type(sample);
    await page.keyboard.press('ArrowDown');
    await page.keyboard.press('Enter');
    await page.waitForTimeout(100);
    assert.equal(await page.evaluate(() => document.getElementById('cf-drawer').hidden), false, name + ': drawer opens from search');
    if (shots) await page.screenshot({path: path.join(outdir, name + '-drawer.png')});
    await page.keyboard.press('Escape');
    assert.equal(await page.evaluate(() => document.getElementById('cf-drawer').hidden), true, name + ': Escape closes the drawer');
    const filter = await page.evaluate(() => { const field = window.CF.filterFields()[0]; return field ? [field.key, Object.keys(field.counts)[0]] : null; });
    if (filter) {
      await page.evaluate(([field, value]) => window.CF.setFilter(field, [value]), filter);
      assert.equal(await page.evaluate(() => document.getElementById('cf-filterbar').hidden), false, name + ': filter bar visible');
      await page.evaluate(() => window.CF.clearFilters());
    }
  }
  await page.evaluate(() => window.CF.setTheme('dark'));
  assert.equal(await page.evaluate(() => document.documentElement.dataset.theme), 'dark');
  if (shots) {
    for (const id of ['overview', 'qc', 'phylogeny', 'snp', 'network']) {
      await page.evaluate(target => { location.hash = target; }, id);
      await page.waitForTimeout(80);
      await page.screenshot({path: path.join(outdir, name + '-dark-' + id + '.png'), fullPage: true});
    }
  }
  await page.evaluate(() => window.CF.setTheme('light'));
  const downloads = {};
  const exportsAvailable = await page.evaluate(() => window.CF.exports.list().filter(item => item.id !== 'print' && window.CF.exports.status(item).ok).map(item => item.id));
  for (const id of exportsAvailable) {
    const [download] = await Promise.all([page.waitForEvent('download', {timeout: 5000}), page.evaluate(target => window.CF.exports.run(target), id)]);
    const target = path.join(outdir, name + '-' + download.suggestedFilename());
    await download.saveAs(target);
    downloads[id] = fs.statSync(target).size;
    assert.ok(downloads[id] > 0, name + ': export ' + id + ' is not empty');
  }
  const svgOk = await page.evaluate(() => {
    const figure = Array.from(document.querySelectorAll('.cf-figure')).find(node => node.querySelector('svg:not(.cf-icon)'));
    if (!figure) return null;
    const text = window.CF.buildFigureSvg(figure);
    return !new DOMParser().parseFromString(text, 'image/svg+xml').querySelector('parsererror') && !/<script/i.test(text);
  });
  assert.notEqual(svgOk, false, name + ': standalone SVG export is valid');
  await page.setViewportSize({width: 390, height: 844});
  for (const id of ['overview', 'samples', 'phylogeny', 'amr', 'snp']) {
    await page.evaluate(target => { location.hash = target; }, id);
    await page.waitForTimeout(80);
    assert.equal(await page.evaluate(() => document.documentElement.scrollWidth > innerWidth + 1), false, name + ': mobile ' + id + ' fits the viewport');
    if (shots) await page.screenshot({path: path.join(outdir, name + '-mobile-' + id + '.png'), fullPage: true});
  }
  await page.emulateMedia({media: 'print'});
  assert.equal(await page.evaluate(() => Array.from(document.querySelectorAll('main > .cf-page')).every(node => getComputedStyle(node).display !== 'none')), true, name + ': print shows every page');
  let pdfBytes = null;
  if (name === 'chromium') {
    await page.setViewportSize({width: 1240, height: 1754});
    await page.pdf({path: path.join(outdir, 'report-A4.pdf'), format: 'A4', printBackground: true});
    await page.pdf({path: path.join(outdir, 'report-Letter.pdf'), format: 'Letter', printBackground: true});
    pdfBytes = fs.statSync(path.join(outdir, 'report-A4.pdf')).size;
  }
  assert.deepEqual(errors, [], name + ': no script errors');
  assert.deepEqual(requests, [], name + ': no network requests');
  await browser.close();
  return {browser: name, pages: pages.length, loadMs, exports: downloads, pdfBytes};
}

(async () => {
  const results = [];
  for (const name of browsers) results.push(await check(name));
  console.log(JSON.stringify(results, null, 2));
  console.log('Browser checks passed: ' + results.map(item => item.browser).join(', '));
})().catch(error => { console.error(error); process.exit(1); });
