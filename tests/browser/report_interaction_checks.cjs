/* Additional file:// QA: viewport coverage, numeric cell inspection, figure
 * downloads, static fallback, and cohort-size benchmarks. Development only.
 * NODE_PATH=/tmp/cf-browser/node_modules node tests/browser/report_interaction_checks.cjs REPORT OUTDIR [SCALE_DIR]
 */
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const {chromium} = require('playwright-core');
const [,, reportArg, outArg, scaleArg] = process.argv;
const report=path.resolve(reportArg), out=path.resolve(outArg);
fs.mkdirSync(out,{recursive:true});
(async()=>{
 const browser=await chromium.launch({channel:'chrome'});
 const context=await browser.newContext({offline:true,acceptDownloads:true,reducedMotion:'reduce'});
 const page=await context.newPage();const errors=[],requests=[];
 page.on('pageerror',e=>errors.push(e.message));page.on('request',r=>{if(/^https?:/.test(r.url()))requests.push(r.url());});
 await page.goto('file://'+report);await page.evaluate(()=>window.CF.initAll());
 const pages=await page.evaluate(()=>Array.from(document.querySelectorAll('main>.cf-page')).map(n=>n.id));
 for(const [width,height] of [[1920,1080],[1440,900],[1280,800],[1024,768],[768,1024],[390,844]]){
  await page.setViewportSize({width,height});
  for(const id of pages){
   await page.evaluate(id=>window.CF.go(id),id);
   await page.waitForFunction(id=>window.CF.state.page===id,id);
   await page.evaluate(()=>new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve))));
   assert.equal(await page.evaluate(()=>document.documentElement.scrollWidth>innerWidth+1),false,`${width}x${height}: ${id} fits`);
  }
 }
 await page.setViewportSize({width:1440,height:900});
 await page.evaluate(()=>window.CF.go('phylogeny'));
 await page.waitForFunction(()=>window.CF.state.page==='phylogeny');
 const layout=page.locator('[aria-label="Tree layout"]');
 if(await layout.count()){
  await layout.selectOption('circular');
  assert.equal(await page.locator('#cohort-tree-circular').isVisible(),true);
  await page.screenshot({path:path.join(out,'chromium-circular-tree.png'),fullPage:true});
  await page.evaluate(()=>window.CF.setTheme('dark'));
  await page.screenshot({path:path.join(out,'chromium-dark-circular-tree.png'),fullPage:true});
  await page.evaluate(()=>window.CF.setTheme('light'));
 }
 await page.evaluate(()=>window.CF.go('snp'));
 await page.waitForFunction(()=>window.CF.state.page==='snp');
 const cell=page.locator('#snp .cf-distance-pair').first();
 if(await cell.count()){
  await cell.click();
  assert.equal(await page.locator('#cf-drawer').isVisible(),true);
  assert.match(await page.locator('#cf-drawer').innerText(),/transmission/);
  await page.screenshot({path:path.join(out,'chromium-pair-inspection.png')});
  await page.keyboard.press('Escape');
 }
 // PNG and SVG exports must produce usable local files, not decorative buttons.
 await page.evaluate(()=>window.CF.go('workflow'));
 await page.waitForFunction(()=>window.CF.state.page==='workflow');
 const fig=page.locator('#pipeline-schematic');
 const [png]=await Promise.all([page.waitForEvent('download'),page.evaluate(()=>window.CF.exportPng(document.getElementById('pipeline-schematic')))]);
 const pngPath=path.join(out,png.suggestedFilename());await png.saveAs(pngPath);
 assert.equal(fs.readFileSync(pngPath).subarray(0,8).toString('hex'),'89504e470d0a1a0a');
 const nojs=await browser.newContext({offline:true,javaScriptEnabled:false,viewport:{width:1440,height:900}});
 const staticPage=await nojs.newPage();await staticPage.goto('file://'+report);
 assert.equal(await staticPage.locator('main>.cf-page:visible').count(),25);
 assert.ok(await staticPage.locator('svg.tree-svg:visible').count()>=2);
 assert.match(await staticPage.locator('body').innerText(),/Synthetic demonstration/);
 await nojs.close();
 const benchmarks=[];
 if(scaleArg){
  for(const n of [0,1,2,3,10,100,500]){
   const file=path.resolve(scaleArg,`cohort-${n}.html`);if(!fs.existsSync(file))continue;
   const start=Date.now();await page.goto('file://'+file,{waitUntil:'load'});const loadMs=Date.now()-start;
   const measure=await page.evaluate(()=>{
    const cf=window.CF;const field=cf.filterFields()[0];const start=performance.now();
    if(field)cf.setFilter(field.key,[Object.keys(field.counts)[0]]);
    const filterMs=performance.now()-start;cf.clearFilters();
    const t=performance.now();cf.initAll();const initMs=performance.now()-t;
    return {samples:cf.sampleIds.length,filterMs,initMs,failed:Array.from(document.querySelectorAll('[data-cf-error]')).map(n=>n.dataset.cfError)};
   });assert.equal(measure.samples,n);assert.deepEqual(measure.failed,[]);benchmarks.push({n,bytes:fs.statSync(file).size,loadMs,...measure});
  }
 }
 assert.deepEqual(errors,[]);assert.deepEqual(requests,[]);
 fs.writeFileSync(path.join(out,'interaction-results.json'),JSON.stringify({viewports:6,pages:pages.length,offlineRequests:requests,benchmarks},null,2));
 await context.close();await browser.close();console.log(JSON.stringify({viewports:6,pages:pages.length,benchmarks},null,2));
})().catch(e=>{console.error(e);process.exit(1);});
