'use strict';
const plotSpecs = JSON.parse(document.getElementById('plot-data')?.textContent || '{}');
const sections = [...document.querySelectorAll('main > section')];
const toc = document.getElementById('toc');
const progress = document.getElementById('progress');
const fullscreen = document.getElementById('fullscreen');
let presentation = new URLSearchParams(location.search).get('view') === 'presentation';
let current = 0;
let slides = [];

// Reports built with plotly="none" carry a pre-rendered PNG per variant instead of a plotting library.
function renderStaticImage(canvas, variant) {
  if(!/^data:image\/(png|jpeg|webp);base64,/.test(variant.image || '')) throw new Error('Plot view has no embedded image: '+variant.label);
  const image = document.createElement('img');
  image.src = variant.image;
  image.alt = variant.caption || variant.label;
  canvas.replaceChildren(image);
}

function renderPlot(host, variant) {
  const canvas = host.querySelector('.plot-canvas');
  // Python copies of this base layout: report_shell.PLOT_LAYOUT and presentation import_html.BASE_LAYOUT.
  const layout = {paper_bgcolor:'white', plot_bgcolor:'white', font:{family:'Arial, sans-serif',size:14,color:'#172e38'}, margin:{l:72,r:28,t:35,b:65}, autosize:true, hovermode:'closest', ...variant.layout};
  if('image' in variant) renderStaticImage(canvas, variant);
  else Plotly.react(canvas, variant.data, layout, {responsive:true,displaylogo:false,scrollZoom:false,modeBarButtonsToRemove:['lasso2d','select2d','toImage']});
  host.querySelector('.plot-caption').textContent = variant.caption;
  const rows = [['Series',variant.layout.xaxis?.title?.text || 'Category',variant.layout.yaxis?.title?.text || 'Value','Context']];
  variant.data.filter(trace => !trace.meta?.excludeTable).forEach(trace => {
    (trace.x || []).forEach((x,i) => rows.push([trace.name || '',x,trace.y?.[i] ?? '',trace.hovertext?.[i] ?? '']));
  });
  const table = host.querySelector('table');
  table.replaceChildren();
  const head = document.createElement('thead'), body = document.createElement('tbody');
  rows.forEach((row,i) => {const tr=document.createElement('tr'); row.forEach(value=>{const cell=document.createElement(i===0?'th':'td');cell.textContent=String(value);if(i===0)cell.scope='col';tr.append(cell)});(i===0?head:body).append(tr)});
  table.append(head,body);
  const download=host.querySelector('[download]');
  if(download.dataset.blob)URL.revokeObjectURL(download.dataset.blob);
  const csv=rows.map(row=>row.map(value=>'"'+String(value).replaceAll('"','""')+'"').join(',')).join('\n');
  const url=URL.createObjectURL(new Blob([csv],{type:'text/csv;charset=utf-8'}));
  download.href=url;download.dataset.blob=url;
}

document.querySelectorAll('[data-plot]').forEach(host => {
  const variants=plotSpecs[host.dataset.plot], keys=Object.keys(variants);
  if(!keys.length)throw new Error('Empty figure: '+host.dataset.plot);
  host.classList.add('plot-shell');
  const tools=document.createElement('div');tools.className='plot-tools';
  const label=document.createElement('label');label.textContent='View';label.htmlFor='select-'+host.dataset.plot;
  const select=document.createElement('select');select.id=label.htmlFor;
  keys.forEach(key=>{const option=document.createElement('option');option.value=key;option.textContent=variants[key].label;select.append(option)});
  tools.append(label,select);if(keys.length===1)tools.hidden=true;
  const canvas=document.createElement('div');canvas.className='plot-canvas';canvas.setAttribute('role','img');canvas.setAttribute('aria-label',variants[keys[0]].label+'; exact values in the data table below');
  const caption=document.createElement('p');caption.className='plot-caption';caption.setAttribute('aria-live','polite');
  const details=document.createElement('details');details.className='doc-only';
  const summary=document.createElement('summary');summary.textContent='View data and download CSV';
  const download=document.createElement('a');download.download=host.dataset.plot+'.csv';download.textContent='Download displayed data';
  const wrap=document.createElement('div');wrap.className='table-wrap';const table=document.createElement('table');wrap.append(table);details.append(summary,download,wrap);
  host.append(tools,canvas,caption,details);
  select.addEventListener('change',()=>{canvas.setAttribute('aria-label',variants[select.value].label+'; exact values in the data table below');renderPlot(host,variants[select.value])});
  renderPlot(host,variants[keys[0]]);
});

function applyMode(anchor) {
  document.body.classList.toggle('presentation',presentation);
  document.getElementById('document-mode').setAttribute('aria-pressed',String(!presentation));
  document.getElementById('presentation-mode').setAttribute('aria-pressed',String(presentation));
  document.getElementById('slide-controls').hidden=!presentation;
  progress.hidden=!presentation;
  fullscreen.hidden=!presentation;
  const cited=new Set([...document.querySelectorAll('[data-cite]')].filter(el=>presentation ? !el.closest('.prose,.doc-only,[data-document-only]') : !el.closest('[data-slide-only]')).map(el=>el.dataset.cite));
  document.querySelectorAll('[data-reference]').forEach(el=>el.hidden=!cited.has(el.dataset.reference));
  slides=sections.filter(section=>!presentation || (!section.matches('.prose,.doc-only,[data-document-only]') && (section.id!=='bibliography'||cited.size>0)));
  const selected=slides.findIndex(section=>section.id===anchor);
  current=selected>=0?selected:Math.min(current,slides.length-1);
  toc.replaceChildren();
  slides.forEach(section=>{const link=document.createElement('a');link.href='#'+section.id;link.textContent=section.querySelector('h2').textContent;toc.append(link)});
  showSlide(false);
  if(!presentation)spy();
}

function markToc(id) {
  toc.querySelectorAll('a').forEach(link=>link.setAttribute('aria-current',String(link.getAttribute('href')==='#'+id)));
}

function spy() {
  // At the page bottom, short trailing sections can never reach the top; mark the last.
  const bottom=innerHeight+scrollY>=document.documentElement.scrollHeight-2;
  const last=bottom?sections.at(-1):sections.filter(section=>section.getBoundingClientRect().top<=120).at(-1);
  markToc((last||sections[0]).id);
}

function showSlide(scroll=true) {
  sections.forEach(section=>section.hidden=presentation && section!==slides[current]);
  if(presentation && slides[current]){
    document.getElementById('slide-status').textContent=`${current+1} / ${slides.length} · ${slides[current].querySelector('h2').textContent}`;
    if(location.protocol!=='file:')history.replaceState(null,'',`?view=presentation#${slides[current].id}`);
    progress.style.width=`${(current+1)/slides.length*100}%`;
    markToc(slides[current].id);
  }
  document.getElementById('previous').disabled=current===0;
  document.getElementById('next').disabled=current===slides.length-1;
  requestAnimationFrame(()=>document.querySelectorAll('.plot-canvas').forEach(el=>{if(el.getClientRects().length&&typeof Plotly!=='undefined')Plotly.Plots.resize(el)}));
  if(scroll)window.scrollTo({top:0,behavior:'instant'});
}

function enterPresentation() {
  const visible=sections.filter(section=>section.getBoundingClientRect().bottom>120).sort((a,b)=>Math.abs(a.getBoundingClientRect().top-100)-Math.abs(b.getBoundingClientRect().top-100))[0];
  presentation=true;applyMode(visible?.id);window.scrollTo(0,0);
}
function leavePresentation() {
  if(document.fullscreenElement)document.exitFullscreen();
  const anchor=slides[current]?.id;presentation=false;applyMode(anchor);if(location.protocol!=='file:')history.replaceState(null,'',location.pathname+'#'+anchor);document.getElementById(anchor)?.scrollIntoView();
}
function toggleFullscreen() {
  // Embedded views (iframes, previews) often block fullscreen; say so instead of doing nothing.
  const blocked=()=>{document.getElementById('slide-status').textContent='Le plein écran est bloqué ici. Ouvrir le fichier dans un onglet dédié puis appuyer sur F.'};
  if(document.fullscreenElement)document.exitFullscreen();
  else if(document.fullscreenEnabled)document.documentElement.requestFullscreen().catch(blocked);
  else blocked();
}
function step(delta) {
  current=Math.max(0,Math.min(slides.length-1,current+delta));showSlide();
}

document.getElementById('presentation-mode').addEventListener('click',enterPresentation);
document.getElementById('document-mode').addEventListener('click',()=>{if(presentation)leavePresentation()});
document.getElementById('previous').addEventListener('click',()=>step(-1));
document.getElementById('next').addEventListener('click',()=>step(1));
fullscreen.addEventListener('click',toggleFullscreen);
document.addEventListener('fullscreenchange',()=>fullscreen.setAttribute('aria-pressed',String(!!document.fullscreenElement)));
addEventListener('scroll',()=>{if(!presentation)spy()},{passive:true});
document.addEventListener('keydown',event=>{
  if(event.altKey||event.ctrlKey||event.metaKey||event.target.closest('input,select,textarea,[contenteditable]'))return;
  const key=event.key.length===1?event.key.toLowerCase():event.key;
  if(key==='p'){event.preventDefault();presentation?leavePresentation():enterPresentation();return}
  if(!presentation)return;
  if(key==='Escape'){leavePresentation();return}
  if(key==='f'){event.preventDefault();toggleFullscreen();return}
  if(event.target.closest('button,a,summary'))return;
  const moves={ArrowRight:1,PageDown:1,' ':event.shiftKey?-1:1,ArrowLeft:-1,PageUp:-1,Home:-Infinity,End:Infinity};
  if(key in moves){event.preventDefault();step(moves[key])}
});
document.addEventListener('click',event=>{
  const link=event.target.closest('a[href^="#"]');if(!link)return;
  const target=document.getElementById(link.getAttribute('href').slice(1));if(!target)return;
  if(presentation){const index=slides.indexOf(target.closest('section'));if(index>=0){event.preventDefault();current=index;showSlide();target.scrollIntoView({block:'start'})}}
});
// Click-to-advance, except on controls, figures, tables and while text is selected.
document.addEventListener('click',event=>{
  if(!presentation||event.defaultPrevented||!event.target.closest('main')||event.target.closest('a,button,input,select,textarea,label,summary,details,img,figure,table,.plot-shell,.table-wrap')||String(getSelection()))return;
  step(1);
});
applyMode(location.hash.slice(1));

/* Small accessible filters; they change visibility, never medical conclusions. */
document.documentElement.lang='fr';
const frenchLabels={
  'document-mode':'Document','presentation-mode':'Présentation','fullscreen':'Plein écran',
  'previous':'← Précédent','next':'Suivant →'
};
Object.entries(frenchLabels).forEach(([id,label])=>{const el=document.getElementById(id);if(el)el.textContent=label});
document.querySelectorAll('[data-filter-scope]').forEach(scope=>{
  const input=scope.querySelector('[data-gene-search]');
  const select=scope.querySelector('[data-tier-select]');
  const cards=[...scope.querySelectorAll('[data-gene-card]')];
  const count=scope.querySelector('[data-filter-count]');
  const normalize=text=>text.normalize('NFD').replace(/[\u0300-\u036f]/g,'').toLowerCase();
  const update=()=>{
    const query=normalize(input?.value||'');const tier=select?.value||'all';
    let visible=0;
    cards.forEach(card=>{
      const matches=normalize(card.dataset.search||card.textContent).includes(query)&&
        (tier==='all'||(card.dataset.tiers||'').split(' ').includes(tier));
      card.hidden=!matches;visible+=matches?1:0;
    });
    if(count)count.textContent=`${visible} fiche(s) affichée(s) sur ${cards.length}`;
  };
  input?.addEventListener('input',update);select?.addEventListener('change',update);
  scope.querySelector('[data-reset-filters]')?.addEventListener('click',()=>{
    if(input)input.value='';if(select)select.value='all';update();
  });
  scope.querySelector('[data-open-cards]')?.addEventListener('click',()=>cards.forEach(card=>{if(!card.hidden)card.open=true}));
  scope.querySelector('[data-close-cards]')?.addEventListener('click',()=>cards.forEach(card=>card.open=false));
  update();
});
addEventListener('beforeprint',()=>document.querySelectorAll('details').forEach(el=>{
  el.dataset.openBeforePrint=String(el.open);el.open=true;
}));
addEventListener('afterprint',()=>document.querySelectorAll('details[data-open-before-print]').forEach(el=>{
  el.open=el.dataset.openBeforePrint==='true';delete el.dataset.openBeforePrint;
}));
