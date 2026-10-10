/* Small accessible filters; they change visibility, never medical conclusions. */
document.documentElement.lang='fr';
const frenchLabels={
  'document-mode':'Document','presentation-mode':'Présentation','fullscreen':'Plein écran',
  'previous':'← Précédent','next':'Suivant →'
};
Object.entries(frenchLabels).forEach(([id,label])=>{const el=document.getElementById(id);if(el)el.textContent=label});
document.querySelectorAll('.plot-tools label').forEach(el=>el.textContent='Vue');
document.querySelectorAll('.plot-shell a[download]').forEach(el=>el.textContent='Télécharger les données affichées');
document.querySelectorAll('.plot-shell > details > summary').forEach(el=>el.textContent='Voir les valeurs et télécharger le tableau');
document.querySelectorAll('.plot-canvas').forEach(el=>el.setAttribute('aria-label','Profil génétique et références publiques ; valeurs exactes disponibles dans le tableau'));
// Open every enclosing detail when following a finding or a source citation.
const revealAnchor=()=>{
  const target=document.getElementById(decodeURIComponent(location.hash.slice(1)));
  if(!target)return;
  if(target.tagName==='DETAILS')target.open=true;
  for(let parent=target.parentElement;parent;parent=parent.parentElement){
    if(parent.tagName==='DETAILS')parent.open=true;
  }
  requestAnimationFrame(()=>target.scrollIntoView({block:'start'}));
};
addEventListener('hashchange',revealAnchor);revealAnchor();
const countReferences=()=>{
  const details=document.querySelector('.bibliography-details');
  if(details){
    details.querySelector('summary').textContent=`Consulter les ${details.querySelectorAll('li:not([hidden])').length} références de cette vue`;
    if(document.body.classList.contains('presentation'))details.open=true;
  }
};
document.getElementById('presentation-mode')?.addEventListener('click',()=>{
  document.querySelectorAll('.bibliography-details').forEach(el=>el.open=true);
  countReferences();
});
document.getElementById('document-mode')?.addEventListener('click',countReferences);countReferences();
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
