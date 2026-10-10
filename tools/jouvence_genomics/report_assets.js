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
