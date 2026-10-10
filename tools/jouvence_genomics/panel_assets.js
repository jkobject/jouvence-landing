/* Progressive, offline panel controls. Analysis values are never recomputed. */
(function () {
  'use strict';
  const csvColumns = ['chrom', 'pos', 'ref', 'alt', 'rsid', 'gene', 'genotype',
    'effect_allele', 'dosage', 'weight', 'contribution', 'status', 'quality_issues'];

  // Text cells that spreadsheet software could execute are explicitly text.
  // Numeric source values, including negative contributions, retain their type.
  function csvCell(value) {
    if (value === null || value === undefined) return '""';
    let text = typeof value === 'object' ? JSON.stringify(value) : String(value);
    if (typeof value !== 'number' && (/^[\s\uFEFF]*[=+\-@]/u.test(text) || /^[\t\r\n]/.test(text))) text = "'" + text;
    return '"' + text.replace(/"/g, '""') + '"';
  }

  function panelCsv(rows) {
    const extraColumns = [...new Set(rows.flatMap(row => Object.keys(row)))].filter(key => !csvColumns.includes(key));
    const columns = [...csvColumns, ...extraColumns];
    return '\uFEFF' + [columns.map(csvCell).join(','),
      ...rows.map(row => columns.map(key => csvCell(row[key])).join(','))].join('\r\n') + '\r\n';
  }

  // A small testable export is available only to a local Node test runner.
  if (typeof module !== 'undefined' && module.exports) module.exports = {csvCell, panelCsv, csvColumns};
  if (typeof document === 'undefined') return;
  const normalize = text => text.normalize('NFD').replace(/[\u0300-\u036f]/g, '').toLowerCase();
  const displayValue = value => value === null || value === undefined || value === '' || Array.isArray(value) && value.length === 0
    ? '—' : typeof value === 'object' ? JSON.stringify(value) : String(value);

  document.querySelectorAll('[data-panel-dashboard]').forEach(dashboard => {
    if (dashboard.dataset.panelInitialized === 'true') return;
    dashboard.dataset.panelInitialized = 'true';
    dashboard.classList.add('jp-ready');
    const cards = [...dashboard.querySelectorAll('[data-panel-card]')];
    const buttons = [...dashboard.querySelectorAll('[data-panel-category]')];
    const input = dashboard.querySelector('[data-panel-search-input]');
    const count = dashboard.querySelector('[data-panel-count]');
    const includeAudit = dashboard.querySelector('[data-panel-include-audit]');
    let category = 'all';
    const update = () => {
      const query = normalize(input?.value || '');
      let visible = 0;
      cards.forEach(card => {
        const matches = (includeAudit?.checked || card.dataset.panelAuditOnly !== 'true')
          && (category === 'all' || card.dataset.panelCategoryValue === category)
          && normalize(card.dataset.panelSearch || '').includes(query);
        card.hidden = !matches;
        visible += Number(matches);
      });
      dashboard.querySelectorAll('[data-panel-category-group]').forEach(group => {
        group.hidden = ![...group.querySelectorAll('[data-panel-card]')].some(card => !card.hidden);
      });
      if (count) count.textContent = `${visible} panel(s) affiché(s) sur ${cards.length}`;
      const empty = dashboard.querySelector('[data-panel-no-results]');
      if (empty) empty.hidden = visible !== 0;
    };
    buttons.forEach(button => button.addEventListener('click', () => {
      category = button.dataset.panelCategory;
      buttons.forEach(other => other.setAttribute('aria-pressed', String(other === button)));
      update();
    }));
    input?.addEventListener('input', update);
    includeAudit?.addEventListener('change', update);
    update();

    cards.forEach(card => {
      // Parse only on first export/expansion; large panels do not block startup.
      let rows;
      const getRows = () => rows ??= JSON.parse(card.dataset.panelExport);
      const columns = JSON.parse(card.dataset.panelColumns);
      card.querySelector('[data-panel-download]')?.addEventListener('click', () => {
        const blob = new Blob([panelCsv(getRows())], {type: 'text/csv;charset=utf-8'});
        const link = document.createElement('a');
        const url = URL.createObjectURL(blob);
        link.href = url;
        link.download = (card.dataset.panelId.replace(/[^\p{L}\p{N}._-]/gu, '_').slice(0, 100) || 'panel') + '.csv';
        document.body.append(link);
        link.click();
        link.remove();
        setTimeout(() => URL.revokeObjectURL(url), 1000);
      });
      const expand = card.querySelector('[data-panel-expand]');
      expand?.addEventListener('click', () => {
        const table = card.querySelector('tbody');
        const all = getRows();
        const shown = card.querySelector('[data-panel-rows-count]');
        if (expand.getAttribute('aria-expanded') === 'true') {
          while (table.rows.length > 20) table.deleteRow(-1);
          expand.setAttribute('aria-expanded', 'false');
          expand.textContent = `Afficher toutes les lignes (${all.length})`;
          shown.textContent = `20 ligne(s) affichée(s) sur ${all.length}`;
          return;
        }
        let next = 20;
        expand.disabled = true;
        expand.textContent = 'Affichage des lignes…';
        // Build in small batches so a long variant list does not freeze reading.
        const appendBatch = () => {
          const fragment = document.createDocumentFragment();
          const end = Math.min(next + 250, all.length);
          for (; next < end; next++) {
            const row = document.createElement('tr');
            columns.forEach(key => {
              const cell = document.createElement('td');
              cell.dataset.panelCell = key;
              cell.textContent = displayValue(all[next][key]);
              row.append(cell);
            });
            fragment.append(row);
          }
          table.append(fragment);
          shown.textContent = `${next} ligne(s) affichée(s) sur ${all.length}`;
          if (next < all.length) requestAnimationFrame(appendBatch);
          else {
            expand.disabled = false;
            expand.setAttribute('aria-expanded', 'true');
            expand.textContent = 'Revenir aux 20 premières lignes';
          }
        };
        requestAnimationFrame(appendBatch);
      });
    });

    const reducedMotion = window.matchMedia('(prefers-reduced-motion: reduce)');
    if (!reducedMotion.matches && 'IntersectionObserver' in window) {
      const observer = new IntersectionObserver(entries => {
        entries.forEach(entry => {
          if (entry.isIntersecting && !reducedMotion.matches) {
            entry.target.animate?.([{opacity: .7, transform: 'translateY(7px)'}, {opacity: 1, transform: 'translateY(0)'}],
              {duration: 320, easing: 'ease-out'});
            observer.unobserve(entry.target);
          }
        });
      }, {threshold: .08});
      cards.forEach(card => observer.observe(card));
      reducedMotion.addEventListener('change', event => {if (event.matches) observer.disconnect();});
    }
  });
})();
