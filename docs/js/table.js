(function () {
  async function init() {
    const tableEl = document.getElementById('filo-table');
    if (!tableEl || typeof gridjs === 'undefined') return;      // not the browse page
    if (tableEl.dataset.rendered) return;                       // guard double-render
    tableEl.dataset.rendered = '1';

    const errEl = document.getElementById('filo-error');
    let data;
    try {
      const res = await fetch('data.json', { cache: 'no-cache' });
      if (!res.ok) throw new Error('HTTP ' + res.status);
      data = await res.json();
    } catch (e) {
      if (errEl) errEl.textContent = 'Could not load data.json (' + e.message + ').';
      return;
    }

    const m = data.meta || {};
    const statsEl = document.getElementById('filo-stats');
    if (statsEl) {
      const tile = (n, l) => `<div class="filo-stat"><b>${n}</b><span>${l}</span></div>`;
      statsEl.innerHTML =
        tile(m.n_structures ?? '–', 'structures') +
        tile(m.n_entities ?? '–', 'entities') +
        tile((m.species || []).length, 'species') +
        tile(m.n_protein ?? '–', 'protein chains');
    }

    const footEl = document.getElementById('filo-foot');
    if (footEl) {
      footEl.innerHTML = 'Generated ' + (m.generated || '') +
        ' · <a href="data.json">data.json</a>';
    }

    const columns = (data.columns || []).map((name, i) => {
      if (i === 0) {
        return {
          name,
          formatter: (cell) => gridjs.html(
            `<a href="https://www.rcsb.org/structure/${String(cell).toUpperCase()}"` +
            ` target="_blank" rel="noopener">${cell}</a>`),
        };
      }
      return { name };
    });

    new gridjs.Grid({
      columns,
      data: data.rows,
      search: true,
      sort: true,
      resizable: true,
      pagination: { limit: 25, summary: true },
      style: { table: { 'white-space': 'nowrap', 'font-size': '13px' } },
    }).render(tableEl);
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
