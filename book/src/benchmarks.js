// Draw the benchmarks chapter from benchmarks.csv: one SVG chart per dataset
// and a table of every result. Colors come from the book theme (custom.css).
(function () {
  'use strict';

  var PROGRAMS = ['fastp', 'fasterp'];
  var SVG = 'http://www.w3.org/2000/svg';

  function parse(text) {
    var lines = text.trim().split('\n');
    var keys = lines[0].split(',');
    return lines.slice(1).map(function (line) {
      var values = line.split(','), row = {};
      keys.forEach(function (k, i) { row[k] = values[i]; });
      ['workers', 'mean', 'stddev', 'min', 'max', 'output_bytes'].forEach(function (k) { row[k] = +row[k]; });
      return row;
    });
  }

  function unique(rows, key) {
    var seen = [];
    rows.forEach(function (r) { if (seen.indexOf(r[key]) < 0) seen.push(r[key]); });
    return seen;
  }

  function el(name, attrs, parent, text) {
    var node = document.createElementNS(SVG, name);
    Object.keys(attrs).forEach(function (k) { node.setAttribute(k, attrs[k]); });
    if (text !== undefined) node.textContent = text;
    if (parent) parent.appendChild(node);
    return node;
  }

  function niceMax(value) {
    var step = Math.pow(10, Math.floor(Math.log10(value)) - 1);
    var ticks = [1, 2, 2.5, 5, 10].map(function (m) { return m * step; });
    for (var i = 0; i < ticks.length; i++) if (value / ticks[i] <= 5) return { max: Math.ceil(value / ticks[i]) * ticks[i], step: ticks[i] };
  }

  // One panel: worker counts down the side, a fastp and a fasterp bar for each.
  // Drawn at its real pixel width W so text stays at the page's size.
  function panel(rows, workload, scale, W) {
    var workers = unique(rows, 'workers');
    var left = 40, right = 44, top = 8, group = 44, bar = 15, bottom = 24;
    var H = top + group * workers.length + bottom;
    var plot = W - left - right;
    var x = function (v) { return left + (v / scale.max) * plot; };
    var svg = el('svg', { width: W, height: H, viewBox: '0 0 ' + W + ' ' + H, class: 'bench-panel', role: 'img' });
    el('title', {}, svg, workload);

    for (var t = 0; t <= scale.max + 1e-9; t += scale.step) {
      el('line', { x1: x(t), x2: x(t), y1: top, y2: H - bottom, class: 'bench-grid' }, svg);
      el('text', { x: x(t), y: H - 8, class: 'bench-tick', 'text-anchor': 'middle' }, svg, +t.toFixed(2));
    }

    workers.forEach(function (w, g) {
      var y0 = top + g * group + (group - 2 * bar - 2) / 2;
      el('text', { x: left - 8, y: y0 + bar + 1, class: 'bench-tick', 'text-anchor': 'end', 'dominant-baseline': 'middle' }, svg, '-w ' + w);
      var means = {}, edge = 0;
      PROGRAMS.forEach(function (program, p) {
        var r = rows.filter(function (r) { return r.workers === w && r.program === program; })[0];
        if (!r) return;
        var y = y0 + p * (bar + 2);
        means[program] = r.mean;
        edge = Math.max(edge, r.max);
        var g2 = el('g', { class: 'bench-mark' }, svg);
        el('title', {}, g2, program + ' ' + r.version + ', -w ' + w + ': ' + r.mean.toFixed(2) + ' s (min ' + r.min.toFixed(2) + ', max ' + r.max.toFixed(2) + ')');
        el('rect', { x: left, y: y, width: Math.max(1, x(r.mean) - left), height: bar, rx: 2, class: 'bench-' + program }, g2);
        el('line', { x1: x(r.min), x2: x(r.max), y1: y + bar / 2, y2: y + bar / 2, class: 'bench-whisker' }, g2);
        // a wider invisible target so the tooltip is easy to hit
        el('rect', { x: left, y: y, width: plot, height: bar, class: 'bench-hit' }, g2);
      });
      if (means.fastp && means.fasterp) {
        el('text', { x: x(edge) + 6, y: y0 + bar + 1, class: 'bench-label', 'dominant-baseline': 'middle' }, svg,
          (means.fastp / means.fasterp).toFixed(1) + 'x');
      }
    });
    return svg;
  }

  function chart(root, rows) {
    root.textContent = '';
    var gap = 16, workloads = unique(rows, 'workload').length;
    var cols = Math.max(1, Math.min(workloads, Math.floor((root.clientWidth + gap) / (200 + gap))));
    var width = Math.floor((root.clientWidth - gap * (cols - 1)) / cols);
    unique(rows, 'dataset').forEach(function (dataset) {
      var data = rows.filter(function (r) { return r.dataset === dataset; });
      var scale = niceMax(Math.max.apply(null, data.map(function (r) { return r.max; })));
      var section = document.createElement('figure');
      section.className = 'bench-dataset';
      var caption = document.createElement('figcaption');
      caption.textContent = dataset;
      section.appendChild(caption);
      var grid = document.createElement('div');
      grid.className = 'bench-grid-panels';
      grid.style.gridTemplateColumns = 'repeat(' + cols + ', ' + width + 'px)';
      unique(data, 'workload').forEach(function (workload) {
        var cell = document.createElement('div');
        var title = document.createElement('div');
        title.className = 'bench-panel-title';
        title.textContent = workload;
        cell.appendChild(title);
        cell.appendChild(panel(data.filter(function (r) { return r.workload === workload; }), workload, scale, width));
        grid.appendChild(cell);
      });
      section.appendChild(grid);
      root.appendChild(section);
    });
  }

  // One table per dataset so the columns fit the page.
  function table(root, rows) {
    unique(rows, 'dataset').forEach(function (dataset) {
      var heading = document.createElement('p');
      heading.innerHTML = '<strong></strong>';
      heading.firstChild.textContent = dataset;
      root.appendChild(heading);
      var t = document.createElement('table');
      t.innerHTML = '<thead><tr><th>Workload</th><th>Workers</th>' +
        '<th>fastp (s)</th><th>fasterp (s)</th><th>Speedup</th><th>Output size vs fastp</th>' +
        '<th>Identical output</th></tr></thead>';
      var body = document.createElement('tbody');
      rows.filter(function (r) { return r.dataset === dataset && r.program === 'fastp'; }).forEach(function (a) {
        var b = rows.filter(function (r) {
          return r.program === 'fasterp' && r.dataset === a.dataset && r.workload === a.workload && r.workers === a.workers;
        })[0];
        if (!b) return;
        var tr = document.createElement('tr');
        [a.workload, a.workers,
          a.mean.toFixed(2) + ' ± ' + a.stddev.toFixed(2),
          b.mean.toFixed(2) + ' ± ' + b.stddev.toFixed(2),
          (a.mean / b.mean).toFixed(2) + 'x',
          // Gzip output size, recorded since rows gained the column.
          a.output_bytes > 0 && b.output_bytes > 0
            ? (b.output_bytes >= a.output_bytes ? '+' : '') +
              (100 * (b.output_bytes / a.output_bytes - 1)).toFixed(1) + '%'
            : '',
          // Digest of the decompressed records, so equal means byte for byte.
          a.output_sha256 && b.output_sha256
            ? (a.output_sha256 === b.output_sha256 ? 'yes' : 'no')
            : ''].forEach(function (v, i) {
          var td = document.createElement('td');
          td.textContent = v;
          if (i >= 1) td.style.textAlign = 'right';
          tr.appendChild(td);
        });
        body.appendChild(tr);
      });
      t.appendChild(body);
      root.appendChild(t);
    });
  }

  function meta(root, rows) {
    var parts = unique(rows, 'machine').map(function (m) { return m; });
    var versions = PROGRAMS.map(function (p) {
      return p + ' ' + unique(rows.filter(function (r) { return r.program === p; }), 'version').join('/');
    });
    root.textContent = versions.join(', ') + ' on ' + parts.join('; ') + ', measured ' +
      unique(rows, 'date').join(', ') + '. Mean of ' + unique(rows, 'runs').join('/') + ' runs after one warmup.';
  }

  var chartRoot = document.getElementById('bench-chart');
  if (!chartRoot) return;
  var prefix = (typeof path_to_root !== 'undefined' && path_to_root) ? path_to_root : '';
  // Revalidate on every load so new rows show up without a hard reload.
  fetch(prefix + 'benchmarks.csv', { cache: 'no-cache' })
    .then(function (response) { return response.text(); })
    .then(function (text) {
      var rows = parse(text);
      chart(chartRoot, rows);
      var pending;
      window.addEventListener('resize', function () {
        clearTimeout(pending);
        pending = setTimeout(function () { chart(chartRoot, rows); }, 100);
      });
      table(document.getElementById('bench-table'), rows);
      meta(document.getElementById('bench-meta'), rows);
    })
    .catch(function (error) {
      chartRoot.textContent = 'Could not load benchmarks.csv: ' + error;
    });
})();
