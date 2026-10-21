// render.js - the plot: one column per window (the candidate's region, each junction partner's), one row per
// member (child, father, mother), coverage over reads, in samplot's picture (pairs and splits as lines at
// their insert size, coloured by event type) or as a pileup (rows of alignments with their clips). Reads
// whose ends fall in two windows are drawn across the gutter. The child's informative reads (dk) carry a
// halo. Colour by event, by MAPQ, by informative-or-not, or by taxon. Hit-testing for the tooltip.

export const EVENT_COLORS = {
  'Deletion/Normal': '#1b1b1f', Duplication: '#d62728', Inversion: '#1f77b4', InterChrm: '#2ca02c',
};
export const DK_COLOR = '#e69f00';
export const TAXON_COLORS = { human: '#8a8a8a', nonhuman: '#cc79a7', unclassified: '#bbbbbb' };
export const REGION_COLORS = { own: 'rgba(0,114,178,0.13)', partner: 'rgba(213,94,0,0.13)' };
const ROLE_LABELS = { child: 'child', father: 'father', mother: 'mother' };

function hashJitter(name) {
  let h = 2166136261;
  for (let i = 0; i < name.length; i++) { h ^= name.charCodeAt(i); h = Math.imul(h, 16777619); }
  return ((h >>> 0) % 1000) / 1000 - 0.5;   // [-0.5, 0.5)
}

export function mapqColor(mapq) {
  const q = Math.max(0, Math.min(60, mapq ?? 0)) / 60;
  const v = Math.round(205 - q * 180);
  return `rgb(${v},${v},${v + Math.round(q * 20)})`;
}

export function formatPos(n) { return n.toLocaleString('en-US'); }
export function formatBp(n) {
  if (n >= 1e6) return (n / 1e6).toFixed(2) + ' Mb';
  if (n >= 1e3) return (n / 1e3).toFixed(1) + ' kb';
  return n + ' bp';
}

export class Plot {
  constructor(canvas, options = {}) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.options = Object.assign({
      left: 92, right: 14, top: 34, bottom: 30, gutter: 28, labelH: 16, coverageH: 48, readsH: 170, rowGap: 12,
      minPanelW: 120, pileRow: 7, maxPileRows: 48, font: '12px -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif',
      ink: '#1b1b1f', muted: '#6b6e76', line: '#d6d8de', background: '#ffffff', covHigh: 'rgba(90,90,96,0.55)',
      covLow: 'rgba(150,150,156,0.35)',
    }, options);
    this.hits = [];
    this.layoutInfo = null;
  }

  setTheme(theme) { Object.assign(this.options, theme); }

  /** Geometry of the panels for `windows` and `members` (an array of roles) at the canvas' CSS width. */
  layout(windows, members, state) {
    const o = this.options;
    const width = this.canvas.clientWidth || 900;
    const n = windows.length;
    const panelsW = width - o.left - o.right - o.gutter * (n - 1);
    const total = windows.reduce((s, w) => s + (w.end - w.start), 0);
    let widths = windows.map((w) => Math.max(o.minPanelW, panelsW * (w.end - w.start) / total));
    const sum = widths.reduce((a, b) => a + b, 0);
    widths = widths.map((x) => x * panelsW / sum);
    const panels = [];
    let x = o.left;
    for (let i = 0; i < n; i++) {
      panels.push({ x, w: widths[i], window: windows[i] });
      x += widths[i] + o.gutter;
    }
    const readsH = state.mode === 'pileup'
      ? Math.max(60, Math.min(o.maxPileRows, state.pileRows || 20) * o.pileRow + 10)
      : o.readsH;
    const rows = [];
    let y = o.top;
    for (const role of members) {
      rows.push({ role, y, labelY: y + 12, covY: y + o.labelH, covH: o.coverageH, readsY: y + o.labelH + o.coverageH + 4, readsH });
      y += o.labelH + o.coverageH + 4 + readsH + o.rowGap;
    }
    const height = y - o.rowGap + o.bottom;
    this.layoutInfo = { width, height, panels, rows, readsH };
    return this.layoutInfo;
  }

  xOf(panel, pos) {
    const w = panel.window;
    const f = (pos - w.start) / (w.end - w.start);
    return panel.x + Math.max(-0.02, Math.min(1.02, f)) * panel.w;
  }

  /**
   * Draw everything.
   * @param data { windows, members: { role: memberData }, regions: [{chrom,start,end,own,...}], breakpoints }
   * @param state { mode: 'samplot'|'pileup', colorBy: 'event'|'mapq'|'informative'|'taxon', hideProper, minMapq,
   *                highlightDk, hover }
   */
  draw(data, state) {
    const o = this.options;
    const roles = Object.keys(data.members);
    const L = this.layout(data.windows, roles, state);
    const dpr = window.devicePixelRatio || 1;
    this.canvas.width = Math.round(L.width * dpr);
    this.canvas.height = Math.round(L.height * dpr);
    this.canvas.style.height = L.height + 'px';
    const ctx = this.ctx;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.fillStyle = o.background;
    ctx.fillRect(0, 0, L.width, L.height);
    ctx.font = o.font;
    this.hits = [];

    // global scales: coverage and insert size shared by the members, so the parents read against the child
    let maxCov = 1;
    const inserts = [];
    for (const role of roles) {
      const m = data.members[role];
      if (!m) continue;
      maxCov = Math.max(maxCov, m.maxCoverage);
      for (const p of m.pairs) if (p.insert && !p.cross) inserts.push(p.insert);
      for (const s of m.splits) if (s.insert && !s.cross) inserts.push(s.insert);
    }
    inserts.sort((a, b) => a - b);
    let maxInsert = inserts.length ? inserts[Math.min(inserts.length - 1, Math.floor(inserts.length * 0.995))] : 1000;
    maxInsert = Math.max(maxInsert, 600);
    const scales = { maxCov, maxInsert };

    // window titles and axes
    for (const panel of L.panels) this._drawWindowHeader(panel, data);
    for (const row of L.rows) {
      const m = data.members[row.role];
      this._drawRowLabel(row, m, data);
      for (let wi = 0; wi < L.panels.length; wi++) {
        const panel = L.panels[wi];
        this._drawRegions(panel, row, data.regions, data.breakpoints, wi);
        if (m) this._drawCoverage(panel, row, m, wi, scales, state);
        this._drawAxis(panel, row, wi === 0, scales, state);
      }
      if (!m) { this._drawMissing(row, L); continue; }
      if (state.mode === 'pileup') this._drawPileup(row, m, L, state, data);
      else this._drawSamplot(row, m, L, state, scales, data);
      if (m) this._drawCaptions(row, m, L);
    }
    this._drawXAxis(L);
  }

  // ---- pieces ----------------------------------------------------------------------------------------

  _drawWindowHeader(panel, data) {
    const ctx = this.ctx, o = this.options, w = panel.window;
    ctx.fillStyle = o.ink;
    ctx.textAlign = 'left';
    ctx.textBaseline = 'alphabetic';
    ctx.font = `bold ${o.font}`;
    const title = `${w.chrom}:${formatPos(w.start + 1)}-${formatPos(w.end)}`;
    ctx.fillText(title, panel.x, 14);
    ctx.font = o.font;
    ctx.fillStyle = o.muted;
    ctx.fillText(formatBp(w.end - w.start), panel.x + ctx.measureText(title).width + 8, 14);
  }

  _drawRowLabel(row, m, data) {
    const ctx = this.ctx, o = this.options;
    ctx.fillStyle = o.ink;
    ctx.font = `bold ${o.font}`;
    ctx.textAlign = 'left';
    ctx.fillText(ROLE_LABELS[row.role] || row.role, 6, row.labelY);
    ctx.font = o.font;
    ctx.fillStyle = o.muted;
    const sample = data.samples && data.samples[row.role];
    if (sample && sample !== row.role) ctx.fillText(sample, 6, row.labelY + 13);
  }

  _drawMissing(row, L) {
    const ctx = this.ctx, o = this.options;
    ctx.fillStyle = o.muted;
    ctx.textAlign = 'left';
    ctx.fillText('no reads for this member', o.left + 6, row.readsY + 16);
  }

  _drawRegions(panel, row, regions, breakpoints, wi) {
    const ctx = this.ctx, o = this.options, w = panel.window;
    const top = row.covY, bottom = row.readsY + row.readsH;
    for (const r of regions || []) {
      if (r.chrom !== w.chrom || r.end <= w.start || r.start >= w.end) continue;
      const x1 = this.xOf(panel, Math.max(r.start, w.start)), x2 = this.xOf(panel, Math.min(r.end, w.end));
      ctx.fillStyle = r.own ? REGION_COLORS.own : REGION_COLORS.partner;
      ctx.fillRect(x1, top, Math.max(2, x2 - x1), bottom - top);
      if (row.role === 'child') {
        // the region's label at the foot of the reads area, where nothing else is drawn
        ctx.fillStyle = r.own ? '#0072B2' : '#D55E00';
        ctx.font = `10px ${o.font.split('px ')[1]}`;
        ctx.textAlign = 'left';
        const label = `${r.own ? 'candidate' : 'partner'}${r.sv_type ? ' ' + r.sv_type : ''}${r.class ? ' (' + r.class + ')' : ''}`;
        ctx.fillText(label, Math.max(panel.x, Math.min(x1, panel.x + panel.w - ctx.measureText(label).width)), bottom - 3);
        ctx.font = o.font;
      }
    }
    // the breakpoints the clipped reads agree on, as thin marks
    for (const b of (breakpoints || [])) {
      if (b.role !== row.role || b.wi !== wi) continue;
      const x = this.xOf(panel, b.pos);
      ctx.strokeStyle = 'rgba(213,94,0,0.8)';
      ctx.setLineDash([3, 3]);
      ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(x, top); ctx.lineTo(x, bottom); ctx.stroke();
      ctx.setLineDash([]);
    }
  }

  _drawCoverage(panel, row, m, wi, scales, state) {
    const ctx = this.ctx, o = this.options;
    const cov = m.coverage[wi];
    if (!cov) return;
    const n = cov.high.length, y0 = row.covY + row.covH, h = row.covH - 2;
    const step = Math.max(1, Math.floor(n / Math.max(200, panel.w * 2)));
    const px = (i) => panel.x + (i / n) * panel.w;
    const py = (d) => y0 - Math.min(1, d / scales.maxCov) * h;
    // high MAPQ
    ctx.beginPath(); ctx.moveTo(panel.x, y0);
    for (let i = 0; i < n; i += step) ctx.lineTo(px(i), py(cov.high[i]));
    ctx.lineTo(panel.x + panel.w, y0); ctx.closePath();
    ctx.fillStyle = o.covHigh; ctx.fill();
    // low MAPQ, stacked
    ctx.beginPath(); ctx.moveTo(panel.x, y0);
    for (let i = 0; i < n; i += step) ctx.lineTo(px(i), py(cov.high[i] + cov.low[i]));
    for (let i = n - 1; i >= 0; i -= step) ctx.lineTo(px(i), py(cov.high[i]));
    ctx.closePath(); ctx.fillStyle = o.covLow; ctx.fill();
    // informative reads' coverage: the k-mer evidence layer
    if (state.highlightDk !== false) {
      let any = false;
      for (let i = 0; i < n; i += step) if (cov.dk[i] > 0) { any = true; break; }
      if (any) {
        ctx.beginPath(); ctx.moveTo(panel.x, y0);
        for (let i = 0; i < n; i += step) ctx.lineTo(px(i), py(cov.dk[i]));
        ctx.lineTo(panel.x + panel.w, y0); ctx.closePath();
        ctx.fillStyle = 'rgba(230,159,0,0.55)'; ctx.fill();
      }
    }
    ctx.strokeStyle = o.line; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(panel.x, y0 + 0.5); ctx.lineTo(panel.x + panel.w, y0 + 0.5); ctx.stroke();
  }

  _drawAxis(panel, row, first, scales, state) {
    const ctx = this.ctx, o = this.options;
    if (!first) return;
    ctx.fillStyle = o.muted;
    ctx.textAlign = 'right';
    ctx.font = `10px ${o.font.split('px ')[1]}`;
    ctx.fillText(String(Math.round(scales.maxCov)) + '×', o.left - 4, row.covY + 9);
    ctx.fillText('0', o.left - 4, row.covY + row.covH - 1);
    if (state.mode !== 'pileup') {
      ctx.fillText(formatBp(scales.maxInsert), o.left - 4, row.readsY + 9);
      ctx.fillText('insert', o.left - 4, row.readsY + row.readsH / 2);
      ctx.fillText('0', o.left - 4, row.readsY + row.readsH - 8);
    } else {
      ctx.fillText('reads', o.left - 4, row.readsY + 9);
    }
    ctx.font = o.font;
  }

  _drawXAxis(L) {
    const ctx = this.ctx, o = this.options;
    const y = L.height - o.bottom + 6;
    ctx.strokeStyle = o.line; ctx.fillStyle = o.muted;
    ctx.textAlign = 'center'; ctx.font = `10px ${o.font.split('px ')[1]}`;
    for (const panel of L.panels) {
      const w = panel.window, span = w.end - w.start;
      const nice = [50, 100, 200, 250, 500, 1000, 2000, 2500, 5000, 10000, 20000, 25000, 50000, 100000, 250000, 500000, 1000000];
      let tick = nice.find((t) => span / t <= 6) || span;
      ctx.beginPath(); ctx.moveTo(panel.x, y + 0.5); ctx.lineTo(panel.x + panel.w, y + 0.5); ctx.stroke();
      for (let p = Math.ceil(w.start / tick) * tick; p <= w.end; p += tick) {
        const x = this.xOf(panel, p);
        ctx.beginPath(); ctx.moveTo(x, y); ctx.lineTo(x, y + 4); ctx.stroke();
        ctx.fillText(formatPos(p), x, y + 15);
      }
    }
    ctx.font = o.font;
  }

  _colorOf(item, state) {
    if (state.colorBy === 'mapq') return mapqColor(item.mapq);
    if (state.colorBy === 'informative') return item.dk ? DK_COLOR : '#9a9a9a';
    if (state.colorBy === 'taxon') {
      const t = item.mol && item.mol.taxon;
      if (!t) return TAXON_COLORS.unclassified;
      return t.human ? TAXON_COLORS.human : TAXON_COLORS.nonhuman;
    }
    return EVENT_COLORS[item.event] || EVENT_COLORS['Deletion/Normal'];
  }

  _endX(panel, L, end) {
    // x of an end: in its own panel, or clamped to the edge of the nearest panel when outside every window
    if (end.wi >= 0) return { x1: this.xOf(L.panels[end.wi], end.pos), x2: this.xOf(L.panels[end.wi], end.end), outside: false };
    return null;
  }

  _drawSamplot(row, m, L, state, scales, data) {
    const ctx = this.ctx;
    const yTop = row.readsY + 10, yBottom = row.readsY + row.readsH - 6;
    const yOf = (insert, name, cross) => {
      if (insert === null || insert === undefined || cross && !insert) return yTop + 2 + hashJitter(name) * 6;
      const f = Math.min(1, insert / scales.maxInsert);
      const jitter = hashJitter(name) * 0.05 * f;
      return yBottom - (f + jitter) * (yBottom - yTop);
    };
    const items = [];
    for (const p of m.pairs) {
      if (state.hideProper && p.proper && !p.dk) continue;
      if (p.mapq < (state.minMapq || 0) && !p.dk) continue;
      items.push(p);
    }
    for (const s of m.splits) {
      if (s.mapq < (state.minMapq || 0) && !s.dk) continue;
      items.push(s);
    }
    // informative last, so they sit on top
    items.sort((a, b) => (a.dk ? 1 : 0) - (b.dk ? 1 : 0));
    for (const item of items) {
      const ends = item.kind === 'pair' ? item.ends : item.parts;
      const y = yOf(item.insert, item.name, item.cross);
      const color = this._colorOf(item, state);
      const dim = item.kind === 'pair' && item.proper && !item.dk;
      const alpha = dim ? 0.18 : (item.dk ? 1 : 0.5);
      const placed = ends.map((e) => (e.wi >= 0 ? { x1: this.xOf(L.panels[e.wi], e.pos), x2: this.xOf(L.panels[e.wi], e.end), e } : null));
      const present = placed.filter(Boolean);
      if (!present.length) continue;
      // connector
      const xs = present.map((q) => [q.x1, q.x2]).flat();
      const xa = Math.min(...xs), xb = Math.max(...xs);
      ctx.lineWidth = 1;
      if (item.dk && state.highlightDk !== false) {
        ctx.strokeStyle = DK_COLOR; ctx.globalAlpha = 0.35; ctx.lineWidth = 5;
        ctx.beginPath(); ctx.moveTo(xa, y); ctx.lineTo(xb, y); ctx.stroke();
        ctx.lineWidth = 1;
      }
      ctx.globalAlpha = alpha;
      ctx.strokeStyle = color;
      if (item.kind === 'split') ctx.setLineDash([3, 3]);
      ctx.beginPath(); ctx.moveTo(xa, y); ctx.lineTo(xb, y); ctx.stroke();
      ctx.setLineDash([]);
      // an end outside every window: a stub with an open marker at the panel edge it points to
      placed.forEach((q, i) => {
        if (q) return;
        const e = ends[i];
        const other = present[0];
        const panel = L.panels[other.e.wi];
        const left = e.chrom === panel.window.chrom && e.pos < panel.window.start;
        const right = e.chrom === panel.window.chrom && e.pos >= panel.window.end;
        const x = left ? panel.x : (right ? panel.x + panel.w : (other.x1 < panel.x + panel.w / 2 ? panel.x : panel.x + panel.w));
        ctx.beginPath(); ctx.moveTo(other.x1 < x ? other.x2 : other.x1, y); ctx.lineTo(x, y); ctx.stroke();
        ctx.beginPath(); ctx.arc(x, y, 2.5, 0, Math.PI * 2); ctx.stroke();
      });
      // the reads themselves: thick segments with a strand tick
      ctx.lineWidth = item.dk ? 3.5 : 2.5;
      ctx.fillStyle = color;
      for (const q of present) {
        ctx.beginPath(); ctx.moveTo(q.x1, y); ctx.lineTo(Math.max(q.x2, q.x1 + 2), y); ctx.stroke();
        if (item.kind === 'split') { ctx.beginPath(); ctx.arc(q.e.strand === '+' ? q.x2 : q.x1, y, 2.2, 0, Math.PI * 2); ctx.fill(); }
        else ctx.fillRect((q.e.strand === '+' ? q.x2 : q.x1) - 1.5, y - 2.5, 3, 5);
      }
      ctx.globalAlpha = 1;
      this.hits.push({ x1: xa - 3, x2: xb + 3, y1: y - 3, y2: y + 3, item, role: row.role });
    }
    // unmapped mates: a hollow marker at the mapped read's outer end, in the top band (no insert size)
    for (const u of m.unmappedMates || []) {
      if (u.mapq < (state.minMapq || 0) && !u.dk) continue;
      const x = this.xOf(L.panels[u.wi], u.pos);
      const y = yTop + 6 + hashJitter(u.name) * 8;
      ctx.globalAlpha = u.dk ? 1 : 0.6;
      ctx.strokeStyle = u.dk ? DK_COLOR : '#9467bd';
      ctx.lineWidth = 1.5;
      ctx.beginPath(); ctx.arc(x, y, 3, 0, Math.PI * 2); ctx.stroke();
      ctx.globalAlpha = 1;
      this.hits.push({ x1: x - 4, x2: x + 4, y1: y - 4, y2: y + 4, item: u, role: row.role });
    }
  }

  _drawPileup(row, m, L, state, data) {
    const ctx = this.ctx, o = this.options;
    const pitch = o.pileRow, barH = pitch - 2;
    for (let wi = 0; wi < L.panels.length; wi++) {
      const panel = L.panels[wi];
      const w = panel.window;
      // one row per molecule in this window, spanning its alignments here
      const spans = [];
      for (const mol of m.molecules.values()) {
        const alns = [...mol.r1, ...mol.r2].filter((a) => a.chrom === w.chrom && a.end > w.start && a.pos < w.end);
        if (!alns.length) continue;
        const dk = mol.dk;
        const mapq = mol.mapq ?? 0;
        if (mapq < (state.minMapq || 0) && !dk) continue;
        const proper = alns.every((a) => a.record && a.record.isProperPair && !a.supplementary && a.record.supplementaryAlignments.length === 0 && a.clips[0] < 20 && a.clips[1] < 20);
        if (state.hideProper && proper && !dk) continue;
        const s = Math.min(...alns.map((a) => a.pos - (a.clips ? a.clips[0] : 0)));
        const e = Math.max(...alns.map((a) => a.end + (a.clips ? a.clips[1] : 0)));
        spans.push({ mol, alns, s, e, dk, proper, mapq });
      }
      spans.sort((a, b) => (b.dk ? 1 : 0) - (a.dk ? 1 : 0) || a.s - b.s);
      const rowEnds = [];
      const placed = [];
      for (const sp of spans) {
        let r = rowEnds.findIndex((end) => end < sp.s);
        if (r < 0) { r = rowEnds.length; rowEnds.push(-Infinity); }
        if (r >= o.maxPileRows) continue;
        rowEnds[r] = sp.e + Math.max(10, (w.end - w.start) / panel.w * 2);
        placed.push([sp, r]);
      }
      const y0 = row.readsY + 4;
      for (const [sp, r] of placed) {
        const y = y0 + r * pitch;
        if (y + barH > row.readsY + row.readsH) continue;
        const alpha = sp.proper && !sp.dk ? 0.35 : 1;
        ctx.globalAlpha = alpha;
        // a line joining the molecule's alignments in this window
        const x1 = this.xOf(panel, Math.max(sp.s, w.start)), x2 = this.xOf(panel, Math.min(sp.e, w.end));
        ctx.strokeStyle = '#b5b7bd'; ctx.lineWidth = 1;
        ctx.beginPath(); ctx.moveTo(x1, y + barH / 2); ctx.lineTo(x2, y + barH / 2); ctx.stroke();
        for (const a of sp.alns) {
          const ax1 = this.xOf(panel, Math.max(a.pos, w.start)), ax2 = this.xOf(panel, Math.min(a.end, w.end));
          const item = { kind: 'alignment', mol: sp.mol, aln: a, event: sp.mol.event || 'Deletion/Normal', dk: sp.dk, mapq: a.mapq, name: sp.mol.name };
          const color = state.colorBy === 'event' ? (a.record && !a.record.isProperPair ? '#1b1b1f' : '#7d8089') : this._colorOf(item, state);
          if (sp.dk && state.highlightDk !== false) {
            ctx.fillStyle = DK_COLOR; ctx.fillRect(ax1 - 1, y - 1, Math.max(2, ax2 - ax1) + 2, barH + 2);
          }
          ctx.fillStyle = color;
          ctx.fillRect(ax1, y, Math.max(1.5, ax2 - ax1), barH);
          // soft clips as orange stubs (the clipped bases, to scale, capped)
          if (a.clips) {
            const bp = (w.end - w.start) / panel.w;
            const [cl, cr] = a.clips;
            ctx.fillStyle = '#D55E00';
            if (cl >= 1 && a.pos >= w.start) ctx.fillRect(ax1 - Math.min(40, cl / bp), y, Math.min(40, cl / bp), barH);
            if (cr >= 1 && a.end <= w.end) ctx.fillRect(ax2, y, Math.min(40, cr / bp), barH);
          }
          // strand
          ctx.fillStyle = '#ffffff';
          if (a.strand === '+') ctx.fillRect(ax2 - 2, y + 1, 1, barH - 2); else ctx.fillRect(ax1 + 1, y + 1, 1, barH - 2);
          this.hits.push({ x1: ax1, x2: Math.max(ax2, ax1 + 2), y1: y, y2: y + barH, item, role: row.role });
        }
        ctx.globalAlpha = 1;
      }
      if (rowEnds.length > o.maxPileRows) {
        ctx.fillStyle = o.muted; ctx.textAlign = 'right';
        ctx.fillText(`+${spans.length - placed.length} more`, panel.x + panel.w - 2, row.readsY + row.readsH - 4);
      }
    }
    // cross-window links of split reads and pairs, drawn over the rows as thin arcs between the parts
    ctx.lineWidth = 1;
    for (const item of [...m.splits, ...m.pairs]) {
      if (!item.cross) continue;
      const ends = item.kind === 'pair' ? item.ends : item.parts;
      const here = ends.filter((e) => e.wi >= 0);
      if (here.length < 2) continue;
      const y = row.readsY + row.readsH - 3;
      const xa = this.xOf(L.panels[here[0].wi], here[0].pos), xb = this.xOf(L.panels[here[here.length - 1].wi], here[here.length - 1].pos);
      ctx.strokeStyle = item.dk ? DK_COLOR : (EVENT_COLORS[item.event] || '#1b1b1f');
      ctx.globalAlpha = item.dk ? 0.8 : 0.35;
      if (item.kind === 'split') ctx.setLineDash([3, 3]);
      ctx.beginPath(); ctx.moveTo(xa, y); ctx.quadraticCurveTo((xa + xb) / 2, y - 24, xb, y); ctx.stroke();
      ctx.setLineDash([]);
      ctx.globalAlpha = 1;
    }
  }

  _drawCaptions(row, m, L) {
    const ctx = this.ctx, o = this.options;
    ctx.fillStyle = o.muted; ctx.textAlign = 'right';
    ctx.font = `10px ${o.font.split('px ')[1]}`;
    for (let wi = 0; wi < L.panels.length; wi++) {
      const s = m.stats[wi];
      const panel = L.panels[wi];
      const text = `${s.discordant} disc · ${s.split} split · ${s.clipped} clip` + (s.informative ? ` · ${s.informative} inf` : '');
      ctx.fillText(text, panel.x + panel.w - 2, row.labelY);
    }
    ctx.font = o.font;
  }

  /** The item under (x, y) in CSS pixels, or null. */
  hitTest(x, y) {
    let best = null, bestD = Infinity;
    for (const h of this.hits) {
      if (x < h.x1 || x > h.x2 || y < h.y1 || y > h.y2) continue;
      const d = Math.abs((h.y1 + h.y2) / 2 - y);
      if (d < bestD) { bestD = d; best = h; }
    }
    return best;
  }
}
