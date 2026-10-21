// app.js - the viewer: the candidate table, the plot of the selected candidate over its bundle, the verdicts.
// Data come from a manifest.json (svreview site) next to the mini-BAMs: over HTTP (a static host, or
// `svreview serve`, which also takes the verdicts), or from a folder opened in the browser.

import { BamFile, RemoteSource, FileSource } from './bam.js';
import { buildMember, captionOf, consensusBreakpoint } from './reads.js';
import { Plot, EVENT_COLORS, DK_COLOR, TAXON_COLORS, formatPos, formatBp } from './render.js';

const VERDICTS = ['de_novo', 'inherited', 'artifact', 'unsure'];
const VERDICT_LABELS = { de_novo: 'de novo', inherited: 'inherited', artifact: 'artifact', unsure: 'unsure' };
const VERDICT_COLUMNS = ['candidate', 'trio', 'chrom', 'start', 'end', 'sv_type', 'verdict', 'note', 'reviewer', 'reviewed_at'];
const ROLES = ['child', 'father', 'mother'];

const $ = (sel) => document.querySelector(sel);

class DataSource {
  /** `base`: a URL (ending in /) or a Map of relative path -> File. */
  constructor(base) {
    this.base = base;
    this.isFiles = base instanceof Map;
    this.bams = new Map();
  }
  async json(path) {
    if (this.isFiles) {
      const f = this.base.get(path);
      if (!f) throw new Error(`${path} is not in the folder`);
      return JSON.parse(await f.text());
    }
    const res = await fetch(this.base + path, { cache: 'no-cache' });
    if (!res.ok) throw new Error(`${this.base + path}: HTTP ${res.status}`);
    return res.json();
  }
  async text(path) {
    if (this.isFiles) {
      const f = this.base.get(path);
      return f ? f.text() : null;
    }
    const res = await fetch(this.base + path, { cache: 'no-cache' });
    return res.ok ? res.text() : null;
  }
  bam(path) {
    let b = this.bams.get(path);
    if (!b) {
      if (this.isFiles) {
        const f = this.base.get(path), i = this.base.get(path + '.bai') || this.base.get(path.replace(/\.bam$/, '.bai'));
        if (!f || !i) throw new Error(`${path} or its .bai is not in the folder`);
        b = new BamFile(new FileSource(f), new FileSource(i));
      } else {
        b = new BamFile(new RemoteSource(this.base + path), new RemoteSource(this.base + path + '.bai'));
      }
      this.bams.set(path, b);
    }
    return b;
  }
  async postVerdicts(rows) {
    if (this.isFiles) return false;
    const body = rows.map((r) => VERDICT_COLUMNS.map((c) => String(r[c] ?? '').replace(/[\t\n]/g, ' ')).join('\t')).join('\n') + '\n';
    try {
      const res = await fetch(this.base + 'verdicts.tsv', { method: 'POST', headers: { 'Content-Type': 'text/tab-separated-values' }, body });
      return res.ok;
    } catch (e) {
      return false;
    }
  }
}

/** A member's junction evidence as molecules, each counted once (a pair in two windows is one pair). */
function memberTotals(m) {
  const clipped = new Set(m.clips.map((c) => c.mol.name));
  let informative = 0;
  for (const mol of m.molecules.values()) if (mol.dk) informative++;
  return { discordant: m.pairs.filter((p) => !p.proper).length, split: m.splits.length, clipped: clipped.size, informative };
}

class App {
  constructor() {
    this.manifest = null;
    this.source = null;
    this.dataKey = null;
    this.candidates = [];
    this.selected = null;
    this.verdicts = new Map();
    this.plot = new Plot($('#plot'));
    this.state = { mode: 'samplot', colorBy: 'event', hideProper: false, minMapq: 0, highlightDk: true, pileRows: 24 };
    this.loadToken = 0;
    this.plotData = null;
    this.filter = { tiers: new Set([1, 2, 3, 4]), text: '', verdict: 'all' };
    this.sort = { key: 'rank', dir: 1 };
    this._wire();
  }

  // ---- loading ---------------------------------------------------------------------------------------

  async openUrl(dataUrl) {
    let base = dataUrl;
    if (!base.endsWith('/')) base += '/';
    this.source = new DataSource(base);
    this.dataKey = base;
    await this._load();
  }

  async openFiles(fileList) {
    const files = new Map();
    let root = null;
    for (const f of fileList) {
      const rel = f.webkitRelativePath || f.name;
      const parts = rel.split('/');
      if (root === null && parts.length > 1) root = parts[0];
      const key = parts.length > 1 ? parts.slice(1).join('/') : rel;
      files.set(key, f);
    }
    if (!files.has('manifest.json')) throw new Error('the folder has no manifest.json (svreview site writes one)');
    this.source = new DataSource(files);
    this.dataKey = 'folder:' + (root || 'local');
    await this._load();
  }

  async _load() {
    this.setStatus('loading the manifest…');
    this.manifest = await this.source.json('manifest.json');
    if (this.manifest.contract !== 1) this.setStatus(`manifest contract ${this.manifest.contract}: this viewer reads 1`, true);
    this.candidates = this.manifest.candidates.map((c) => Object.assign({}, c));
    document.title = `svreview · ${this.manifest.title || ''}`;
    $('#title').textContent = this.manifest.title || 'svreview';
    $('#subtitle').textContent = `${this.candidates.length} candidates · ${Object.keys(this.manifest.trios).length} trio(s) · built ${this.manifest.created || ''}`;
    this._loadVerdicts();
    const served = await this.source.text('verdicts.tsv');
    if (served) this._mergeVerdictText(served, 'server');
    this.renderTable();
    const want = new URLSearchParams(location.search).get('c');
    const first = this.visibleCandidates()[0];
    const pick = (want && this.candidates.find((c) => c.id === want)) || first;
    if (pick) await this.select(pick.id);
    else this.setStatus('no candidate to show');
    $('#app').classList.remove('empty');
  }

  // ---- the table -------------------------------------------------------------------------------------

  visibleCandidates() {
    const text = this.filter.text.trim().toLowerCase();
    const out = this.candidates.filter((c) => {
      if (!this.filter.tiers.has(c.review_tier ?? 4)) return false;
      const v = this.verdicts.get(c.id);
      if (this.filter.verdict === 'unreviewed' && v) return false;
      if (this.filter.verdict === 'reviewed' && !v) return false;
      if (text) {
        const hay = `${c.id} ${c.trio} ${c.sv_type} ${c.class} ${(c.flags || []).join(' ')} ${c.chrom}:${c.start}-${c.end} ${v ? v.verdict + ' ' + v.note : ''}`.toLowerCase();
        if (!hay.includes(text)) return false;
      }
      return true;
    });
    const key = this.sort.key, dir = this.sort.dir;
    const val = (c) => {
      switch (key) {
        case 'rank': return [c.review_tier ?? 9, c.rank ?? 1e9];
        case 'locus': return [c.chrom, c.start];
        case 'type': return [c.sv_type, c.rank ?? 1e9];
        case 'reads': return [-(c.evidence?.breakpoint_reads ?? 0), -(c.evidence?.reads ?? 0)];
        case 'trio': return [c.trio, c.rank ?? 1e9];
        case 'verdict': return [this.verdicts.get(c.id)?.verdict || 'zz', c.rank ?? 1e9];
        default: return [c.rank ?? 1e9];
      }
    };
    out.sort((a, b) => {
      const va = val(a), vb = val(b);
      for (let i = 0; i < va.length; i++) {
        if (va[i] < vb[i]) return -dir;
        if (va[i] > vb[i]) return dir;
      }
      return 0;
    });
    return out;
  }

  renderTable() {
    const tbody = $('#candidates tbody');
    tbody.innerHTML = '';
    const rows = this.visibleCandidates();
    for (const c of rows) {
      const tr = document.createElement('tr');
      tr.dataset.id = c.id;
      if (this.selected && this.selected.id === c.id) tr.classList.add('selected');
      const v = this.verdicts.get(c.id);
      const flags = (c.flags || []).join(' ');
      const ev = c.evidence || {};
      const nj = (c.junctions || []).length;
      tr.innerHTML = `
        <td class="tier t${c.review_tier ?? 4}">${c.review_tier ?? '·'}</td>
        <td class="num">${c.rank ?? '·'}</td>
        <td class="mono">${c.trio}</td>
        <td class="mono locus">${c.chrom}:${formatPos(c.start)}-${formatPos(c.end)}</td>
        <td><span class="type ${c.sv_type}">${c.sv_type}</span>${nj ? `<span class="junc" title="${nj} junction(s)">${nj}</span>` : ''}</td>
        <td class="num" title="breakpoint reads / reads">${ev.breakpoint_reads ?? '·'}/${ev.reads ?? '·'}</td>
        <td class="flags" title="${flags}">${flags || (c.bundle ? '' : '<span class="warn">no bundle</span>')}</td>
        <td class="verdict ${v ? v.verdict : ''}">${v ? VERDICT_LABELS[v.verdict] : ''}</td>`;
      tr.addEventListener('click', () => this.select(c.id));
      tbody.appendChild(tr);
    }
    $('#table-count').textContent = `${rows.length} of ${this.candidates.length}`;
    const reviewed = [...this.verdicts.keys()].filter((id) => this.candidates.some((c) => c.id === id)).length;
    $('#reviewed-count').textContent = `${reviewed} reviewed`;
  }

  // ---- selection and the plot ------------------------------------------------------------------------

  async select(id) {
    const c = this.candidates.find((x) => x.id === id);
    if (!c) return;
    this.selected = c;
    const token = ++this.loadToken;
    const url = new URL(location.href);
    url.searchParams.set('c', id);
    history.replaceState(null, '', url);
    for (const tr of document.querySelectorAll('#candidates tbody tr')) tr.classList.toggle('selected', tr.dataset.id === id);
    const sel = document.querySelector('#candidates tbody tr.selected');
    if (sel) sel.scrollIntoView({ block: 'nearest' });
    this.renderHeader(c);
    this.renderVerdictBox(c);
    if (!c.bundle || !c.members || !Object.keys(c.members).length) {
      this.plotData = null;
      this.setStatus(c.missing || 'this candidate has no bundle: nothing to draw');
      this.plot.ctx.clearRect(0, 0, this.plot.canvas.width, this.plot.canvas.height);
      return;
    }
    this.setStatus('reading the bundle…');
    try {
      const members = {};
      const breakpoints = [];
      const taxa = await this._taxaOf(c);
      for (const role of ROLES) {
        const m = c.members[role];
        if (!m) { members[role] = null; continue; }
        const bam = this.source.bam(m.bam);
        const perWindow = [];
        for (const w of c.windows) perWindow.push(await bam.records(w.chrom, w.start, w.end));
        if (token !== this.loadToken) return;
        const built = buildMember(c.windows, perWindow, { taxa: taxa && taxa[role] });
        members[role] = built;
        for (let wi = 0; wi < c.windows.length; wi++) {
          const b = consensusBreakpoint(built.clips, wi);
          if (b) breakpoints.push({ role, wi, pos: b.pos, n: b.n });
        }
      }
      this.plotData = { windows: c.windows, members, regions: c.regions, breakpoints, samples: this.manifest.trios[c.trio]?.members || {} };
      this.redraw();
      this.renderMemberSummary(c, members);
      this.setStatus('');
    } catch (e) {
      console.error(e);
      this.setStatus(`could not read the bundle: ${e.message}`, true);
    }
  }

  async _taxaOf(c) {
    // optional per-read taxa sidecars: <bundle>/<role>.taxa.tsv (read, taxid, name, human) - absent in most runs
    if (!c.taxa) return null;
    const out = {};
    for (const role of ROLES) {
      if (!c.taxa[role]) continue;
      const text = await this.source.text(c.taxa[role]);
      if (!text) continue;
      const map = new Map();
      for (const line of text.split('\n')) {
        if (!line || line.startsWith('read\t')) continue;
        const f = line.split('\t');
        map.set(f[0], { taxid: f[1], name: f[2], human: f[3] === '1' || f[3] === 'true' });
      }
      out[role] = map;
    }
    return out;
  }

  redraw() {
    if (!this.plotData) return;
    const dark = matchMedia('(prefers-color-scheme: dark)').matches;
    this.plot.setTheme(dark
      ? { ink: '#e8e8ec', muted: '#a3a5ad', line: '#2f323a', background: '#15171c', covHigh: 'rgba(190,190,200,0.55)', covLow: 'rgba(120,120,130,0.4)' }
      : { ink: '#1b1b1f', muted: '#6b6e76', line: '#d6d8de', background: '#ffffff', covHigh: 'rgba(90,90,96,0.55)', covLow: 'rgba(150,150,156,0.35)' });
    this.plot.draw(this.plotData, this.state);
    this.renderLegend();
  }

  renderHeader(c) {
    const ev = c.evidence || {};
    const trio = this.manifest.trios[c.trio] || {};
    const nhfStatus = trio.nhf?.status || 'none';
    const own = (c.regions || []).find((r) => r.own) || {};
    const screen = own.nhf_screen;
    let nhfText;
    if (nhfStatus !== 'screened') nhfText = `screen ${nhfStatus}${trio.nhf?.reason ? ': ' + trio.nhf.reason : ''}`;
    else if (!screen) nhfText = 'region not screened';
    else nhfText = `${(screen.fractions?.nonhuman ?? 0).toFixed(3)} non-human of ${screen.reads} reads` + (screen.top_taxon && screen.top_taxon !== '.' ? ` · top ${screen.top_taxon}` : '');
    const nhfClass = nhfStatus !== 'screened' ? 'unknown' : (screen && (screen.fractions?.nonhuman ?? 0) >= 0.05 ? 'bad' : 'ok');
    const junctions = (c.junctions || []).map((j) => `<span class="chip j">${j.chrom}:${formatPos(j.start)}-${formatPos(j.end)} <b>${j.strands}</b> ${j.sv_type} · ${j.reads} reads</span>`).join(' ');
    const flags = [...(c.locus_flags || []).map((f) => `<span class="chip flag">${f}</span>`),
      ...(c.partner_flags || []).map((f) => `<span class="chip flag">partner:${f}</span>`),
      ...(c.flags || []).filter((f) => !(c.locus_flags || []).includes(f) && !f.startsWith('partner:')).map((f) => `<span class="chip flag">${f}</span>`)];
    const uniqueFlags = [...new Set(flags)].join(' ');
    $('#cand-title').innerHTML = `<span class="tier t${c.review_tier ?? 4}">tier ${c.review_tier ?? '·'}</span> rank ${c.rank ?? '·'} · <span class="mono">${c.trio}</span> · <span class="mono">${c.chrom}:${formatPos(c.start + 1)}-${formatPos(c.end)}</span> · <span class="type ${c.sv_type}">${c.sv_type}</span> ${c.class}${c.shared_with ? ` <span class="muted">(bundle shared with ${c.shared_with})</span>` : ''}`;
    $('#cand-evidence').innerHTML = [
      ['reads', ev.reads], ['k-mers', ev.unique_kmers], ['split', ev.split_reads], ['discordant', ev.discordant_pairs],
      ['unmapped mates', ev.unmapped_mates], ['breakpoint', ev.breakpoint_reads], ['indel ≥50', ev.large_indel_reads],
      ['max clip', ev.max_clip_len], ['max MAPQ', ev.max_mapq], ['MAPQ 0', ev.mapq0_reads],
    ].filter(([, v]) => v !== undefined && v !== null).map(([k, v]) => `<span class="stat"><b>${v}</b> ${k}</span>`).join('');
    $('#cand-junctions').innerHTML = junctions || '<span class="muted">no junction: clipped reads alone</span>';
    $('#cand-flags').innerHTML = uniqueFlags || '<span class="muted">no flags</span>';
    $('#cand-nhf').innerHTML = `<span class="chip nhf ${nhfClass}">non-human fraction: ${nhfText}</span>`;
    const inf = c.members?.child?.informative;
    $('#cand-informative').innerHTML = inf !== null && inf !== undefined
      ? `<span class="chip inf">${inf} informative reads of the child in the windows</span>`
      : `<span class="chip muted">informative BAM not joined (no dk tags)</span>`;
    $('#cand-windows').innerHTML = (c.windows || []).map((w) => `<span class="chip w">${w.chrom}:${formatPos(w.start + 1)}-${formatPos(w.end)} (${formatBp(w.end - w.start)})</span>`).join(' ');
  }

  renderMemberSummary(c, members) {
    const parts = [];
    for (const role of ROLES) {
      const m = members[role];
      if (!m) continue;
      const tot = memberTotals(m);
      parts.push(`<span class="member"><b>${role}</b> ${tot.discordant} discordant · ${tot.split} split · ${tot.clipped} clipped${tot.informative ? ` · <span class="dk">${tot.informative} informative</span>` : ''}</span>`);
    }
    const child = members.child, father = members.father, mother = members.mother;
    let reading = '';
    if (child && father && mother) {
      const sum = (m) => { const t = memberTotals(m); return t.discordant + t.split + t.clipped; };
      const cs = sum(child), ps = sum(father) + sum(mother);
      if (cs >= 2 && ps === 0) reading = 'the junction-crossing reads are the child\'s only: consistent with de novo (or a somatic change of the child\'s cells)';
      else if (cs >= 2 && ps >= 2) reading = 'a parent shows the same kind of reads at these windows: look for an inherited variant';
      else if (cs < 2) reading = 'few junction-crossing reads in the child at the drawn MAPQ: weak';
    }
    $('#member-summary').innerHTML = parts.join(' ') + (reading ? `<div class="reading">${reading}</div>` : '');
  }

  renderLegend() {
    const s = this.state;
    const sw = (color, dashed) => `<span class="sw" style="background:${color};${dashed ? 'opacity:.6' : ''}"></span>`;
    let items = [];
    if (s.colorBy === 'event') {
      items = [[sw(EVENT_COLORS['Deletion/Normal']), 'deletion / normal (+ −)'], [sw(EVENT_COLORS.Duplication), 'duplication (− +)'],
        [sw(EVENT_COLORS.Inversion), 'inversion (+ + or − −)'], [sw(EVENT_COLORS.InterChrm), 'between chromosomes']];
    } else if (s.colorBy === 'mapq') {
      items = [[sw('rgb(205,205,205)'), 'MAPQ 0'], [sw('rgb(145,145,152)'), 'MAPQ 20'], [sw('rgb(25,25,45)'), 'MAPQ 60']];
    } else if (s.colorBy === 'informative') {
      items = [[sw(DK_COLOR), 'informative (proband-unique k-mers)'], [sw('#9a9a9a'), 'other reads']];
    } else {
      items = [[sw(TAXON_COLORS.human), 'human lineage'], [sw(TAXON_COLORS.nonhuman), 'non-human'], [sw(TAXON_COLORS.unclassified), 'unclassified / no sidecar']];
    }
    items.push([`<span class="sw halo"></span>`, 'halo: the child\'s informative reads (dk)']);
    items.push([`<span class="sw line"></span>`, s.mode === 'pileup' ? 'bar: an alignment; orange stub: soft clip' : 'solid: read pair; dotted: split read']);
    items.push([`<span class="sw" style="background:rgba(0,114,178,0.3)"></span>`, 'candidate region']);
    items.push([`<span class="sw" style="background:rgba(213,94,0,0.3)"></span>`, 'partner region']);
    items.push([`<span class="sw" style="background:rgba(230,159,0,0.6)"></span>`, 'coverage of informative reads']);
    if (s.mode !== 'pileup') items.push([`<span class="sw ring"></span>`, 'ring: a read whose mate is unmapped']);
    $('#legend').innerHTML = items.map(([k, v]) => `<span class="li">${k}${v}</span>`).join('');
  }

  // ---- verdicts --------------------------------------------------------------------------------------

  get storageKey() { return `svreview:verdicts:${this.dataKey}`; }

  _loadVerdicts() {
    this.verdicts = new Map();
    try {
      const raw = localStorage.getItem(this.storageKey);
      if (raw) for (const [k, v] of Object.entries(JSON.parse(raw))) this.verdicts.set(k, v);
    } catch (e) { /* storage unavailable: verdicts live for the page */ }
    try { $('#reviewer').value = localStorage.getItem('svreview:reviewer') || ''; } catch (e) { /* ignore */ }
  }

  _saveLocal() {
    try { localStorage.setItem(this.storageKey, JSON.stringify(Object.fromEntries(this.verdicts))); } catch (e) { /* ignore */ }
  }

  _mergeVerdictText(text, origin) {
    const lines = text.split('\n').filter((l) => l.trim());
    if (!lines.length) return 0;
    let header = lines[0].split('\t');
    let start = 1;
    if (!header.includes('verdict')) { header = VERDICT_COLUMNS; start = 0; }
    let n = 0;
    for (const line of lines.slice(start)) {
      const f = line.split('\t');
      const row = {};
      header.forEach((h, i) => { row[h] = f[i] ?? ''; });
      const id = `${row.trio}/${row.chrom}_${row.start}_${row.end}`;
      if (!row.verdict) { this.verdicts.delete(id); continue; }
      if (!VERDICTS.includes(row.verdict)) continue;
      const cur = this.verdicts.get(id);
      if (cur && cur.reviewed_at && row.reviewed_at && cur.reviewed_at > row.reviewed_at) continue;
      this.verdicts.set(id, { verdict: row.verdict, note: row.note || '', reviewer: row.reviewer || '', reviewed_at: row.reviewed_at || '', origin });
      n++;
    }
    this._saveLocal();
    return n;
  }

  verdictRow(c, v) {
    return { candidate: `${c.chrom}:${c.start}-${c.end}`, trio: c.trio, chrom: c.chrom, start: c.start, end: c.end, sv_type: c.sv_type,
      verdict: v.verdict, note: v.note, reviewer: v.reviewer, reviewed_at: v.reviewed_at };
  }

  async setVerdict(verdict) {
    const c = this.selected;
    if (!c) return;
    const note = $('#note').value.trim();
    const reviewer = $('#reviewer').value.trim();
    try { localStorage.setItem('svreview:reviewer', reviewer); } catch (e) { /* ignore */ }
    const now = new Date().toISOString().replace(/\.\d+Z$/, '');
    let row;
    if (verdict === null) {
      this.verdicts.delete(c.id);
      row = Object.assign(this.verdictRow(c, { verdict: '', note: '', reviewer, reviewed_at: now }));
    } else {
      const v = { verdict, note, reviewer, reviewed_at: now };
      this.verdicts.set(c.id, v);
      row = this.verdictRow(c, v);
    }
    this._saveLocal();
    this.renderVerdictBox(c);
    this.renderTable();
    const ok = await this.source.postVerdicts([row]);
    $('#verdict-sync').textContent = ok ? 'written to verdicts.tsv on the server' : 'kept in this browser (export the TSV to hand it over)';
    $('#verdict-sync').className = ok ? 'ok' : 'local';
  }

  renderVerdictBox(c) {
    const v = this.verdicts.get(c.id);
    for (const b of document.querySelectorAll('#verdict-buttons button')) b.classList.toggle('active', !!v && b.dataset.verdict === v.verdict);
    $('#note').value = v ? v.note : '';
    $('#verdict-meta').textContent = v ? `${VERDICT_LABELS[v.verdict]} · ${v.reviewer || 'anonymous'} · ${v.reviewed_at}` : 'no verdict yet';
  }

  exportVerdicts() {
    const rows = [VERDICT_COLUMNS.join('\t')];
    for (const c of this.candidates) {
      const v = this.verdicts.get(c.id);
      if (!v) continue;
      const r = this.verdictRow(c, v);
      rows.push(VERDICT_COLUMNS.map((k) => String(r[k] ?? '').replace(/[\t\n]/g, ' ')).join('\t'));
    }
    const blob = new Blob([rows.join('\n') + '\n'], { type: 'text/tab-separated-values' });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = `verdicts.${(this.manifest.title || 'svreview').replace(/[^\w.-]+/g, '_')}.tsv`;
    a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  }

  async importVerdicts(file) {
    const n = this._mergeVerdictText(await file.text(), 'import');
    this.renderTable();
    if (this.selected) this.renderVerdictBox(this.selected);
    this.setStatus(`${n} verdicts read from ${file.name}`);
  }

  // ---- navigation and UI -----------------------------------------------------------------------------

  step(delta) {
    const rows = this.visibleCandidates();
    if (!rows.length) return;
    const i = this.selected ? rows.findIndex((c) => c.id === this.selected.id) : -1;
    const j = Math.max(0, Math.min(rows.length - 1, i + delta));
    if (rows[j]) this.select(rows[j].id);
  }

  setStatus(text, error = false) {
    const el = $('#status');
    el.textContent = text;
    el.classList.toggle('error', !!error);
  }

  _wire() {
    const canvas = this.plot.canvas;
    const tip = $('#tooltip');
    canvas.addEventListener('mousemove', (ev) => {
      const rect = canvas.getBoundingClientRect();
      const hit = this.plot.hitTest(ev.clientX - rect.left, ev.clientY - rect.top);
      if (!hit) { tip.hidden = true; return; }
      tip.innerHTML = this._describe(hit);
      tip.hidden = false;
      const x = ev.clientX + 14, y = ev.clientY + 14;
      tip.style.left = Math.min(x, window.innerWidth - tip.offsetWidth - 10) + 'px';
      tip.style.top = Math.min(y, window.innerHeight - tip.offsetHeight - 10) + 'px';
    });
    canvas.addEventListener('mouseleave', () => { tip.hidden = true; });
    window.addEventListener('resize', () => this.redraw());
    matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => this.redraw());
    $('#mode').addEventListener('change', (e) => { this.state.mode = e.target.value; this.redraw(); });
    $('#colorBy').addEventListener('change', (e) => { this.state.colorBy = e.target.value; this.redraw(); });
    $('#minMapq').addEventListener('change', (e) => { this.state.minMapq = parseInt(e.target.value, 10) || 0; this.redraw(); });
    $('#hideProper').addEventListener('change', (e) => { this.state.hideProper = e.target.checked; this.redraw(); });
    $('#highlightDk').addEventListener('change', (e) => { this.state.highlightDk = e.target.checked; this.redraw(); });
    $('#pileRows').addEventListener('change', (e) => { this.state.pileRows = parseInt(e.target.value, 10) || 24; this.redraw(); });
    for (const b of document.querySelectorAll('#verdict-buttons button[data-verdict]')) {
      b.addEventListener('click', () => this.setVerdict(b.dataset.verdict));
    }
    $('#clear-verdict').addEventListener('click', () => this.setVerdict(null));
    $('#note').addEventListener('change', () => { const v = this.selected && this.verdicts.get(this.selected.id); if (v) this.setVerdict(v.verdict); });
    $('#export').addEventListener('click', () => this.exportVerdicts());
    $('#import').addEventListener('change', (e) => { if (e.target.files[0]) this.importVerdicts(e.target.files[0]); e.target.value = ''; });
    $('#open-folder').addEventListener('change', (e) => { if (e.target.files.length) this.openFiles(e.target.files).catch((err) => this.setStatus(err.message, true)); });
    $('#filter-text').addEventListener('input', (e) => { this.filter.text = e.target.value; this.renderTable(); });
    $('#filter-verdict').addEventListener('change', (e) => { this.filter.verdict = e.target.value; this.renderTable(); });
    for (const cb of document.querySelectorAll('#tier-filter input')) {
      cb.addEventListener('change', () => {
        this.filter.tiers = new Set([...document.querySelectorAll('#tier-filter input:checked')].map((x) => parseInt(x.value, 10)));
        this.renderTable();
      });
    }
    for (const th of document.querySelectorAll('#candidates th[data-sort]')) {
      th.addEventListener('click', () => {
        if (this.sort.key === th.dataset.sort) this.sort.dir = -this.sort.dir; else { this.sort.key = th.dataset.sort; this.sort.dir = 1; }
        this.renderTable();
      });
    }
    $('#prev').addEventListener('click', () => this.step(-1));
    $('#next').addEventListener('click', () => this.step(1));
    document.addEventListener('keydown', (ev) => {
      if (ev.target.matches('input, textarea, select')) { if (ev.key === 'Escape') ev.target.blur(); return; }
      if (ev.key === 'j' || ev.key === 'ArrowDown') { ev.preventDefault(); this.step(1); }
      else if (ev.key === 'k' || ev.key === 'ArrowUp') { ev.preventDefault(); this.step(-1); }
      else if (ev.key >= '1' && ev.key <= '4') this.setVerdict(VERDICTS[parseInt(ev.key, 10) - 1]);
      else if (ev.key === '0') this.setVerdict(null);
      else if (ev.key === 'n') { ev.preventDefault(); $('#note').focus(); }
      else if (ev.key === 'p') { $('#mode').value = this.state.mode === 'pileup' ? 'samplot' : 'pileup'; this.state.mode = $('#mode').value; this.redraw(); }
      else if (ev.key === 'e') this.exportVerdicts();
      else if (ev.key === '?') $('#help').toggleAttribute('hidden');
    });
    $('#help-button').addEventListener('click', () => $('#help').toggleAttribute('hidden'));
    $('#help-close').addEventListener('click', () => $('#help').setAttribute('hidden', ''));
  }

  _describe(hit) {
    const it = hit.item;
    const mol = it.mol;
    const line = (e) => `${e.chrom}:${formatPos(e.pos + 1)}-${formatPos(e.end)} ${e.strand} MAPQ ${e.mapq ?? '?'} <span class="mono">${e.cigar}</span>${e.stub ? ' (mate, from the record)' : ''}${e.fromSA ? ' (from SA)' : ''}`;
    let body = '';
    if (it.kind === 'pair') body = `<b>read pair</b> · ${it.event}${it.insert ? ` · insert ${formatPos(it.insert)} bp` : ''}${it.proper ? ' · proper' : ' · discordant'}<br>` + it.ends.map(line).join('<br>');
    else if (it.kind === 'split') body = `<b>split read</b> · ${it.event}${it.insert ? ` · span ${formatPos(it.insert)} bp` : ''}<br>` + it.parts.map((p) => line(p) + (p.clipSide ? ` · clipped ${p.clipSide}` : '')).join('<br>');
    else if (it.kind === 'unmapped') body = `<b>unmapped mate</b> of ${line(it.aln)}: sequence the reference lacks, placed here`;
    else body = `<b>alignment</b> ${line(it.aln)}${it.aln.clips ? ` · clips ${it.aln.clips[0]}/${it.aln.clips[1]}` : ''}${it.aln.record && it.aln.record.tlen ? ` · TLEN ${it.aln.record.tlen}` : ''}`;
    const tags = [];
    if (it.dk) tags.push('<span class="dk">informative (dk)</span>');
    if (mol && mol.taxon) tags.push(`taxon ${mol.taxon.name || mol.taxon.taxid}`);
    return `<div class="tt-name mono">${hit.role} · ${it.name}</div>${body}${tags.length ? '<br>' + tags.join(' · ') : ''}`;
  }
}

async function main() {
  const app = new App();
  window.svreview = app;
  const params = new URLSearchParams(location.search);
  const data = params.get('data') || (location.protocol === 'file:' ? null : 'demo');
  if (data) {
    try { await app.openUrl(data); } catch (e) { console.error(e); app.setStatus(`could not open ${data}: ${e.message}`, true); }
  } else {
    app.setStatus('open a site folder (svreview site) to begin');
  }
}

main();
