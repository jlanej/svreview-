// reads.js - from BAM records to what the plot draws: molecules (a read name's alignments), paired ends,
// split alignments, clips at a breakpoint and coverage, over the candidate's windows (its own region's and
// each junction partner's). The classification follows samplot (Belyeu et al. 2018) - a pair's event type
// from the strands of its two ends ordered by position, its insert size as the outer distance - extended
// to several windows, so that a pair or a split whose ends fall in two windows is drawn across them.

import { clippedSide, cigarReferenceLength } from './bam.js';

export const EVENTS = {
  NORMAL: 'Deletion/Normal', DEL: 'Deletion/Normal', DUP: 'Duplication', INV: 'Inversion', BND: 'InterChrm',
};
export const CLIP_MIN = 20;   // the tool's breakpoint_reads rule: a soft clip of at least 20 bp

/** Which window (index) a position falls in, or -1. */
export function windowOf(windows, chrom, pos) {
  for (let i = 0; i < windows.length; i++) {
    const w = windows[i];
    if (w.chrom === chrom && pos >= w.start && pos < w.end) return i;
  }
  return -1;
}

/** samplot's event type from two strands ordered by position (true = forward). */
export function eventType(firstForward, secondForward) {
  if (firstForward && !secondForward) return EVENTS.NORMAL;
  if (!firstForward && secondForward) return EVENTS.DUP;
  return EVENTS.INV;
}

function isInformative(r) {
  return r.tags.dk === 1 || r.tags.dk === '1' || r.tags.DV !== undefined;
}

/**
 * Build a member's plot data from the records fetched in each window.
 * @param windows [{chrom, start, end}] (0-based, half-open)
 * @param recordsPerWindow [[BamRecord]] one array per window, as BamFile.records() returns them
 * @param opts { minMapq (reads below are not drawn, default 0), taxa (Map name -> label) }
 */
export function buildMember(windows, recordsPerWindow, opts = {}) {
  const minMapq = opts.minMapq ?? 0;
  const taxa = opts.taxa || null;
  const coverage = windows.map((w) => ({
    high: new Float32Array(w.end - w.start), low: new Float32Array(w.end - w.start), dk: new Float32Array(w.end - w.start),
  }));
  const molecules = new Map();
  const seen = new Set();
  const stats = windows.map(() => ({ reads: 0, discordant: 0, split: 0, clipped: 0, informative: 0, unmapped: 0 }));
  for (let wi = 0; wi < windows.length; wi++) {
    const w = windows[wi];
    for (const r of recordsPerWindow[wi]) {
      const key = `${r.name}/${r.flag & 0xc0}/${r.refId}:${r.pos}/${r.flag & 0x900}/${r.cigarString}`;
      if (seen.has(key)) continue;        // a record fetched in two overlapping windows: once
      seen.add(key);
      if (r.isSecondary || r.isDuplicate || r.isQCFail) continue;
      const dk = isInformative(r);
      if (r.isUnmapped) {
        stats[wi].unmapped++;
        if (dk) stats[wi].informative++;
        const mol = getMolecule(molecules, r.name);
        mol.unmapped.push({ record: r, wi, pos: r.pos, dk });
        if (dk) mol.dk = true;
        continue;
      }
      // coverage over the window (every window the alignment touches)
      for (let i = 0; i < windows.length; i++) {
        const ww = windows[i];
        if (ww.chrom !== r.chrom) continue;
        const cov = coverage[i];
        for (const [bs, be] of r.blocks) {
          const s = Math.max(bs, ww.start), e = Math.min(be, ww.end);
          for (let p = s; p < e; p++) {
            const idx = p - ww.start;
            if (r.mapq >= 20) cov.high[idx]++; else cov.low[idx]++;
            if (dk) cov.dk[idx]++;
          }
        }
      }
      if (r.mapq < minMapq) continue;
      if (windowOf(windows, r.chrom, r.pos) === wi || windowOf(windows, r.chrom, Math.max(r.pos, r.end - 1)) === wi) {
        stats[wi].reads++;
        if (dk) stats[wi].informative++;
      }
      const mol = getMolecule(molecules, r.name);
      if (dk) mol.dk = true;
      if (taxa && taxa.has(r.name)) mol.taxon = taxa.get(r.name);
      const aln = {
        record: r, wi: windowOf(windows, r.chrom, r.pos), chrom: r.chrom, pos: r.pos, end: r.end, strand: r.strand,
        mapq: r.mapq, cigar: r.cigarString, clips: r.softClips, supplementary: r.isSupplementary, dk,
        readNumber: r.readNumber,
      };
      const list = r.readNumber === 2 ? mol.r2 : mol.r1;
      list.push(aln);
      if (mol.mapq === null || r.mapq > mol.mapq) mol.mapq = r.mapq;
    }
  }
  const pairs = [], splits = [], clips = [], unmappedMates = [];
  for (const mol of molecules.values()) {
    // a mapped read whose mate is unmapped (placed beside it): the mate carries sequence the reference lacks
    for (const u of mol.unmapped) {
      const mate = [...mol.r1, ...mol.r2].find((a) => !a.supplementary && a.readNumber !== u.record.readNumber);
      if (mate && mate.wi >= 0) unmappedMates.push({ kind: 'unmapped', mol, aln: mate, wi: mate.wi, pos: mate.strand === '+' ? mate.end : mate.pos,
        dk: mol.dk || u.dk, mapq: mate.mapq, name: mol.name, event: 'Unmapped mate' });
    }
    // pairs: the two primary alignments
    const p1 = mol.r1.find((a) => !a.supplementary), p2 = mol.r2.find((a) => !a.supplementary);
    if (p1 && p2) {
      const pair = makePair(mol, p1, p2, windows);
      if (pair) {
        pairs.push(pair);
        if (!p1.record.isProperPair) {
          if (p1.wi >= 0) stats[p1.wi].discordant++;
          if (p2.wi >= 0 && p2.wi !== p1.wi) stats[p2.wi].discordant++;
        }
      }
    } else if (p1 || p2) {
      const p = p1 || p2;
      const r = p.record;
      if (r.isPaired && !r.isMateUnmapped && r.nextChrom !== null) {
        // the mate is not in the bundle (outside every window and not fetched): draw to its position
        const pair = makePair(mol, p, mateStub(p, windows), windows);
        if (pair) {
          pairs.push(pair);
          if (!r.isProperPair && p.wi >= 0) stats[p.wi].discordant++;
        }
      }
    }
    // splits: a primary with supplementary alignments (in the file, or only in its SA tag)
    for (const list of [mol.r1, mol.r2]) {
      const primary = list.find((a) => !a.supplementary);
      if (!primary) continue;
      const parts = [primary];
      const sa = primary.record.supplementaryAlignments;
      for (const s of sa) {
        const found = list.find((a) => a.supplementary && a.chrom === s.chrom && a.pos === s.pos);
        if (found) parts.push(found);
        else parts.push({ record: null, wi: windowOf(windows, s.chrom, s.pos), chrom: s.chrom, pos: s.pos,
          end: s.pos + cigarReferenceLength(s.cigar), strand: s.strand, mapq: s.mapq, cigar: s.cigar,
          clips: null, supplementary: true, dk: primary.dk, readNumber: primary.readNumber, fromSA: true });
      }
      if (parts.length > 1) {
        const split = makeSplit(mol, parts, windows);
        splits.push(split);
        if (primary.wi >= 0) stats[primary.wi].split++;
      }
      // clips at a breakpoint
      const [left, right] = primary.clips;
      if (left >= CLIP_MIN && primary.wi >= 0) {
        clips.push({ mol, aln: primary, wi: primary.wi, pos: primary.pos, side: 'left', len: left, dk: primary.dk });
        stats[primary.wi].clipped++;
      }
      if (right >= CLIP_MIN) {
        const wi = windowOf(windows, primary.chrom, primary.end - 1);
        if (wi >= 0) {
          clips.push({ mol, aln: primary, wi, pos: primary.end, side: 'right', len: right, dk: primary.dk });
          stats[wi].clipped++;
        }
      }
    }
  }
  let maxCoverage = 0;
  for (const cov of coverage) {
    for (let i = 0; i < cov.high.length; i++) maxCoverage = Math.max(maxCoverage, cov.high[i] + cov.low[i]);
  }
  return { windows, coverage, maxCoverage, molecules, pairs, splits, clips, unmappedMates, stats };
}

function getMolecule(molecules, name) {
  let mol = molecules.get(name);
  if (!mol) {
    mol = { name, r1: [], r2: [], unmapped: [], dk: false, mapq: null, taxon: null };
    molecules.set(name, mol);
  }
  return mol;
}

function mateStub(p, windows) {
  const r = p.record;
  const mc = r.tags.MC;
  const len = mc ? cigarReferenceLength(mc) : (r.end - r.pos);
  return { record: null, wi: windowOf(windows, r.nextChrom, r.nextPos), chrom: r.nextChrom, pos: r.nextPos,
    end: r.nextPos + len, strand: r.isMateReverse ? '-' : '+', mapq: r.tags.MQ ?? null, cigar: mc || '*', clips: null,
    supplementary: false, dk: false, readNumber: p.readNumber === 1 ? 2 : 1, stub: true };
}

function makePair(mol, a, b, windows) {
  const sameChrom = a.chrom === b.chrom;
  let first = a, second = b;
  if (sameChrom ? a.pos > b.pos : a.chrom > b.chrom) { first = b; second = a; }
  const ends = [first, second];
  let event, insert;
  if (!sameChrom) {
    event = EVENTS.BND;
    insert = null;
  } else {
    event = eventType(first.strand === '+', second.strand === '+');
    const tlen = Math.abs(a.record ? a.record.tlen : (b.record ? b.record.tlen : 0));
    insert = tlen || (Math.max(first.end, second.end) - Math.min(first.pos, second.pos));
  }
  const proper = (a.record && a.record.isProperPair) || (b.record && b.record.isProperPair) || false;
  const cross = first.wi !== second.wi;                       // an end in another window, or outside every window
  const linked = first.wi >= 0 && second.wi >= 0 && cross;     // both ends drawn, in two windows
  if (first.wi < 0 && second.wi < 0) return null;
  return { kind: 'pair', mol, ends, event, insert, proper, cross, linked, dk: mol.dk, mapq: Math.min(
    a.mapq ?? 255, b.mapq ?? 255), name: mol.name };
}

function makeSplit(mol, parts, windows) {
  parts.sort((x, y) => (x.chrom === y.chrom ? x.pos - y.pos : (x.chrom < y.chrom ? -1 : 1)));
  const first = parts[0], last = parts[parts.length - 1];
  const sameChrom = parts.every((p) => p.chrom === first.chrom);
  let event, span;
  if (!sameChrom) { event = EVENTS.BND; span = null; } else {
    event = eventType(first.strand === '+', last.strand === '+');
    span = last.end - first.pos;
  }
  const cross = parts.some((p) => p.wi !== first.wi);
  const placed = new Set(parts.map((p) => p.wi).filter((wi) => wi >= 0));
  const linked = placed.size > 1;
  // which side of each part is clipped says where the junction is on it
  for (const p of parts) p.clipSide = clippedSide(p.cigar);
  return { kind: 'split', mol, parts, event, insert: span, cross, linked, dk: mol.dk, mapq: Math.max(...parts.map((p) => p.mapq ?? 0)),
    name: mol.name };
}

/** A few numbers per window, for the panel's caption: reads, discordant pairs, split reads, clipped reads,
 *  informative reads. */
export function captionOf(stats) {
  return `${stats.reads} reads: ${stats.discordant} discordant, ${stats.split} split, ${stats.clipped} clipped` +
    (stats.informative ? `, ${stats.informative} informative` : '');
}

/** The breakpoint that a region's clipped reads agree on: the most frequent clip position (within 5 bp),
 *  as the tool's breakpoint_reads counts it. Returns { pos, n, side } or null. */
export function consensusBreakpoint(clips, wi) {
  const byPos = new Map();
  for (const c of clips) {
    if (c.wi !== wi) continue;
    const bin = Math.round(c.pos / 5) * 5;
    byPos.set(bin, (byPos.get(bin) || 0) + 1);
  }
  let best = null;
  for (const [bin, n] of byPos) if (!best || n > best.n) best = { pos: bin, n };
  return best && best.n >= 2 ? best : null;
}
