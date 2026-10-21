// check_bam.mjs - the checks of the browser's BAM reader against pysam's reading of the same files (the
// fixtures), shared by the node test (bam.test.mjs) and the browser runner (browser.html).

import { BamFile, BytesSource, sequenceOf, countBlocks } from '../../docs/js/bam.js';
import { parseIndex, reg2bins, chunksForRegion } from '../../docs/js/bai.js';
import { buildMember, windowOf, eventType, EVENTS } from '../../docs/js/reads.js';

export const FIXTURES = ['chr1_54860_55115', 'chr1_9860_10000'];

function eq(a, b, what) {
  if (a !== b) throw new Error(`${what}: got ${JSON.stringify(a)}, expected ${JSON.stringify(b)}`);
}

/**
 * @param loadBytes async (relativePath) => Uint8Array  (relative to tests/js/fixtures)
 * @param loadJson  async (relativePath) => object
 * @returns [{ name, ok, detail }]
 */
export async function runChecks(loadBytes, loadJson) {
  const results = [];
  const check = async (name, fn) => {
    try {
      const detail = await fn();
      results.push({ name, ok: true, detail: detail || '' });
    } catch (e) {
      results.push({ name, ok: false, detail: e.message });
    }
  };

  await check('reg2bins: the five-level scheme', () => {
    const bins = reg2bins(0, 1);
    eq(bins.join(','), '0,1,9,73,585,4681', 'bins of [0,1)');
    const b2 = reg2bins(16384 * 3, 16384 * 3 + 1);
    eq(b2[b2.length - 1], 4681 + 3, 'leaf bin of 49152');
    return bins.join(',');
  });

  for (const fx of FIXTURES) {
    const expected = await loadJson(`${fx}.expected.json`);
    for (const role of ['child', 'father', 'mother']) {
      const member = expected.members[role];
      const bamBytes = await loadBytes(member.bam);
      const baiBytes = await loadBytes(member.bam + '.bai');
      await check(`${fx} ${role}: BGZF blocks`, async () => {
        const n = await countBlocks(bamBytes);
        if (n < 2) throw new Error(`only ${n} blocks`);
        return `${n} blocks`;
      });
      await check(`${fx} ${role}: header`, async () => {
        const bam = new BamFile(new BytesSource(bamBytes), new BytesSource(baiBytes));
        const h = await bam.header();
        eq(h.refs.length, member.refs.length, 'reference count');
        member.refs.forEach((r, i) => { eq(h.refs[i].name, r.name, `ref ${i} name`); eq(h.refs[i].length, r.length, `ref ${i} length`); });
        if (!h.text.includes('@HD')) throw new Error('no @HD line in the header text');
        return `${h.refs.length} refs`;
      });
      await check(`${fx} ${role}: index`, async () => {
        const index = parseIndex(baiBytes);
        eq(index.refs.length, member.refs.length, 'index refs');
        const w = member.windows[0];
        const refId = member.refs.findIndex((r) => r.name === w.chrom);
        const chunks = chunksForRegion(index, refId, w.start, w.end);
        if (!chunks.length) throw new Error('no chunks for the first window');
        return `${chunks.length} chunk(s)`;
      });
      for (const w of member.windows) {
        await check(`${fx} ${role}: records of ${w.chrom}:${w.start}-${w.end}`, async () => {
          const bam = new BamFile(new BytesSource(bamBytes), new BytesSource(baiBytes));
          const got = await bam.records(w.chrom, w.start, w.end);
          eq(got.length, w.reads.length, 'record count');
          for (let i = 0; i < got.length; i++) {
            const g = got[i], e = w.reads[i];
            const at = `record ${i} (${e.name})`;
            eq(g.name, e.name, `${at} name`);
            eq(g.flag, e.flag, `${at} flag`);
            eq(g.chrom, e.chrom, `${at} chrom`);
            eq(g.pos, e.pos, `${at} pos`);
            eq(g.end, e.end, `${at} end`);
            eq(g.mapq, e.mapq, `${at} mapq`);
            eq(g.cigarString, e.cigar, `${at} cigar`);
            eq(g.tlen, e.tlen, `${at} tlen`);
            eq(g.nextChrom, e.mate_chrom, `${at} mate chrom`);
            eq(g.nextPos, e.mate_pos, `${at} mate pos`);
            eq(g.lSeq, e.lseq, `${at} seq length`);
            if (e.seq) eq(sequenceOf(g), e.seq, `${at} sequence`);
            for (const [k, v] of Object.entries(e.tags)) eq(g.tags[k], v, `${at} tag ${k}`);
          }
          return `${got.length} records`;
        });
      }
    }
    await check(`${fx}: classification over the windows`, async () => {
      const windows = expected.windows;
      const child = expected.members.child;
      const bamBytes = await loadBytes(child.bam);
      const baiBytes = await loadBytes(child.bam + '.bai');
      const bam = new BamFile(new BytesSource(bamBytes), new BytesSource(baiBytes));
      const per = [];
      for (const w of windows) per.push(await bam.records(w.chrom, w.start, w.end));
      const m = buildMember(windows, per, {});
      if (m.maxCoverage < 10) throw new Error(`max coverage ${m.maxCoverage}: the demo is 30x`);
      const cross = m.pairs.filter((p) => p.linked).length + m.splits.filter((s) => s.linked).length;
      if (windows.length > 1 && cross === 0) throw new Error('no pair or split links the two windows');
      const splits = m.splits.length, discordant = m.pairs.filter((p) => !p.proper).length;
      if (splits < 2 || discordant < 2) throw new Error(`${splits} splits, ${discordant} discordant pairs: the junction should show`);
      if (fx === 'chr1_54860_55115') {
        if (!m.pairs.some((p) => p.event === EVENTS.BND)) throw new Error('no inter-chromosomal pair at the translocation');
      } else if (!m.pairs.some((p) => p.event === EVENTS.DEL && !p.proper && p.linked)) {
        throw new Error('no deletion-type discordant pair across the windows');
      }
      eq(windowOf(windows, windows[0].chrom, windows[0].start), 0, 'windowOf start');
      eq(windowOf(windows, windows[0].chrom, windows[0].end), windows.length > 1 && windows[1].chrom === windows[0].chrom && windows[1].start === windows[0].end ? 1 : -1, 'windowOf end');
      eq(eventType(true, false), EVENTS.NORMAL, 'event +-');
      eq(eventType(false, true), EVENTS.DUP, 'event -+');
      eq(eventType(true, true), EVENTS.INV, 'event ++');
      return `${m.pairs.length} pairs (${discordant} discordant), ${splits} splits, ${cross} linking the windows, ${m.clips.length} clips`;
    });
  }
  return results;
}
