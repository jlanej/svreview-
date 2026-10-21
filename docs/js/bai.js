// bai.js - the BAM index: which compressed chunks of the file hold the records of a region.
// Virtual offsets are kept as plain numbers, coffset * 65536 + uoffset (exact below 128 GB of BAM).

const METADATA_BIN = 37450;

export function voffset(coffset, uoffset) { return coffset * 65536 + uoffset; }
export function coffsetOf(v) { return Math.floor(v / 65536); }
export function uoffsetOf(v) { return v % 65536; }

function readU64(dv, p) {
  const lo = dv.getUint32(p, true), hi = dv.getUint32(p + 4, true);
  return hi * 4294967296 + lo; // < 2^53 for any offset a BAM index holds
}

/** Parse a .bai file: { refs: [{ bins: Map(bin -> [[beg, end]]), intervals: [voffset], meta }] }. */
export function parseIndex(bytes) {
  const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  if (bytes[0] !== 66 || bytes[1] !== 65 || bytes[2] !== 73 || bytes[3] !== 1) {
    throw new Error('not a BAM index (no BAI\\1 magic)');
  }
  let p = 4;
  const nRef = dv.getInt32(p, true); p += 4;
  const refs = [];
  for (let r = 0; r < nRef; r++) {
    const nBin = dv.getInt32(p, true); p += 4;
    const bins = new Map();
    let meta = null;
    for (let b = 0; b < nBin; b++) {
      const bin = dv.getUint32(p, true); p += 4;
      const nChunk = dv.getInt32(p, true); p += 4;
      const chunks = [];
      for (let c = 0; c < nChunk; c++) {
        chunks.push([readU64(dv, p), readU64(dv, p + 8)]);
        p += 16;
      }
      if (bin === METADATA_BIN) {
        meta = chunks;
      } else {
        bins.set(bin, chunks);
      }
    }
    const nIntv = dv.getInt32(p, true); p += 4;
    const intervals = new Array(nIntv);
    for (let i = 0; i < nIntv; i++) { intervals[i] = readU64(dv, p); p += 8; }
    refs.push({ bins, intervals, meta });
  }
  return { refs };
}

/** The bins overlapping [beg, end) (0-based, half-open): the five-level scheme of the SAM specification. */
export function reg2bins(beg, end) {
  const list = [0];
  end -= 1;
  for (let k = 1 + (beg >> 26); k <= 1 + (end >> 26); k++) list.push(k);
  for (let k = 9 + (beg >> 23); k <= 9 + (end >> 23); k++) list.push(k);
  for (let k = 73 + (beg >> 20); k <= 73 + (end >> 20); k++) list.push(k);
  for (let k = 585 + (beg >> 17); k <= 585 + (end >> 17); k++) list.push(k);
  for (let k = 4681 + (beg >> 14); k <= 4681 + (end >> 14); k++) list.push(k);
  return list;
}

/** The chunks (as [beg, end] virtual offsets, merged, in file order) to read for the records overlapping
 *  [beg, end) of reference `refId`; [] when the index has nothing there. */
export function chunksForRegion(index, refId, beg, end) {
  const ref = index.refs[refId];
  if (!ref || end <= beg) return [];
  const minOffset = ref.intervals.length ? ref.intervals[Math.min(beg >> 14, ref.intervals.length - 1)] : 0;
  const chunks = [];
  for (const bin of reg2bins(beg, end)) {
    const list = ref.bins.get(bin);
    if (!list) continue;
    for (const [cb, ce] of list) {
      if (ce > minOffset) chunks.push([Math.max(cb, minOffset), ce]);
    }
  }
  chunks.sort((a, b) => a[0] - b[0]);
  const merged = [];
  for (const c of chunks) {
    const last = merged[merged.length - 1];
    if (last && c[0] <= last[1]) {
      if (c[1] > last[1]) last[1] = c[1];
    } else {
      merged.push([c[0], c[1]]);
    }
  }
  return merged;
}
