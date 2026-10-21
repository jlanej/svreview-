// bam.js - a small BAM reader for the browser: the header, and the records of a region through the .bai,
// over HTTP range requests (a static site, GitHub Pages included), a local File, or bytes in memory.
// It reads what a review needs - name, flag, position, MAPQ, CIGAR, mate, template length, sequence and
// the tags - and nothing else; the format is the SAM specification's section 4.

import { inflateBlocks, blockSize, concat, MAX_BLOCK } from './bgzf.js';
import { parseIndex, chunksForRegion, coffsetOf, uoffsetOf } from './bai.js';

export const CIGAR_OPS = 'MIDNSHP=X';
const SEQ_NT = '=ACMGRSVTWYHKDBN';
export const FLAG = {
  PAIRED: 0x1, PROPER: 0x2, UNMAPPED: 0x4, MATE_UNMAPPED: 0x8, REVERSE: 0x10, MATE_REVERSE: 0x20,
  READ1: 0x40, READ2: 0x80, SECONDARY: 0x100, QCFAIL: 0x200, DUPLICATE: 0x400, SUPPLEMENTARY: 0x800,
};

// ---- byte sources ---------------------------------------------------------------------------------------

/** A URL read by range requests. A server that answers 200 to a range request (no range support) hands over
 *  the whole file once, which is then sliced in memory: fine for the mini-BAMs of a bundle. */
export class RemoteSource {
  constructor(url, fetchFn) {
    this.url = url;
    this.size = null;
    this.whole = null;
    this.fetch = fetchFn || ((...a) => fetch(...a));
  }
  async read(offset, length) {
    if (this.whole) return this.whole.subarray(offset, offset + length);
    const res = await this.fetch(this.url, { headers: { Range: `bytes=${offset}-${offset + length - 1}` } });
    if (res.status === 206) {
      const range = res.headers.get('Content-Range');
      const m = range && /\/(\d+)$/.exec(range);
      if (m) this.size = parseInt(m[1], 10);
      return new Uint8Array(await res.arrayBuffer());
    }
    if (res.status === 200) {
      this.whole = new Uint8Array(await res.arrayBuffer());
      this.size = this.whole.length;
      return this.whole.subarray(offset, offset + length);
    }
    throw new Error(`${this.url}: HTTP ${res.status}`);
  }
  async readAll() {
    if (!this.whole) {
      const res = await this.fetch(this.url);
      if (!res.ok) throw new Error(`${this.url}: HTTP ${res.status}`);
      this.whole = new Uint8Array(await res.arrayBuffer());
      this.size = this.whole.length;
    }
    return this.whole;
  }
}

/** A File from an <input type="file"> or a drop. */
export class FileSource {
  constructor(file) { this.file = file; this.size = file.size; }
  async read(offset, length) {
    return new Uint8Array(await this.file.slice(offset, offset + length).arrayBuffer());
  }
  async readAll() { return new Uint8Array(await this.file.arrayBuffer()); }
}

/** Bytes already in memory (tests). */
export class BytesSource {
  constructor(bytes) { this.bytes = bytes; this.size = bytes.length; }
  async read(offset, length) { return this.bytes.subarray(offset, offset + length); }
  async readAll() { return this.bytes; }
}

// ---- records --------------------------------------------------------------------------------------------

export class BamRecord {
  constructor(fields) { Object.assign(this, fields); }
  get isPaired() { return (this.flag & FLAG.PAIRED) !== 0; }
  get isProperPair() { return (this.flag & FLAG.PROPER) !== 0; }
  get isUnmapped() { return (this.flag & FLAG.UNMAPPED) !== 0; }
  get isMateUnmapped() { return (this.flag & FLAG.MATE_UNMAPPED) !== 0; }
  get isReverse() { return (this.flag & FLAG.REVERSE) !== 0; }
  get isMateReverse() { return (this.flag & FLAG.MATE_REVERSE) !== 0; }
  get isRead1() { return (this.flag & FLAG.READ1) !== 0; }
  get isRead2() { return (this.flag & FLAG.READ2) !== 0; }
  get isSecondary() { return (this.flag & FLAG.SECONDARY) !== 0; }
  get isQCFail() { return (this.flag & FLAG.QCFAIL) !== 0; }
  get isDuplicate() { return (this.flag & FLAG.DUPLICATE) !== 0; }
  get isSupplementary() { return (this.flag & FLAG.SUPPLEMENTARY) !== 0; }
  /** 1 for read 1, 2 for read 2, 0 for an unpaired read: with the name, the molecule's end. */
  get readNumber() { return this.isRead1 ? 1 : (this.isRead2 ? 2 : 0); }
  get strand() { return this.isReverse ? '-' : '+'; }
  /** The CIGAR as a string ('*' when there is none). */
  get cigarString() {
    if (!this.cigar.length) return '*';
    let s = '';
    for (const [len, op] of this.cigar) s += len + op;
    return s;
  }
  /** Soft (and hard) clipping on each side of the alignment, in read orientation as aligned: [left, right]. */
  get clips() {
    const c = this.cigar;
    if (!c.length) return [0, 0];
    let left = 0, right = 0, i = 0, j = c.length - 1;
    while (i < c.length && (c[i][1] === 'S' || c[i][1] === 'H')) { left += c[i][0]; i++; }
    while (j >= 0 && (c[j][1] === 'S' || c[j][1] === 'H')) { right += c[j][0]; j--; }
    return [left, right];
  }
  /** Soft clips only (bases present in seq): [left, right]. */
  get softClips() {
    const c = this.cigar;
    if (!c.length) return [0, 0];
    const left = c[0][1] === 'S' ? c[0][0] : 0;
    const right = c.length > 1 && c[c.length - 1][1] === 'S' ? c[c.length - 1][0] : 0;
    return [left, right];
  }
  /** The reference-aligned blocks of the CIGAR, as [start, end) pairs (deletions and skips break them). */
  get blocks() {
    const out = [];
    let pos = this.pos, cur = null;
    for (const [len, op] of this.cigar) {
      if (op === 'M' || op === '=' || op === 'X') {
        if (cur && cur[1] === pos) cur[1] = pos + len; else { cur = [pos, pos + len]; out.push(cur); }
        pos += len;
      } else if (op === 'D' || op === 'N') {
        pos += len;
        cur = null;
      }
    }
    return out;
  }
  /** The SA tag's entries: [{ chrom, pos (0-based), strand, cigar, mapq, nm }]. */
  get supplementaryAlignments() {
    const sa = this.tags.SA;
    if (!sa) return [];
    const out = [];
    for (const entry of sa.split(';')) {
      if (!entry) continue;
      const f = entry.split(',');
      if (f.length < 6) continue;
      out.push({ chrom: f[0], pos: parseInt(f[1], 10) - 1, strand: f[2], cigar: f[3], mapq: parseInt(f[4], 10),
        nm: parseInt(f[5], 10) });
    }
    return out;
  }
}

/** The reference span of a CIGAR string (M, D, N, =, X). */
export function cigarReferenceLength(cigarString) {
  let n = 0;
  const re = /(\d+)([MIDNSHP=X])/g;
  let m;
  while ((m = re.exec(cigarString)) !== null) {
    if ('MDN=X'.includes(m[2])) n += parseInt(m[1], 10);
  }
  return n;
}

/** Which side of a split alignment is clipped, from its CIGAR string: 'left', 'right', 'both' or 'none'. */
export function clippedSide(cigarString) {
  const left = /^\d+[SH]/.test(cigarString), right = /\d+[SH]$/.test(cigarString);
  return left && right ? 'both' : (left ? 'left' : (right ? 'right' : 'none'));
}

const decoder = new TextDecoder('utf-8');

function readTags(bytes, dv, p, end) {
  const tags = {};
  while (p + 3 <= end) {
    const name = String.fromCharCode(bytes[p], bytes[p + 1]);
    const type = String.fromCharCode(bytes[p + 2]);
    p += 3;
    let v;
    switch (type) {
      case 'A': v = String.fromCharCode(bytes[p]); p += 1; break;
      case 'c': v = dv.getInt8(p); p += 1; break;
      case 'C': v = bytes[p]; p += 1; break;
      case 's': v = dv.getInt16(p, true); p += 2; break;
      case 'S': v = dv.getUint16(p, true); p += 2; break;
      case 'i': v = dv.getInt32(p, true); p += 4; break;
      case 'I': v = dv.getUint32(p, true); p += 4; break;
      case 'f': v = dv.getFloat32(p, true); p += 4; break;
      case 'Z': case 'H': {
        let q = p;
        while (q < end && bytes[q] !== 0) q++;
        v = decoder.decode(bytes.subarray(p, q));
        p = q + 1;
        break;
      }
      case 'B': {
        const sub = String.fromCharCode(bytes[p]);
        const n = dv.getInt32(p + 1, true);
        p += 5;
        const size = { c: 1, C: 1, s: 2, S: 2, i: 4, I: 4, f: 4 }[sub];
        const arr = new Array(n);
        for (let i = 0; i < n; i++) {
          switch (sub) {
            case 'c': arr[i] = dv.getInt8(p); break;
            case 'C': arr[i] = bytes[p]; break;
            case 's': arr[i] = dv.getInt16(p, true); break;
            case 'S': arr[i] = dv.getUint16(p, true); break;
            case 'i': arr[i] = dv.getInt32(p, true); break;
            case 'I': arr[i] = dv.getUint32(p, true); break;
            default: arr[i] = dv.getFloat32(p, true);
          }
          p += size;
        }
        v = arr;
        break;
      }
      default:
        throw new Error(`unknown tag type ${type} for ${name}`);
    }
    tags[name] = v;
  }
  return tags;
}

/** Parse one record at `p` (the start of its block_size field); returns { record, next } or null when the
 *  bytes end before the record does. */
export function parseRecord(bytes, dv, p, refNames) {
  if (p + 4 > bytes.length) return null;
  const blockSize = dv.getInt32(p, true);
  const end = p + 4 + blockSize;
  if (end > bytes.length) return null;
  let q = p + 4;
  const refId = dv.getInt32(q, true); q += 4;
  const pos = dv.getInt32(q, true); q += 4;
  const lReadName = bytes[q]; q += 1;
  const mapq = bytes[q]; q += 1;
  q += 2; // bin
  const nCigar = dv.getUint16(q, true); q += 2;
  const flag = dv.getUint16(q, true); q += 2;
  const lSeq = dv.getInt32(q, true); q += 4;
  const nextRefId = dv.getInt32(q, true); q += 4;
  const nextPos = dv.getInt32(q, true); q += 4;
  const tlen = dv.getInt32(q, true); q += 4;
  const name = decoder.decode(bytes.subarray(q, q + lReadName - 1)); q += lReadName;
  const cigar = new Array(nCigar);
  let refLen = 0;
  for (let i = 0; i < nCigar; i++) {
    const v = dv.getUint32(q, true); q += 4;
    const len = v >>> 4, op = CIGAR_OPS[v & 0xf];
    cigar[i] = [len, op];
    if (op === 'M' || op === 'D' || op === 'N' || op === '=' || op === 'X') refLen += len;
  }
  const seqOffset = q;
  q += (lSeq + 1) >> 1;
  const qualOffset = q;
  q += lSeq;
  const tags = readTags(bytes, dv, q, end);
  const unmapped = (flag & FLAG.UNMAPPED) !== 0;
  const record = new BamRecord({
    name, flag, refId, pos, mapq, cigar, tlen, nextRefId, nextPos, lSeq, tags,
    chrom: refId >= 0 ? refNames[refId] : null,
    nextChrom: nextRefId >= 0 ? refNames[nextRefId] : null,
    end: unmapped || refLen === 0 ? pos + 1 : pos + refLen,
    _bytes: bytes, _seqOffset: seqOffset, _qualOffset: qualOffset,
  });
  return { record, next: end };
}

/** The read's sequence, decoded on demand. */
export function sequenceOf(r) {
  const n = r.lSeq, b = r._bytes, o = r._seqOffset;
  let s = '';
  for (let i = 0; i < n; i++) {
    const byte = b[o + (i >> 1)];
    s += SEQ_NT[(i & 1) ? (byte & 0xf) : (byte >> 4)];
  }
  return s;
}

/** Base qualities as numbers (Phred), or null when absent (0xff). */
export function qualitiesOf(r) {
  const n = r.lSeq, b = r._bytes, o = r._qualOffset;
  if (n > 0 && b[o] === 0xff) return null;
  const out = new Array(n);
  for (let i = 0; i < n; i++) out[i] = b[o + i];
  return out;
}

// ---- the file -------------------------------------------------------------------------------------------

export class BamFile {
  /** @param source a RemoteSource / FileSource / BytesSource for the BAM, `indexSource` the same for its .bai */
  constructor(source, indexSource) {
    this.source = source;
    this.indexSource = indexSource;
    this._header = null;
    this._index = null;
  }

  /** { text, refs: [{ name, length }], refNames: [...] } */
  async header() {
    if (this._header) return this._header;
    let compressed = await this.source.read(0, 65536 * 2);
    let blocks = await inflateBlocks(compressed, 0);
    let data = concat(blocks.map((b) => b.data));
    let parsed = this._parseHeader(data);
    let got = compressed.length;
    while (parsed === null) {
      const more = await this.source.read(got, 65536 * 4);
      if (!more.length) throw new Error('BAM header is truncated');
      compressed = concat([compressed, more]);
      got += more.length;
      blocks = await inflateBlocks(compressed, 0);
      data = concat(blocks.map((b) => b.data));
      parsed = this._parseHeader(data);
    }
    this._header = parsed;
    return parsed;
  }

  _parseHeader(data) {
    if (data.length < 12) return null;
    if (data[0] !== 66 || data[1] !== 65 || data[2] !== 77 || data[3] !== 1) throw new Error('not a BAM file (no BAM\\1 magic)');
    const dv = new DataView(data.buffer, data.byteOffset, data.byteLength);
    const lText = dv.getInt32(4, true);
    let p = 8 + lText;
    if (p + 4 > data.length) return null;
    const text = decoder.decode(data.subarray(8, 8 + lText)).replace(/\0+$/, '');
    const nRef = dv.getInt32(p, true); p += 4;
    const refs = [];
    for (let i = 0; i < nRef; i++) {
      if (p + 4 > data.length) return null;
      const lName = dv.getInt32(p, true); p += 4;
      if (p + lName + 4 > data.length) return null;
      const name = decoder.decode(data.subarray(p, p + lName - 1)); p += lName;
      const length = dv.getInt32(p, true); p += 4;
      refs.push({ name, length });
    }
    return { text, refs, refNames: refs.map((r) => r.name) };
  }

  async index() {
    if (!this._index) {
      const bytes = await this.indexSource.readAll();
      this._index = parseIndex(bytes);
    }
    return this._index;
  }

  /** The records overlapping [start, end) of `chrom` (0-based, half-open), in file order; a placed unmapped
   *  read counts at its position. Secondary, QC-failed and duplicate records are kept: the caller decides. */
  async records(chrom, start, end) {
    const header = await this.header();
    const refId = header.refNames.indexOf(chrom);
    if (refId < 0) return [];
    const index = await this.index();
    const chunks = chunksForRegion(index, refId, start, end);
    const out = [];
    for (const [vbeg, vend] of chunks) {
      const cbeg = coffsetOf(vbeg), cend = coffsetOf(vend);
      let length = cend - cbeg + MAX_BLOCK;
      if (this.source.size !== null && this.source.size !== undefined) length = Math.min(length, this.source.size - cbeg);
      const compressed = await this.source.read(cbeg, length);
      const blocks = await inflateBlocks(compressed, cbeg, cend);
      if (!blocks.length) continue;
      // the inflated bytes of the chunk: from the first block's uoffset to the last block's
      const data = concat(blocks.map((b) => b.data));
      let stop = data.length;
      if (blocks[blocks.length - 1].coffset === cend) {
        let before = 0;
        for (let i = 0; i < blocks.length - 1; i++) before += blocks[i].data.length;
        stop = Math.min(stop, before + uoffsetOf(vend));
      }
      const dv = new DataView(data.buffer, data.byteOffset, data.byteLength);
      let p = uoffsetOf(vbeg);
      while (p < stop) {
        const parsed = parseRecord(data, dv, p, header.refNames);
        if (parsed === null) break;
        const r = parsed.record;
        p = parsed.next;
        if (r.refId !== refId) { if (r.refId > refId) break; continue; }
        if (r.pos >= end) break;
        if (r.end > start) out.push(r);
      }
    }
    return out;
  }
}

/** Open a BAM by URL (its .bai at `url + '.bai'` unless given). */
export function openRemote(url, indexUrl, fetchFn) {
  return new BamFile(new RemoteSource(url, fetchFn), new RemoteSource(indexUrl || url + '.bai', fetchFn));
}

/** Open a BAM from two File objects. */
export function openFiles(bamFile, indexFile) {
  return new BamFile(new FileSource(bamFile), new FileSource(indexFile));
}

/** Open a BAM from bytes in memory. */
export function openBytes(bamBytes, indexBytes) {
  return new BamFile(new BytesSource(bamBytes), new BytesSource(indexBytes));
}

/** Verify the file's BGZF structure from the start (used by tests): the number of blocks. */
export async function countBlocks(bytes) {
  let pos = 0, n = 0;
  while (pos < bytes.length) {
    const size = blockSize(bytes, pos);
    if (size === null) break;
    pos += size;
    n++;
  }
  return n;
}
