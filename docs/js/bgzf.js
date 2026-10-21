// bgzf.js - BGZF blocks (the gzip variant BAM and its friends are written in), inflated with the browser's
// own DecompressionStream: no library, no build step. A BGZF file is a series of gzip members of at most
// 64 KiB each; a BAM virtual offset is (compressed offset of the block << 16) | (offset in the inflated block).

export const MAX_BLOCK = 65536 + 64; // the largest compressed block plus its header and trailer

/** The compressed size of the block starting at `pos` in `bytes` (its header's BSIZE + 1), or null when
 *  `bytes` does not hold a whole header there. */
export function blockSize(bytes, pos) {
  if (pos + 18 > bytes.length) return null;
  if (bytes[pos] !== 31 || bytes[pos + 1] !== 139 || bytes[pos + 2] !== 8) {
    throw new Error(`not a BGZF block at compressed offset ${pos}`);
  }
  const flg = bytes[pos + 3];
  if (!(flg & 4)) throw new Error(`BGZF block without an extra field at offset ${pos}`);
  const xlen = bytes[pos + 10] | (bytes[pos + 11] << 8);
  let p = pos + 12;
  const stop = p + xlen;
  while (p + 4 <= stop) {
    const si1 = bytes[p], si2 = bytes[p + 1];
    const slen = bytes[p + 2] | (bytes[p + 3] << 8);
    if (si1 === 66 && si2 === 67 && slen === 2) {
      return (bytes[p + 4] | (bytes[p + 5] << 8)) + 1;
    }
    p += 4 + slen;
  }
  throw new Error(`BGZF block without a BC subfield at offset ${pos}`);
}

function headerLength(bytes, pos) {
  return 12 + (bytes[pos + 10] | (bytes[pos + 11] << 8));
}

/** Inflate raw deflate data. */
export async function inflateRaw(data) {
  const stream = new Blob([data]).stream().pipeThrough(new DecompressionStream('deflate-raw'));
  return new Uint8Array(await new Response(stream).arrayBuffer());
}

/** Inflate the one block that starts at `pos`; returns { size, data } (size: compressed bytes). */
export async function inflateBlock(bytes, pos) {
  const size = blockSize(bytes, pos);
  if (size === null || pos + size > bytes.length) return null;
  const hlen = headerLength(bytes, pos);
  const cdata = bytes.subarray(pos + hlen, pos + size - 8);
  const isize = bytes[pos + size - 4] | (bytes[pos + size - 3] << 8) | (bytes[pos + size - 2] << 16) |
    (bytes[pos + size - 1] << 24);
  const data = isize === 0 ? new Uint8Array(0) : await inflateRaw(cdata);
  if (data.length !== (isize >>> 0)) {
    throw new Error(`BGZF block at ${pos}: inflated ${data.length} bytes, the trailer says ${isize >>> 0}`);
  }
  return { size, data };
}

/** Inflate every whole block in `bytes` whose compressed offset (relative to `base`) is at most `lastOffset`.
 *  Returns [{ coffset, data }] in file order. */
export async function inflateBlocks(bytes, base, lastOffset = Infinity) {
  const out = [];
  let pos = 0;
  while (pos < bytes.length && base + pos <= lastOffset) {
    const block = await inflateBlock(bytes, pos);
    if (block === null) break;
    out.push({ coffset: base + pos, data: block.data });
    pos += block.size;
  }
  return out;
}

export function concat(chunks) {
  let n = 0;
  for (const c of chunks) n += c.length;
  const out = new Uint8Array(n);
  let p = 0;
  for (const c of chunks) { out.set(c, p); p += c.length; }
  return out;
}
