// node --test tests/js   (node >= 18: DecompressionStream, Blob and Response are global)
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';
import path from 'node:path';

import { runChecks } from './check_bam.mjs';

const here = path.dirname(fileURLToPath(import.meta.url));
const fixtures = path.join(here, 'fixtures');
const loadBytes = async (p) => new Uint8Array(await readFile(path.join(fixtures, p)));
const loadJson = async (p) => JSON.parse(await readFile(path.join(fixtures, p), 'utf8'));

test('the BAM reader and the read classification against pysam', async () => {
  const results = await runChecks(loadBytes, loadJson);
  for (const r of results) console.log(`${r.ok ? 'ok  ' : 'FAIL'} ${r.name}${r.detail ? ' - ' + r.detail : ''}`);
  const failed = results.filter((r) => !r.ok);
  assert.equal(failed.length, 0, failed.map((r) => `${r.name}: ${r.detail}`).join('\n'));
  assert.ok(results.length > 20);
});
