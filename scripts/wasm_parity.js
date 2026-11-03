// Check that the WebAssembly build writes the same files as the native CLI.
//
//   node scripts/wasm_parity.js <wasm-bindgen nodejs pkg dir> <fasterp binary> <R1.fq> <R2.fq>
//
// Each case runs natively in a temp dir and through `run` from the wasm
// package, then compares every FASTQ output (decompressed) and the JSON report
// without the fields that name the command or version.
const { execFileSync } = require('child_process');
const fs = require('fs');
const os = require('os');
const path = require('path');
const zlib = require('zlib');

const [pkg, fasterp, r1, r2] = process.argv.slice(2).map((p) => path.resolve(p));
const wasm = require(path.join(pkg, 'fasterp.js'));

const se = ['-i', 'in1.fq'];
const pe = ['-i', 'in1.fq', '-I', 'in2.fq'];
const cases = {
  'single-end': [...se, '-o', 'out.fq'],
  'single-end-trim': [...se, '-x', '-y', '-5', '-3', '-l', '30', '--failed_out', 'failed.fq', '-o', 'out.fq.gz'],
  'paired-end': [...pe, '-o', 'out1.fq', '-O', 'out2.fq'],
  'paired-correct-merge': [...pe, '-c', '-m', '--merged_out', 'merged.fq.gz', '-o', 'out1.fq.gz', '-O', 'out2.fq.gz'],
  'paired-dedup-overrep': [...pe, '-D', '-p', '-o', 'out1.fq', '-O', 'out2.fq'],
};

const inputs = new Map([
  ['in1.fq', new Uint8Array(fs.readFileSync(r1))],
  ['in2.fq', new Uint8Array(fs.readFileSync(r2))],
]);

function report(bytes) {
  const json = JSON.parse(Buffer.from(bytes).toString());
  delete json.command;
  delete json.summary.fasterp_version;
  return JSON.stringify(json);
}

function contents(name, bytes) {
  if (name.endsWith('.json')) return report(bytes);
  return name.endsWith('.gz') ? zlib.gunzipSync(bytes) : Buffer.from(bytes);
}

let failed = 0;
for (const [name, args] of Object.entries(cases)) {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'wasm-parity-'));
  fs.copyFileSync(r1, path.join(dir, 'in1.fq'));
  fs.copyFileSync(r2, path.join(dir, 'in2.fq'));
  execFileSync(fasterp, ['-w', '1', ...args], { cwd: dir, stdio: 'ignore' });
  const outputs = wasm.run(args, inputs);
  const expected = fs.readdirSync(dir).filter((f) => !f.startsWith('in') && !f.endsWith('.html'));
  const problems = expected.filter((f) => {
    if (!outputs.has(f)) return true;
    return !Buffer.from(contents(f, outputs.get(f))).equals(Buffer.from(contents(f, fs.readFileSync(path.join(dir, f)))));
  });
  console.log(`${problems.length ? 'FAIL' : 'ok  '} ${name} (${expected.join(', ')})${problems.length ? ': ' + problems.join(', ') : ''}`);
  failed += problems.length ? 1 : 0;
  fs.rmSync(dir, { recursive: true });
}
process.exit(failed ? 1 : 0);
