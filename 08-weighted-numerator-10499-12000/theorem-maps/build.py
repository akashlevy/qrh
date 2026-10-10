"""Compile isolated new modules against the original pinned QRH dependency cache."""
from pathlib import Path
import argparse, json, os, subprocess, sys

ROOT = Path(__file__).resolve().parent
CONFIG = json.loads((ROOT / 'lean_environment.json').read_text())
ENV = os.environ.copy()
ENV['LEAN_PATH'] = CONFIG['LEAN_PATH']

def compile_module(module):
    source = (ROOT / module.replace('.', '/')).with_suffix('.lean')
    output = source.with_suffix('.olean')
    log = ROOT / 'logs' / (module + '.log')
    print('Compiling', module, flush=True)
    with log.open('w') as f:
        result = subprocess.run([CONFIG['lean'], '-DautoImplicit=false', '-M', '8192',
                                 '-j', '2', '-R', str(ROOT), '-o', str(output), str(source)],
                                cwd=ROOT, env=ENV, stdout=f, stderr=subprocess.STDOUT)
    print(log.read_text(), end='')
    if result.returncode:
        sys.exit(result.returncode)
    print('PASS', module, flush=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('modules', nargs='+')
    args = parser.parse_args()
    for module in args.modules:
        compile_module(module)
