"""Compile the self-contained Chinese LaTeX paper using Tectonic."""
import argparse
from pathlib import Path
import re
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--engine', default='tectonic')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    output = root / 'build' / 'paper'
    output.mkdir(parents=True, exist_ok=True)
    source = root / 'paper' / 'su2_cocycle.tex'
    pdf = output / 'su2_cocycle.pdf'
    pdf.unlink(missing_ok=True)
    command = [args.engine, '--keep-logs', '--outdir', str(output), str(source)]
    subprocess.run(command, cwd=root, check=True)
    if not pdf.is_file():
        raise SystemExit('LaTeX returned without producing the expected PDF')
    log = (output / 'su2_cocycle.log').read_text(errors='replace')
    failures = [pattern for pattern in (
        r'Missing character:', r'There were undefined references',
        r'There were undefined citations', r'Overfull \\[hv]box',
    ) if re.search(pattern, log)]
    if failures:
        raise SystemExit('PDF quality checks failed; inspect build/paper/su2_cocycle.log: '
                         + ', '.join(failures))
    destination = root / 'paper' / pdf.name
    shutil.copy2(pdf, destination)
    print(f'Built {destination.relative_to(root)} ({destination.stat().st_size} bytes)')
    print('Compilation log: build/paper/su2_cocycle.log')


if __name__ == '__main__':
    main()
