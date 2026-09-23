"""Compile one integrated report with one title, contents, and bibliography."""
import argparse
from pathlib import Path
import re
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", default="tectonic")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    source = root / "paper" / "su2_complete_report.tex"
    build = root / "build" / "complete_report"
    build.mkdir(parents=True, exist_ok=True)
    pdf = build / source.with_suffix(".pdf").name
    pdf.unlink(missing_ok=True)
    subprocess.run([args.engine, "--keep-logs", "--outdir", str(build), str(source)],
                   cwd=root, check=True)
    if not pdf.is_file():
        raise SystemExit("LaTeX returned without producing the complete report")
    log = (build / source.with_suffix(".log").name).read_text(errors="replace")
    failures = [pattern for pattern in (
        r"Missing character:", r"There were undefined references",
        r"There were undefined citations", r"Overfull \\[hv]box",
        r"multiply defined", r"destination with the same identifier",
    ) if re.search(pattern, log)]
    if failures:
        raise SystemExit("PDF quality checks failed; inspect " + str(build) + ": "
                         + ", ".join(failures))
    output = root / "paper" / pdf.name
    shutil.copy2(pdf, output)
    print(f"Built {output.relative_to(root)} ({output.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
