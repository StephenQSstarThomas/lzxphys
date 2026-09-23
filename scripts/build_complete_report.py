"""Combine the fully proved SU(2) note and the ADE supplement into one PDF."""
from pathlib import Path
import shutil
import subprocess


def main():
    root = Path(__file__).resolve().parents[1]
    main_pdf = root / "paper" / "su2_cocycle.pdf"
    ade_pdf = root / "paper" / "su2_ade_review.pdf"
    output = root / "paper" / "su2_complete_report.pdf"
    if not main_pdf.is_file() or not ade_pdf.is_file():
        raise SystemExit("Build paper/su2_cocycle.pdf and paper/su2_ade_review.pdf first")
    pdfunite = shutil.which("pdfunite")
    if not pdfunite:
        raise SystemExit("pdfunite is required to combine the two verified PDFs")
    subprocess.run([pdfunite, str(main_pdf), str(ade_pdf), str(output)], check=True)
    print(f"Wrote {output.relative_to(root)} ({output.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
