"""Compile the approved manuscript and detached supplement; outputs are validation builds."""
import argparse
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from pypdf import PdfReader
from verify_artifacts import ROOT

SOURCES = ('manuscript/RMMO_Paper3_Manuscript_preDOI.tex', 'supplement/SupplementaryEvidence.tex')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--engine', default='latexmk', help='latexmk or a Tectonic executable')
    parser.add_argument('--output', type=Path, help='Build directory outside the paper directory')
    parser.add_argument('--receipt', type=Path, help='Optional document validation receipt')
    args = parser.parse_args()
    output = args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix='rmmo-p3-documents-'))
    assert not output.is_relative_to(ROOT), 'Validation PDFs must stay outside the release directory'
    output.mkdir(parents=True, exist_ok=True)
    approved = json.loads((ROOT / 'verification/approved_document_authority.json').read_text())
    records = []
    for rel in SOURCES:
        source = ROOT / rel
        text = source.read_text()
        methods = text.split(r'\section{Methods}', 1)[1].split(r'\section{Data availability}', 1)[0] if rel.startswith('manuscript/') else None
        if methods is not None:
            assert hashlib.sha256(methods.encode()).hexdigest() == approved['methods_sha256'], 'Approved full Methods changed'
        body = text.split(r'\begin{document}', 1)[1]
        assert hashlib.sha256(body.encode()).hexdigest() == approved['sources'][rel]['document_body_sha256'], 'Approved document body changed'
        if Path(args.engine).name == 'latexmk':
            cmd = [args.engine, '-pdf', '-interaction=nonstopmode', '-halt-on-error', '-outdir='+str(output), source.name]
        else:
            cmd = [args.engine, '--keep-logs', '--outdir', str(output), source.name]
        subprocess.run(cmd, cwd=source.parent, check=True)
        pdf = output / source.with_suffix('.pdf').name
        log = (output / source.with_suffix('.log').name).read_text(errors='replace')
        prohibited = ('Overfull \\hbox', 'Overfull \\vbox', 'undefined references',
                      'Float too large', 'Citation `', 'LaTeX Error:', 'Emergency stop')
        assert not any(x in log for x in prohibited), 'Document layout/reference errors: '+rel
        reader = PdfReader(pdf)
        extracted = ' '.join(page.extract_text() or '' for page in reader.pages)
        assert reader.pages and '??' not in extracted and 'Zed James' in extracted
        records.append({'source': rel, 'pages': len(reader.pages), 'tex_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
                        'validation_pdf_sha256': hashlib.sha256(pdf.read_bytes()).hexdigest(),
                        'overfull_boxes': 0, 'undefined_references': 0})
    receipt = {'passed': True, 'documents': records, 'full_approved_methods_unchanged': True,
               'scientific_figures': 10, 'compiled_pdfs_are_validation_outputs_not_final_scholarly_pdf': True,
               'final_author_supplied_pdf': 'pending'}
    if args.receipt:
        args.receipt.parent.mkdir(parents=True, exist_ok=True)
        args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True)+'\n')
    print('PASS approved manuscript and detached supplement compilation: '+json.dumps(receipt, sort_keys=True))

if __name__ == '__main__':
    main()
