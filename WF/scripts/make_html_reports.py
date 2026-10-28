"""
HTML reports of the analysis notebooks for collaborators: the executed notebooks without code (markdown, printed
tables and figures only), with the full-size figures next to them.

- Code cells are hidden (`--no-input`); cells tagged `report-skip` (setup, configuration, raw table dumps) are left
  out completely, including their output.
- Each figure saved by a notebook shows its file name with links (`figure_caption` in the notebooks). The links are
  relative (`figures/<path below FIGURES_DIR>`); the linked files are copied there, so the report folder can be
  zipped or uploaded as a whole.
- The notebooks are exported as they are: run them first (Kernel > Restart & Run All, or an sbatch job with
  `jupyter nbconvert --execute --inplace`).

Usage (from the project folder, sanity_utils environment):
    python WF/scripts/make_html_reports.py --out <folder> [notebook.ipynb ...]
"""

import argparse
import datetime
import html
import os
import re
import shutil

import nbformat
from dotenv import load_dotenv
from nbconvert import HTMLExporter
from nbconvert.preprocessors import TagRemovePreprocessor

PROJECT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
NOTEBOOKS = ["pbody_RNA_Spang.TotalDE.ipynb", "pbody_RNA_Spang.fracSanity.ipynb", "pbody_RNA_Spang.translation.ipynb"]


def export(nb_path, out_dir):
    nb = nbformat.read(nb_path, as_version=4)
    title = next((c.source.splitlines()[0].lstrip("# ").strip() for c in nb.cells
                  if c.cell_type == "markdown" and c.source.startswith("# ")), os.path.basename(nb_path))
    exporter = HTMLExporter(template_name="lab", exclude_input=True, exclude_input_prompt=True,
                            exclude_output_prompt=True)
    exporter.register_preprocessor(TagRemovePreprocessor(remove_cell_tags={"report-skip"}), enabled=True)
    body, _ = exporter.from_notebook_node(nb, resources={"metadata": {"name": title}})
    name = os.path.splitext(os.path.basename(nb_path))[0] + ".html"
    with open(os.path.join(out_dir, name), "w") as fh:
        fh.write(body)
    n_exec = sum(1 for c in nb.cells if c.cell_type == "code" and c.get("execution_count"))
    n_code = sum(1 for c in nb.cells if c.cell_type == "code" and c.source.strip())
    return name, title, body, n_exec, n_code


def copy_figures(body, figures_dir, out_dir):
    copied, missing = 0, []
    for rel in sorted(set(re.findall(r'href="figures/([^"]+)"', body))):
        src = os.path.join(figures_dir, html.unescape(rel))
        dst = os.path.join(out_dir, "figures", html.unescape(rel))
        if os.path.exists(src):
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copy2(src, dst)
            copied += 1
        else:
            missing.append(rel)
    return copied, missing


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("notebooks", nargs="*", default=NOTEBOOKS)
    p.add_argument("--out", required=True, help="report folder (created)")
    p.add_argument("--env", default=os.path.join(PROJECT, "pbody_RNA_Spang.scicore.env"))
    a = p.parse_args()
    load_dotenv(a.env, override=True)
    figures_dir = os.getenv("FIGURES_DIR")
    os.makedirs(a.out, exist_ok=True)

    rows = []
    for nb in a.notebooks:
        nb_path = nb if os.path.isabs(nb) else os.path.join(PROJECT, nb)
        name, title, body, n_exec, n_code = export(nb_path, a.out)
        copied, missing = copy_figures(body, figures_dir, a.out)
        print(f"{name}: {n_exec}/{n_code} code cells executed, {copied} figure files copied"
              + (f", {len(missing)} linked files missing (e.g. {missing[0]})" if missing else ""))
        rows.append(f'<li><a href="{name}">{html.escape(title)}</a> <small>({n_exec}/{n_code} code cells executed)</small></li>')

    with open(os.path.join(a.out, "index.html"), "w") as fh:
        fh.write(f"""<!doctype html><html><head><meta charset="utf-8"><title>P-body RNA analysis reports</title>
<style>body{{font-family:sans-serif;max-width:50em;margin:2em auto;padding:0 1em;line-height:1.5}}</style></head><body>
<h1>P-body RNA analysis reports</h1>
<p>Generated {datetime.date.today().isoformat()}. Each report shows the analysis text, tables and figures without code.
Figure names link to the full-size files (PNG; PDF where available) in <code>figures/</code>.</p>
<ul>{''.join(rows)}</ul></body></html>""")
    print("index:", os.path.join(a.out, "index.html"))


if __name__ == "__main__":
    main()
