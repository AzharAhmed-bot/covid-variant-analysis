import json
from pathlib import Path
import nbformat

cells = []
for part in ["part1", "part2", "part3", "part4"]:
    cells.extend(json.loads(Path(f"build/cells_{part}.json").read_text()))

nb = nbformat.v4.new_notebook()
nb.cells = [
    nbformat.v4.new_markdown_cell(c["source"]) if c["cell_type"] == "markdown"
    else nbformat.v4.new_code_cell(c["source"])
    for c in cells
]
nb.metadata = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python"},
    "colab": {"provenance": [], "toc_visible": True},
}
out = Path("deliverables/spike_variation_analysis.ipynb")
out.parent.mkdir(exist_ok=True)
nbformat.write(nb, str(out))
n_md = sum(1 for c in nb.cells if c.cell_type == "markdown")
n_code = len(nb.cells) - n_md
print(f"wrote {out}: {len(nb.cells)} cells ({n_md} markdown, {n_code} code)")
