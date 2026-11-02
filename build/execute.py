import time
from pathlib import Path
import nbformat
from nbclient import NotebookClient

path = Path("deliverables/spike_variation_analysis.ipynb")
nb = nbformat.read(str(path), as_version=4)
client = NotebookClient(nb, timeout=5400, kernel_name="python3",
                        resources={"metadata": {"path": "."}},
                        allow_errors=False)
t0 = time.time()
print(f"executing {len(nb.cells)} cells", flush=True)
try:
    client.execute()
    print(f"OK in {time.time()-t0:.0f}s", flush=True)
except Exception as exc:
    print(f"FAILED after {time.time()-t0:.0f}s: {type(exc).__name__}", flush=True)
    print(str(exc)[:3000], flush=True)
finally:
    nbformat.write(nb, str(path))
    print("notebook written with outputs", flush=True)
