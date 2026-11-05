# Orthorhombify documentation

Guides live here as **Markdown** (easy to read on GitHub) and are built into an HTML
site with **Sphinx** + **MyST**, published at https://awwal41.github.io/orthorhombify/.

| Guide | When to read it |
|-------|-----------------|
| [Getting started](getting-started.md) | **Start here**: install, convert one crystal, make LAMMPS files |
| [Installation](install.md) | pip, venv, conda, HPC, from source |
| [Orthorhombify](orthorhombify.md) | How the supercell search works and how to tune it |
| [LAMMPS files and charges](lammps.md) | Starting charges, neutralizing, what the input script does |
| [Using it from Python](python-api.md) | The Python API with worked examples |
| [Examples](examples.md) | beta-HMX, MDNT, TNT and FOX-7 |
| [CLI reference](cli.md) | Every option of both commands |
| [`api.rst`](api.rst) | Function reference (built from the docstrings) |

## Build the HTML locally

```bash
pip install -r docs/requirements.txt
sphinx-build -b html docs docs/_build/html
```

Open `docs/_build/html/index.html`.

## Publish

The site is served by GitHub Pages from the `gh-pages` branch. After editing the docs,
build and push the HTML (PowerShell or bash, from the repo root):

```bash
sphinx-build -b html docs docs/_build/html
cd docs/_build/html
git init -b gh-pages
touch .nojekyll        # PowerShell: New-Item .nojekyll
git add -A
git commit -m "Update docs"
git push -f https://github.com/Awwal41/orthorhombify.git gh-pages
```

A `.readthedocs.yaml` is also included if you prefer Read the Docs.
