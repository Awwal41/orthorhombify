# Installation

Orthorhombify is a normal Python package. It needs **Python 3.9 or newer** and pulls
in **NumPy** and **ASE** by itself.

---

## Pick one way

| You have | Run |
|----------|-----|
| `git` installed (most people) | `pip install git+https://github.com/Awwal41/orthorhombify.git` |
| No `git` | `pip install https://github.com/Awwal41/orthorhombify/archive/refs/heads/main.zip` |
| A specific release | `pip install https://github.com/Awwal41/orthorhombify/releases/download/v2.1.0/orthorhombify-2.1.0-py3-none-any.whl` |
| You want to edit the code | see [From source](#from-source) below |

Upgrade later with the same command plus `--upgrade`.

---

## A private Python folder (recommended)

A **venv** keeps these libraries away from the rest of your computer.

**Windows (PowerShell):**

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install git+https://github.com/Awwal41/orthorhombify.git
```

If Windows blocks the script, run this once, then try Activate again:

```powershell
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
```

**macOS / Linux:**

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install git+https://github.com/Awwal41/orthorhombify.git
```

When it worked, the prompt starts with `(.venv)`.

**Conda / Anaconda:** activate your environment, then the same `pip install` line works.

**HPC cluster:** load a Python module first (for example `module load python`), then use
a venv as above, or `pip install --user git+https://github.com/Awwal41/orthorhombify.git`.

---

## Check

```bash
orthorhombify --version
orthorhombify-lammps --version
python -c "import orthorhombify; print(orthorhombify.__version__)"
```

All three should print `2.1.0`.

---

## From source

To change the code, or run the tests:

```bash
git clone https://github.com/Awwal41/orthorhombify.git
cd orthorhombify
pip install -e ".[test]"
pytest
```

`-e` (editable) means your edits take effect without reinstalling.

You can also run it straight from the clone, without installing, as long as NumPy and
ASE are installed:

```bash
python orthorhombify_v2.py my_crystal.cif
python make_lammps_inputs.py output/my_crystal_orth_1.cif
```

---

## If something goes wrong

| What you see | What to do |
|--------------|------------|
| `pip` is not recognized | Use `python -m pip install ...` |
| `git` is not recognized | Use the `.zip` line from the table above |
| `orthorhombify` is not recognized | Open a **new** terminal, or use `python -m orthorhombify` |
| `orthorhombify-lammps` is not recognized | `python -c "from orthorhombify.cli import lammps_main; lammps_main()" FILE` |
| `Requires-Python >=3.9` | Install a newer Python |
| Permission denied (cluster / shared Python) | Add `--user`, or use a venv |

## Uninstall

```bash
pip uninstall orthorhombify
```
