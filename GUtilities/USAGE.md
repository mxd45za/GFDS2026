# Python Plotting

`slice_contour.py` reads FDS slice output through the `.smv` file and generates a PNG together with the plotted data in CSV format.

## Installation

From the GFDS repository root:

```bash
python3 -m venv .gfdstools
source .gfdstools/bin/activate
python -m pip install --upgrade pip
python -m pip install -r GUtilities/requirements.txt
export GFDS_ROOT="$(pwd)"
```

## Examples

Run the utility from the directory containing the selected `.smv` file.

Instantaneous temperature contour:

```bash
python "$GFDS_ROOT/GUtilities/slice_contour.py" case.smv \
  --qua temperature --plane y --time 30 --smooth
```

Time-averaged velocity contour:

```bash
python "$GFDS_ROOT/GUtilities/slice_contour.py" case.smv \
  --qua velocity --time 20 30 --range 0 5 --smooth
```

Time-averaged velocity with equal-length direction arrows:

```bash
python "$GFDS_ROOT/GUtilities/slice_contour.py" case.smv \
  --qua velocity --time 200 220 --smooth --vector \
  --vector-normalize --skip 5 --scale 200 --arrowsize 1 \
  --vector-linewidth 3
```

One value after `--time` selects the nearest saved frame. Two values calculate a time average. For a requested plane coordinate, the utility uses the nearest available slice.

For quantities, planes, coordinate normalization, color ranges, vector controls, obstacle masks, output names, and all other options, run:

```bash
python "$GFDS_ROOT/GUtilities/slice_contour.py" --help
```

## Python Environment

For `ModuleNotFoundError: No module named 'fdsreader'`, activate the environment and install the requirements:

```bash
source "$GFDS_ROOT/.gfdstools/bin/activate"
python -m pip install -r "$GFDS_ROOT/GUtilities/requirements.txt"
```

For `error: externally-managed-environment`, install the packages in the `.gfdstools` environment created above.
