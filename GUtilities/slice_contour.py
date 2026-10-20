#!/usr/bin/env python3
"""Plot instantaneous or time-averaged FDS SLCF output."""

from __future__ import annotations

import argparse
import csv
import math
import re
import struct
import textwrap
import warnings
from pathlib import Path

import fdsreader
from fdsreader import settings as fdsreader_settings
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator, LinearLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable


# fdsreader eagerly loads optional metadata that this slice plotter does not
# use. Keep malformed optional modules and header-only step logs silent; a
# missing or invalid requested SLCF still fails in the plotting path below.
fdsreader_settings.IGNORE_ERRORS = True
warnings.filterwarnings(
    "ignore",
    message=r'loadtxt: input contained no data: ".*_steps\.csv"',
    category=UserWarning,
)


def resolve_smv_path(smv_arg: str | Path) -> Path:
    path = Path(smv_arg)
    candidates = []
    if path.suffix.lower() == ".smv":
        candidates.append(path)
    else:
        candidates.append(path.with_suffix(".smv"))
        candidates.append(path)
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return path


def parse_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in ("1", "true", "t", "yes", "y", "on"):
        return True
    if text in ("0", "false", "f", "no", "n", "off"):
        return False
    raise argparse.ArgumentTypeError(f"invalid boolean value: {value}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot an FDS SLCF slice. Use --time T for one instantaneous frame, "
            "or --time T0 T1 for a time average."
        ),
        epilog=textwrap.dedent(
            """\
            Common --qua names:
              TEMPERATURE
              VELOCITY
              DENSITY
              OXYGEN MASS FRACTION
              SOOT MASS FRACTION
              SOOT DENSITY
              HRRPUV
              INTEGRATED INTENSITY

            Examples:
              python3 slice_contour.py case.smv --qua TEMPERATURE --plane y --time 30
              python3 slice_contour.py case.smv --qua TEMPERATURE --plane y --time 30 40
              python3 slice_contour.py case.smv --qua TEMPERATURE --plane y --time 30 --range 0 45
              python3 slice_contour.py case.smv --qua VELOCITY --plane z 58 --time 30 --range 0 1
              python3 slice_contour.py case.smv --qua RADI --plane y --time 30
              python3 slice_contour.py case.smv --qua "OXYGEN MASS FRACTION" --plane y 1.2 --time 30
              python3 slice_contour.py case.smv --qua "SOOT DENSITY" --plane y --time 30 40
            """
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "smv",
        help="Path to the FDS .smv file.",
    )

    quantity_group = parser.add_argument_group("quantity")
    quantity_group.add_argument(
        "--qua",
        dest="quantity",
        default="TEMPERATURE",
        help="FDS output quantity. Default: TEMPERATURE.",
    )
    quantity_group.add_argument(
        "--quantity",
        dest="quantity",
        help=argparse.SUPPRESS,
    )

    parser.add_argument(
        "--plane",
        nargs="+",
        default=["y", "0"],
        metavar=("DIR", "VALUE"),
        help=(
            "Slice-normal coordinate and optional coordinate. Examples: "
            "--plane y, --plane y 1.2, --plane x 0.5. "
            "Default is --plane y 0. If VALUE is omitted, the nearest slice to 0 is used."
        ),
    )
    parser.add_argument(
        "--plane-value",
        type=float,
        default=None,
        dest="legacy_plane_value",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--t0",
        type=float,
        default=None,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--t1",
        type=float,
        default=None,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--time",
        type=float,
        nargs="+",
        default=None,
        metavar=("T0", "T1"),
        help=(
            "One value means instantaneous frame nearest to T0. "
            "Two values mean time average from T0 to T1."
        ),
    )
    parser.add_argument(
        "--output-prefix",
        default=None,
        help="Output prefix. Default includes the selected time interval.",
    )
    parser.add_argument(
        "--normalize",
        action="store_true",
        help=(
            "Normalize the plotted field. Signed fields are scaled to [-1, 1]; "
            "non-negative magnitude fields are scaled to [0, 1]."
        ),
    )
    parser.add_argument(
        "--g1",
        dest="normalize",
        action="store_true",
        help="Compatibility alias for --normalize.",
    )
    parser.add_argument(
        "--g",
        "--axisg1",
        "--normalize-coordinates",
        dest="normalize_coordinates",
        action="store_true",
        help="Normalize plotted coordinates by the full plotted axis lengths. Alias: --axisg1.",
    )
    parser.add_argument(
        "--axis",
        dest="show_axis",
        action="store_true",
        default=True,
        help="Show x/y axis titles and ticks. Default.",
    )
    parser.add_argument(
        "--noaxis",
        dest="show_axis",
        action="store_false",
        help="Hide x/y axis titles while retaining coordinate tick marks.",
    )
    parser.add_argument(
        "--cmap",
        default="jet",
        help="Matplotlib colormap. Default: jet.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=1200,
        help="PNG output resolution. Default: 1200.",
    )
    parser.add_argument(
        "--range",
        nargs=2,
        type=float,
        metavar=("MIN", "MAX"),
        dest="plot_range",
        default=None,
        help="Fixed color range for the plot. Example: --range 0 45.",
    )
    parser.add_argument(
        "--smooth",
        action="store_true",
        help="Smooth the PNG display with bilinear interpolation. CSV data are unchanged.",
    )
    parser.add_argument(
        "--vector",
        action="store_true",
        help="Overlay velocity vectors on the slice plot.",
    )
    parser.add_argument(
        "--vector-normalize",
        "--direction-only",
        action="store_true",
        help=(
            "Draw equal-length arrows that show direction only; the color field still shows "
            "velocity magnitude. Recommended for long tunnel plots."
        ),
    )
    parser.add_argument(
        "--skip",
        type=int,
        default=1,
        help="Vector sampling stride used with --vector. Default: 1.",
    )
    parser.add_argument(
        "--scale",
        type=float,
        default=80.0,
        help="Vector scale used with --vector. Larger values make arrows shorter. Default: 80.",
    )
    parser.add_argument(
        "--arrowsize",
        type=float,
        default=1.0,
        help="Vector shaft thickness multiplier used with --vector. Default: 1.",
    )
    parser.add_argument(
        "--vector-color",
        default="white",
        help="Vector color, or auto for local black/white contrast. Default: white.",
    )
    parser.add_argument(
        "--vector-edgecolor",
        default="none",
        help="Optional vector outline color. Default: none.",
    )
    parser.add_argument(
        "--vector-linewidth",
        type=float,
        default=0.0,
        help="Vector outline width in points. Default: 0.",
    )
    orientation_group = parser.add_mutually_exclusive_group()
    orientation_group.add_argument(
        "--horizon",
        "-horizon",
        dest="plot_orientation",
        action="store_const",
        const="horizon",
        default="horizon",
        help="Use a horizontal/landscape figure layout. Default.",
    )
    orientation_group.add_argument(
        "--vertical",
        "-vertical",
        dest="plot_orientation",
        action="store_const",
        const="vertical",
        help="Use a vertical/portrait figure layout.",
    )
    parser.add_argument(
        "--no-soot-clean",
        action="store_true",
        help="Do not mask anomalously large SOOT DENSITY values.",
    )
    parser.add_argument(
        "--obst",
        nargs="?",
        const=True,
        default=True,
        type=parse_bool,
        help=(
            "Mask OBST cells in gray. Default: 1. Use --obst 0 to disable. "
            "Geometry from the input file is used first, with velocity fallback."
        ),
    )
    parser.add_argument(
        "--obst-eps",
        type=float,
        default=1.0e-6,
        help="Velocity threshold used by --obst. Default: 1e-6.",
    )
    parser.add_argument(
        "--obst-quantity",
        default="VELOCITY",
        help="Quantity used to build the OBST mask for non-velocity plots. Default: VELOCITY.",
    )
    parser.add_argument(
        "--obst-source",
        choices=("auto", "geometry", "velocity"),
        default="auto",
        help="OBST mask source. Default auto uses input-file geometry, then velocity fallback.",
    )

    args = parser.parse_args()
    # --g1/--normalize affects only the field/color values. Coordinate
    # normalization is deliberately independent and uses --axisg1.
    quantity_aliases = {
        "tempereature": "TEMPERATURE",
        "temperature": "TEMPERATURE",
        "velocity": "VELOCITY",
        "radi": "INTEGRATED INTENSITY",
        "integrated intensity": "INTEGRATED INTENSITY",
        "integrated_intensity": "INTEGRATED INTENSITY",
    }
    q_key = str(args.quantity).strip().lower()
    if q_key in quantity_aliases:
        args.quantity = quantity_aliases[q_key]
    if args.plot_range is not None and args.plot_range[1] <= args.plot_range[0]:
        parser.error("--range expects MIN MAX with MAX > MIN")
    if args.dpi < 100:
        parser.error("--dpi must be >= 100")
    if args.skip < 1:
        parser.error("--skip must be >= 1")
    if args.scale <= 0.0:
        parser.error("--scale must be > 0")
    if args.arrowsize <= 0.0:
        parser.error("--arrowsize must be > 0")
    if args.vector_linewidth < 0.0:
        parser.error("--vector-linewidth must be >= 0")
    if len(args.plane) not in (1, 2):
        parser.error("--plane expects DIR or DIR VALUE, for example --plane y or --plane y 1.2")
    args.plane_dir = args.plane[0].lower()
    if args.plane_dir not in ("x", "y", "z"):
        parser.error("--plane DIR must be x, y, or z")
    if len(args.plane) == 2:
        try:
            args.plane_value = float(args.plane[1])
        except ValueError:
            parser.error("--plane VALUE must be a number, for example --plane y 1.2")
    else:
        args.plane_value = 0.0 if args.legacy_plane_value is None else args.legacy_plane_value

    if args.time is not None and len(args.time) not in (1, 2):
        parser.error("--time expects one value for instantaneous output or two values for average")
    return args


def slice_plane_value(slc, plane: str) -> float:
    extent = slc.extent
    if plane == "x":
        return float(extent.x_start)
    if plane == "y":
        return float(extent.y_start)
    return float(extent.z_start)


def slice_matches_plane(slc, plane: str) -> bool:
    return {"x": 1, "y": 2, "z": 3}[plane] == slc.orientation


def select_slice(
    sim: fdsreader.Simulation,
    quantity: str,
    plane: str,
    value: float | None,
):
    slices = sim.slices.filter_by_quantity(quantity)
    if len(slices) == 0:
        available = sorted({slc.quantity.name for slc in sim.slices})
        raise RuntimeError(
            f"Cannot find SLCF quantity {quantity!r}. Available quantities: {available}"
        )

    plane_slices = [slc for slc in slices if slice_matches_plane(slc, plane)]
    if len(plane_slices) == 0:
        available = [
            f"{slc.quantity.name} at "
            f"x={slice_plane_value(slc, 'x'):.8g}, "
            f"y={slice_plane_value(slc, 'y'):.8g}, "
            f"z={slice_plane_value(slc, 'z'):.8g}"
            for slc in slices
        ]
        raise RuntimeError(
            f"Cannot find a {quantity!r} slice normal to {plane!r}. "
            f"Available {quantity!r} slices: {available}"
        )

    if value is None:
        slc = min(plane_slices, key=lambda item: slice_plane_value(item, plane))
        return slc, slice_plane_value(slc, plane)

    slc = min(plane_slices, key=lambda item: abs(slice_plane_value(item, plane) - value))
    return slc, slice_plane_value(slc, plane)


def load_slice_values(
    slc,
    mode: str,
    args: argparse.Namespace,
    requested_end: float | None,
):
    # fdsreader stores one SLCF as several subslices when the requested plane
    # crosses multiple FDS meshes.  Convert those pieces to one global array
    # before doing time averaging or plotting.  The old single-mesh path is
    # retained because it supports the direct .sf reader used for incomplete
    # cached time ranges.
    multi_mesh = len(slc.subslices) > 1
    if multi_mesh:
        global_result = slc.to_global(return_coordinates=False)
        if isinstance(global_result, tuple) and isinstance(global_result[0], tuple):
            # Cell-centred slices can legitimately produce two global grids
            # at a mesh boundary.  Keep the first physical side, matching the
            # first subslice used by fdsreader for its time axis.
            data = np.asarray(global_result[0][0], dtype=float)
        elif isinstance(global_result, tuple):
            data = np.asarray(global_result[0], dtype=float)
        else:
            data = np.asarray(global_result, dtype=float)
    else:
        data = np.asarray(slc.subslices[0].data, dtype=float)
    times = np.asarray(slc.times, dtype=float)
    direct_info = None
    if (not multi_mesh) and requested_end is not None and times.size and requested_end > times[-1] + 1.0e-6:
        direct_data, direct_times = read_sf_direct(slc.subslices[0].filename, data.shape[1:])
        if direct_times[-1] > times[-1] + 1.0e-6:
            data = direct_data
            times = direct_times
            direct_info = f"direct .sf reader used, full time range [{times[0]:.8g}, {times[-1]:.8g}]"
    elif multi_mesh:
        direct_info = f"merged {len(slc.subslices)} mesh subslices"

    if mode == "legacy_timeavg":
        values, selected_times = average_slice(data, times, args.t0, args.t1)
        return values, selected_times, None, times, direct_info
    if mode == "instant":
        values, selected_times, frame_index = instantaneous_slice(data, times, args.time[0])
        return values, selected_times, frame_index, times, direct_info
    values, selected_times = average_slice(data, times, args.time[0], args.time[1])
    return values, selected_times, None, times, direct_info


def slice_coordinates(slc):
    """Return plotting coordinates for either one mesh or a merged slice."""

    if len(slc.subslices) > 1:
        global_result = slc.to_global(return_coordinates=True)
        if isinstance(global_result[0], tuple):
            return global_result[1][0]
        return global_result[1]
    return slc.subslices[0].get_coordinates()


def vector_component_names(plane: str):
    if plane == "x":
        return "V-VELOCITY", "W-VELOCITY"
    if plane == "y":
        return "U-VELOCITY", "W-VELOCITY"
    return "U-VELOCITY", "V-VELOCITY"


def parse_smv_slcf_records(smv_path: str | Path):
    records = []
    lines = Path(smv_path).read_text(errors="ignore").splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line.startswith("SLCF"):
            i += 1
            continue
        before_comment = line.split("!", 1)[0]
        nums = re.findall(r"[-+]?\d+", before_comment)
        if len(nums) < 7:
            i += 1
            continue
        try:
            bounds = tuple(int(v) for v in nums[-6:])
        except ValueError:
            i += 1
            continue
        orient = None
        if "!" in line:
            comment_nums = re.findall(r"[-+]?\d+", line.split("!", 1)[1])
            if comment_nums:
                orient = int(comment_nums[-1])
        if i + 4 >= len(lines):
            i += 1
            continue
        records.append(
            {
                "bounds": bounds,
                "orientation": orient,
                "filename": lines[i + 1].strip(),
                "quantity": lines[i + 2].strip().upper(),
            }
        )
        i += 5
    return records


def load_smv_vector_component(
    smv_path: str | Path,
    base_filename: str | Path,
    quantity: str,
    mode: str,
    args: argparse.Namespace,
    requested_end: float | None,
    value_shape: tuple[int, ...],
):
    smv_path = Path(smv_path)
    records = parse_smv_slcf_records(smv_path)
    base_name = Path(base_filename).name
    base_record = next((rec for rec in records if Path(rec["filename"]).name == base_name), None)
    if base_record is None:
        raise RuntimeError(f"Cannot find base SLCF record for {base_name!r} in {smv_path.name}.")
    match = next(
        (
            rec
            for rec in records
            if rec["quantity"] == quantity.upper()
            and rec["bounds"] == base_record["bounds"]
            and rec["orientation"] == base_record["orientation"]
        ),
        None,
    )
    if match is None:
        raise RuntimeError(
            f"Cannot find {quantity} vector component matching {base_name!r} in {smv_path.name}."
        )
    filename = smv_path.parent / match["filename"]
    data, times = read_sf_direct(filename, value_shape)
    if mode == "legacy_timeavg":
        values, selected_times = average_slice(data, times, args.t0, args.t1)
    elif mode == "instant":
        values, selected_times, _ = instantaneous_slice(data, times, args.time[0])
    else:
        values, selected_times = average_slice(data, times, args.time[0], args.time[1])
    info = (
        f"{quantity} from {filename.name}, time range "
        f"[{selected_times[0]:.8g}, {selected_times[-1]:.8g}]"
    )
    return values, info


def load_smv_vector_component_multimesh(
    smv_path: str | Path,
    base_slice,
    quantity: str,
    mode: str,
    args: argparse.Namespace,
    requested_end: float | None,
):
    """Read and merge a vector component for every mesh subslice.

    Load VECTOR=.TRUE. companion files into their scalar SubSlices, merge
    them with ``to_global``, then restore the scalar caches.
    """

    key = {"U-VELOCITY": "u", "V-VELOCITY": "v", "W-VELOCITY": "w"}[quantity.upper()]
    saved_data = []
    try:
        for subslice in base_slice.subslices:
            filename = subslice.vector_filenames.get(key)
            if filename is None:
                raise RuntimeError(
                    f"Cannot find {quantity} companion file for {Path(subslice.filename).name}."
                )
            local_shape = np.asarray(subslice.data).shape[1:]
            component_data, component_times = read_sf_direct(
                Path(smv_path).parent / filename, local_shape
            )
            saved_data.append((subslice, subslice._data))
            subslice._data = component_data

        global_result = base_slice.to_global(return_coordinates=False)
        if isinstance(global_result, tuple) and isinstance(global_result[0], tuple):
            data = np.asarray(global_result[0][0], dtype=float)
        elif isinstance(global_result, tuple):
            data = np.asarray(global_result[0], dtype=float)
        else:
            data = np.asarray(global_result, dtype=float)
        times = np.asarray(base_slice.times, dtype=float)

    finally:
        for subslice, old_data in saved_data:
            subslice._data = old_data

    if mode == "legacy_timeavg":
        values, selected_times = average_slice(data, times, args.t0, args.t1)
    elif mode == "instant":
        values, selected_times, _ = instantaneous_slice(data, times, args.time[0])
    else:
        values, selected_times = average_slice(data, times, args.time[0], args.time[1])

    info = f"{quantity} merged {len(base_slice.subslices)} mesh subslices"
    return values, info


def load_vector_overlay(
    sim: fdsreader.Simulation,
    args: argparse.Namespace,
    plane_value: float,
    mode: str,
    requested_end: float | None,
    value_shape: tuple[int, ...],
    base_filename: str | Path | None = None,
    base_slice=None,
):
    q1, q2 = vector_component_names(args.plane_dir)
    try:
        slc1, plane_value_1 = select_slice(sim, q1, args.plane_dir, args.plane_value)
        slc2, plane_value_2 = select_slice(sim, q2, args.plane_dir, args.plane_value)
    except RuntimeError:
        if base_slice is not None and len(base_slice.subslices) > 1:
            v1, info1 = load_smv_vector_component_multimesh(
                args.smv, base_slice, q1, mode, args, requested_end
            )
            v2, info2 = load_smv_vector_component_multimesh(
                args.smv, base_slice, q2, mode, args, requested_end
            )
            return v1, v2, [info1, info2]
        if base_filename is None:
            raise
        v1, info1 = load_smv_vector_component(
            args.smv, base_filename, q1, mode, args, requested_end, value_shape
        )
        v2, info2 = load_smv_vector_component(
            args.smv, base_filename, q2, mode, args, requested_end, value_shape
        )
        return v1, v2, [info1, info2]

    for quantity, component_plane in ((q1, plane_value_1), (q2, plane_value_2)):
        if not math.isclose(component_plane, plane_value, rel_tol=1.0e-9, abs_tol=1.0e-9):
            raise RuntimeError(
                f"Vector plane mismatch: {quantity} at {component_plane:g} m, "
                f"current slice at {plane_value:g} m"
            )

    v1, _, _, _, info1 = load_slice_values(slc1, mode, args, requested_end)
    v2, _, _, _, info2 = load_slice_values(slc2, mode, args, requested_end)
    if v1.shape != value_shape or v2.shape != value_shape:
        raise RuntimeError(
            f"Vector shape mismatch: {q1} {v1.shape}, {q2} {v2.shape}, field {value_shape}"
        )
    infos = [info for info in (info1, info2) if info is not None]
    return v1, v2, infos


def average_slice(data: np.ndarray, times: np.ndarray, t0: float | None, t1: float | None):
    if times.size == 0:
        raise RuntimeError("No SLCF frames were found.")

    start = times[0] if t0 is None else t0
    end = times[-1] if t1 is None else t1
    if end < start:
        raise RuntimeError(f"Invalid averaging interval: t1={end} < t0={start}.")

    mask = (times >= start) & (times <= end)
    if not np.any(mask):
        raise RuntimeError(
            f"No SLCF frames inside interval [{start}, {end}]. "
            f"Available time range is [{times[0]}, {times[-1]}]."
        )

    selected_times = times[mask]
    selected_data = data[mask]
    if selected_times.size == 1:
        return selected_data[0], selected_times

    # FDS may use adaptive time stepping; use time-weighted average.
    duration = selected_times[-1] - selected_times[0]
    if duration <= 0:
        return selected_data.mean(axis=0), selected_times
    return np.trapezoid(selected_data, selected_times, axis=0) / duration, selected_times


def instantaneous_slice(data: np.ndarray, times: np.ndarray, time: float):
    if times.size == 0:
        raise RuntimeError("No SLCF frames were found.")
    index = int(np.argmin(np.abs(times - time)))
    return data[index], np.asarray([times[index]], dtype=float), index


def read_fortran_record(fh):
    raw = fh.read(4)
    if not raw:
        return None
    n = struct.unpack("<i", raw)[0]
    data = fh.read(n)
    end_raw = fh.read(4)
    if len(data) != n or len(end_raw) != 4:
        raise RuntimeError("Unexpected end of slice file.")
    end = struct.unpack("<i", end_raw)[0]
    if end != n:
        raise RuntimeError("Malformed Fortran record in slice file.")
    return data


def reshape_sf_record(data_record: bytes, shape: tuple[int, ...], filename: str | Path):
    values = np.frombuffer(data_record, dtype="<f4")
    expected = int(np.prod(shape))
    if values.size == expected:
        return values.reshape(shape, order="F")

    if len(shape) == 2:
        n1, n2 = shape
        padded_shapes = (
            (n1 + 1, n2 + 1),
            (n1 + 1, n2),
            (n1, n2 + 1),
        )
        for padded_shape in padded_shapes:
            if values.size == int(np.prod(padded_shape)):
                padded = values.reshape(padded_shape, order="F")
                i0 = padded_shape[0] - n1
                j0 = padded_shape[1] - n2
                return padded[i0:i0 + n1, j0:j0 + n2]

    raise RuntimeError(
        f"Cannot reshape {Path(filename).name}: record has {values.size} float32 values, "
        f"expected shape {shape} ({expected} values)."
    )


def read_sf_direct(filename: str | Path, shape: tuple[int, ...]):
    times = []
    frames = []
    with Path(filename).open("rb") as fh:
        for _ in range(4):
            if read_fortran_record(fh) is None:
                raise RuntimeError("Slice file ended before header was complete.")
        while True:
            time_record = read_fortran_record(fh)
            if time_record is None:
                break
            if len(time_record) != 4:
                raise RuntimeError("Unexpected time record size in slice file.")
            data_record = read_fortran_record(fh)
            if data_record is None:
                break
            times.append(struct.unpack("<f", time_record)[0])
            frames.append(reshape_sf_record(data_record, shape, filename))
    if not frames:
        raise RuntimeError("No frames were read from slice file.")
    return np.asarray(frames, dtype=float), np.asarray(times, dtype=float)


def fmt_time(value: float) -> str:
    return f"{value:g}"


def clean_soot_density(values: np.ndarray):
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return values, None
    cleaned = values.copy()
    high_density_count = int(np.count_nonzero(cleaned > 1.0))
    if high_density_count > 0:
        cleaned[cleaned > 1.0] = 0.0

    p99 = float(np.nanpercentile(finite, 99.0))
    vmax = float(np.nanmax(finite))
    if p99 <= 0.0 or vmax <= 1000.0 * p99:
        if high_density_count > 0:
            return cleaned, {
                "cutoff": 1.0,
                "reset": high_density_count,
                "vmax": vmax,
                "mode": "absolute",
            }
        return values, None

    # SOOT DENSITY slice files can include ghost/boundary artifacts that dominate
    # the color scale. Reset only the extreme tail, keeping the physical plume field.
    cutoff = p99
    count = int(np.count_nonzero(cleaned > cutoff))
    cleaned[cleaned > cutoff] = 0.0
    return cleaned, {
        "cutoff": cutoff,
        "reset": high_density_count + count,
        "vmax": vmax,
        "mode": "percentile",
    }


def input_file_from_smv(smv_path: Path) -> Path | None:
    smv_path = resolve_smv_path(smv_path)
    try:
        lines = smv_path.read_text(errors="ignore").splitlines()
    except OSError:
        return None
    for i, line in enumerate(lines):
        if line.strip().upper() == "INPF":
            for item in lines[i + 1 :]:
                name = item.strip()
                if name:
                    candidate = Path(name)
                    if not candidate.is_absolute():
                        candidate = smv_path.parent / candidate
                    return candidate
    return None


def iter_fds_namelists(fds_path: Path, name: str):
    target = "&" + name.upper()
    active = False
    record = []
    for raw in fds_path.read_text(errors="ignore").splitlines():
        line = raw.strip()
        if not active and line.upper().startswith(target):
            active = True
            record = [line]
        elif active:
            record.append(line)
        if active and "/" in line:
            yield " ".join(record)
            active = False
            record = []


def xb_from_record(record: str):
    match = re.search(r"\bXB\s*=\s*([^/]+?)(?=,\s*[A-Z_][A-Z0-9_]*\s*=|/)", record, re.IGNORECASE)
    if match is None:
        return None
    values = []
    for token in match.group(1).split(",")[:6]:
        try:
            values.append(float(token))
        except ValueError:
            return None
    return values if len(values) == 6 else None


def apply_xb_mask(mask, xb, fill_value, plane_dir, plane_value, axis1_name, axis1, axis2_name, axis2):
    bounds = {
        "x": (xb[0], xb[1]),
        "y": (xb[2], xb[3]),
        "z": (xb[4], xb[5]),
    }
    eps = 1.0e-9
    lo, hi = bounds[plane_dir]
    if plane_value < lo - eps or plane_value > hi + eps:
        return
    a1_lo, a1_hi = bounds[axis1_name]
    a2_lo, a2_hi = bounds[axis2_name]
    m1 = (axis1 >= a1_lo - eps) & (axis1 <= a1_hi + eps)
    m2 = (axis2 >= a2_lo - eps) & (axis2 <= a2_hi + eps)
    if np.any(m1) and np.any(m2):
        mask[np.ix_(m1, m2)] = fill_value


def geometry_obst_mask(smv_path, plane_dir, plane_value, axis1_name, axis1, axis2_name, axis2):
    fds_path = input_file_from_smv(Path(smv_path).resolve())
    if fds_path is None or not fds_path.exists():
        return None, None
    mask = np.zeros((axis1.size, axis2.size), dtype=bool)
    for record in iter_fds_namelists(fds_path, "OBST"):
        xb = xb_from_record(record)
        if xb is not None:
            apply_xb_mask(mask, xb, True, plane_dir, plane_value, axis1_name, axis1, axis2_name, axis2)
    for record in iter_fds_namelists(fds_path, "HOLE"):
        xb = xb_from_record(record)
        if xb is not None:
            apply_xb_mask(mask, xb, False, plane_dir, plane_value, axis1_name, axis1, axis2_name, axis2)
    return mask, fds_path


def write_csv(
    path: Path,
    axis1_name: str,
    axis1: np.ndarray,
    axis2_name: str,
    axis2: np.ndarray,
    values: np.ndarray,
    extra_fields: list[tuple[str, np.ndarray]] | None = None,
):
    extra_fields = [] if extra_fields is None else extra_fields
    value_mask = np.ma.getmaskarray(values) if np.ma.isMaskedArray(values) else np.zeros(values.shape, dtype=bool)
    value_data = np.ma.filled(values, np.nan) if np.ma.isMaskedArray(values) else values
    prepared_fields = []
    for name, field in extra_fields:
        if field.shape != values.shape:
            raise RuntimeError(f"CSV field {name!r} shape {field.shape} does not match {values.shape}")
        field_mask = value_mask.copy()
        if np.ma.isMaskedArray(field):
            field_mask |= np.ma.getmaskarray(field)
            field_data = np.ma.filled(field, np.nan)
        else:
            field_data = np.array(field, dtype=float, copy=True)
        field_data[field_mask] = np.nan
        prepared_fields.append((name, field_data))
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([axis1_name, axis2_name, "value"] + [name for name, _ in prepared_fields])
        for i1, v1 in enumerate(axis1):
            for i2, v2 in enumerate(axis2):
                row = [f"{v1:.8g}", f"{v2:.8g}", f"{value_data[i1, i2]:.8g}"]
                row.extend(f"{field[i1, i2]:.8g}" for _, field in prepared_fields)
                writer.writerow(row)


def normalize_field(values: np.ndarray):
    """Normalize one plotted field while preserving its mask."""
    masked = np.ma.isMaskedArray(values)
    if masked:
        mask = np.ma.getmaskarray(values)
        data = np.asarray(np.ma.getdata(values), dtype=float)
        finite = data[~mask & np.isfinite(data)]
    else:
        data = np.asarray(values, dtype=float)
        mask = np.zeros(data.shape, dtype=bool)
        finite = data[np.isfinite(data)]

    if finite.size == 0:
        raise RuntimeError("Cannot normalize a field with no finite values.")

    vmin = float(np.min(finite))
    vmax = float(np.max(finite))
    signed = vmin < 0.0 < vmax

    if signed:
        # Keep zero at zero, while stretching each sign branch independently
        # so that a field with a positive offset is still visually resolved.
        normalized = np.zeros_like(data)
        negative = data < 0.0
        positive = data > 0.0
        if vmin < 0.0:
            normalized[negative] = data[negative] / abs(vmin)
        if vmax > 0.0:
            normalized[positive] = data[positive] / vmax
        scale = float(max(abs(vmin), abs(vmax)))
        mode = "signed [-1,1]"
    else:
        # Magnitude fields use field-wise min/max scaling to [0,1].  This is
        # important for temperature-like fields whose absolute values are large
        # but whose spatial variation is comparatively small.
        magnitude = np.abs(data) if vmax <= 0.0 else data
        mag_finite = magnitude[~mask & np.isfinite(magnitude)]
        mag_min = float(np.min(mag_finite))
        mag_max = float(np.max(mag_finite))
        scale = mag_max - mag_min
        normalized = (magnitude - mag_min) / scale if scale > 0.0 else np.zeros_like(data)
        mode = "magnitude [0,1]"

    if masked:
        normalized = np.ma.array(normalized, mask=mask, copy=False)
    return normalized, {
        "mode": mode,
        "scale": scale,
        "raw_min": vmin,
        "raw_max": vmax,
    }


def normalize_vector_components(vec1: np.ndarray, vec2: np.ndarray):
    """Normalize vector components with one common scale to preserve direction."""
    finite = np.concatenate((vec1[np.isfinite(vec1)], vec2[np.isfinite(vec2)]))
    if finite.size == 0:
        raise RuntimeError("Cannot normalize a vector field with no finite values.")
    scale = float(np.max(np.abs(finite)))
    if scale <= 0.0:
        return np.zeros_like(vec1), np.zeros_like(vec2), 0.0
    return vec1 / scale, vec2 / scale, scale


def normalize_vector_directions(vec1: np.ndarray, vec2: np.ndarray):
    """Convert each nonzero in-plane vector to unit length."""
    magnitude = np.hypot(vec1, vec2)
    valid = np.isfinite(magnitude) & (magnitude > np.finfo(float).eps)
    direction1 = np.zeros_like(vec1, dtype=float)
    direction2 = np.zeros_like(vec2, dtype=float)
    direction1[valid] = vec1[valid] / magnitude[valid]
    direction2[valid] = vec2[valid] / magnitude[valid]
    return direction1, direction2


def normalize_axis_coordinates(axis: np.ndarray):
    """Map cell-center coordinates to dimensionless coordinates in [0, 1]."""
    axis = np.asarray(axis, dtype=float)
    if axis.size < 2:
        return np.zeros_like(axis)
    spacing = float(np.median(np.diff(axis)))
    origin = float(axis[0] - 0.5 * spacing)
    length = float(axis[-1] - axis[0] + spacing)
    if length <= 0.0:
        return np.zeros_like(axis)
    return (axis - origin) / length


def nice_step(raw_step: float) -> float:
    if raw_step <= 0.0 or not math.isfinite(raw_step):
        return 1.0
    exponent = math.floor(math.log10(raw_step))
    fraction = raw_step / (10.0**exponent)
    for candidate in (1.0, 2.0, 2.5, 5.0, 10.0):
        if fraction <= candidate:
            return candidate * (10.0**exponent)
    return 10.0 * (10.0**exponent)


def nice_color_limits(values: np.ndarray, target_ticks: int = 12):
    if np.ma.isMaskedArray(values):
        finite = values.compressed()
        finite = finite[np.isfinite(finite)]
    else:
        finite = values[np.isfinite(values)]
    if finite.size == 0:
        return None, None, None

    vmin = float(np.nanmin(finite))
    vmax = float(np.nanmax(finite))
    if not math.isfinite(vmin) or not math.isfinite(vmax):
        return None, None, None
    if math.isclose(vmin, vmax):
        step = nice_step(abs(vmax) / 5.0 if vmax != 0.0 else 1.0)
        return vmin - step, vmax + step, np.arange(vmin - step, vmax + 1.5 * step, step)

    step = nice_step((vmax - vmin) / max(target_ticks - 1, 1))
    lower = math.floor(vmin / step) * step
    upper = math.ceil(vmax / step) * step
    if math.isclose(upper, vmax, rel_tol=1.0e-12, abs_tol=1.0e-12):
        upper += step

    # For positive fields like soot density, mass fraction, HRRPUV, keep zero
    # as the lower bound when the data are close to zero.
    if vmin >= 0.0 and lower < step:
        lower = 0.0

    ticks = np.arange(lower, upper + 0.5 * step, step)
    return lower, upper, ticks


def fixed_color_limits(vmin: float, vmax: float, target_ticks: int = 12):
    step = nice_step((vmax - vmin) / max(target_ticks - 1, 1))
    ticks = np.arange(vmin, vmax + 0.5 * step, step)
    ticks = ticks[(ticks >= vmin - 1.0e-12) & (ticks <= vmax + 1.0e-12)]
    return vmin, vmax, ticks


def fmt_tick(value: float) -> str:
    if math.isclose(value, round(value), abs_tol=1.0e-10):
        return str(int(round(value)))
    return f"{value:g}"


def axis_spacing(axis: np.ndarray) -> float | None:
    if axis.size < 2:
        return None
    diffs = np.diff(np.asarray(axis, dtype=float))
    diffs = diffs[np.isfinite(diffs) & (diffs > 0.0)]
    if diffs.size == 0:
        return None
    return float(np.median(diffs))


def write_plot(
    path: Path,
    axis1_name: str,
    axis1: np.ndarray,
    axis2_name: str,
    axis2: np.ndarray,
    values: np.ndarray,
    title: str,
    colorbar_label: str,
    cmap: str,
    orientation: str,
    plot_range: tuple[float, float] | None,
    smooth: bool,
    dpi: int,
    vector_overlay=None,
    normalization_mode: str | None = None,
    show_axis: bool = False,
):
    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad(color="0.55")
    if plot_range is None and normalization_mode is not None:
        if normalization_mode.startswith("signed"):
            vmin, vmax, ticks = -1.0, 1.0, np.array([-1.0, 0.0, 1.0])
        else:
            vmin, vmax, ticks = 0.0, 1.0, np.array([0.0, 0.5, 1.0])
    elif plot_range is None:
        vmin, vmax, ticks = nice_color_limits(values)
    else:
        vmin, vmax, ticks = fixed_color_limits(plot_range[0], plot_range[1])

    figsize = (5.0, 7.2) if orientation == "vertical" else (9.2, 5.6)
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=(orientation != "horizon"))
    if orientation == "horizon":
        fig.subplots_adjust(left=0.10, right=0.88, bottom=0.11, top=0.80)
    if smooth:
        image = ax.imshow(
            values.T,
            origin="lower",
            extent=(float(axis1[0]), float(axis1[-1]), float(axis2[0]), float(axis2[-1])),
            interpolation="bilinear",
            cmap=cmap_obj,
            vmin=vmin,
            vmax=vmax,
            aspect="auto",
        )
    else:
        image = ax.pcolormesh(
            axis1,
            axis2,
            values.T,
            shading="auto",
            cmap=cmap_obj,
            vmin=vmin,
            vmax=vmax,
        )
    if vector_overlay is not None:
        vec1, vec2, skip, scale, arrowsize, color, edgecolor, linewidth, scale_units, mask = vector_overlay
        x_grid, y_grid = np.meshgrid(axis1, axis2, indexing="ij")
        q1 = np.ma.array(vec1, mask=mask)
        q2 = np.ma.array(vec2, mask=mask)
        arrow_color = color
        if str(color).lower() == "auto":
            # Use a clean single-color arrow at each location, choosing black
            # or white from the local background luminance.  This avoids the
            # heavy outlined pointer shape while remaining legible across the
            # full colormap.
            sampled_values = np.ma.asarray(values)[::skip, ::skip]
            background_rgba = image.cmap(image.norm(sampled_values.filled(vmin)))
            luminance = (
                0.2126 * background_rgba[..., 0]
                + 0.7152 * background_rgba[..., 1]
                + 0.0722 * background_rgba[..., 2]
            )
            dark = np.array([0.08, 0.08, 0.08, 1.0])
            light = np.array([0.98, 0.98, 0.98, 1.0])
            arrow_color = np.where((luminance >= 0.52)[..., None], dark, light).reshape(-1, 4)
        ax.quiver(
            x_grid[::skip, ::skip],
            y_grid[::skip, ::skip],
            q1[::skip, ::skip],
            q2[::skip, ::skip],
            color=arrow_color,
            edgecolor=edgecolor,
            linewidth=linewidth,
            angles="xy",
            scale_units=scale_units,
            scale=scale,
            width=0.00065 * arrowsize,
            headwidth=3.0,
            headlength=4.2,
            headaxislength=3.7,
            pivot="mid",
            minlength=0.0,
        )
    axis1_label = axis1_name if axis1_name.endswith("*") else f"{axis1_name} (m)"
    axis2_label = axis2_name if axis2_name.endswith("*") else f"{axis2_name} (m)"
    if show_axis:
        ax.set_xlabel(axis1_label, fontsize=10)
        ax.set_ylabel(axis2_label, fontsize=10)
    if axis1_name.endswith("*"):
        ax.set_xlim(0.0, 1.0)
        ax.set_xticks((0.0, 1.0))
    else:
        # SLCF coordinates are physical mesh nodes.  pcolormesh with
        # same-sized coordinate/value arrays treats them as cell centers and
        # expands the outer polygons by half a spacing.  Clip those polygons
        # back to the actual FDS domain instead of showing, for example,
        # z=-0.1..7.1 for a z=0..7 mesh.
        ax.set_xlim(float(axis1[0]), float(axis1[-1]))
    if axis2_name.endswith("*"):
        ax.set_ylim(0.0, 1.0)
        ax.set_yticks((0.0, 1.0))
    else:
        ax.set_ylim(float(axis2[0]), float(axis2[-1]))
    # Keep one minor tick between adjacent major ticks.
    if axis1.size > 1:
        ax.xaxis.set_major_locator(LinearLocator(3))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    if axis2.size > 1:
        ax.yaxis.set_major_locator(LinearLocator(3))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(labelsize=9)
    ax.tick_params(which="minor", length=2.0, width=0.45)
    ax.set_aspect("equal", adjustable="box")
    if orientation == "horizon":
        divider = make_axes_locatable(ax)
        # A long tunnel slice keeps an equal data aspect, so the plotting
        # axes can become only a few pixels high.  A percentage-sized colorbar
        # then collapses into a black line.  Use a fixed physical thickness
        # so the gradient remains visible for long, shallow domains.
        cax = divider.append_axes("top", size=0.18, pad=0.24)
        cbar = fig.colorbar(image, cax=cax, ticks=ticks, orientation="horizontal")
        cbar.ax.xaxis.set_ticks_position("top")
        cbar.ax.xaxis.set_label_position("top")
        cbar.ax.tick_params(labelsize=8, pad=2)
        # Anchor the title to the colorbar itself.  This keeps the title a
        # stable, small distance above the bar for very shallow tunnel
        # slices, independent of the figure's data aspect ratio.
        cax.set_title(title, fontsize=11, pad=8)
    else:
        cbar = fig.colorbar(image, ax=ax, ticks=ticks)
        fig.suptitle(title, fontsize=11, y=0.96)
    if ticks is not None:
        if orientation == "horizon":
            cbar.ax.set_xticklabels([fmt_tick(tick) for tick in ticks])
        else:
            cbar.ax.set_yticklabels([fmt_tick(tick) for tick in ticks])
    # The quantity name is already shown once as the figure title. Do not
    # duplicate it beside or below the colorbar.
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    return vmin, vmax


def main() -> int:
    args = parse_args()
    args.smv = str(resolve_smv_path(args.smv))
    sim = fdsreader.Simulation(args.smv)
    slc, plane_value = select_slice(sim, args.quantity, args.plane_dir, args.plane_value)
    subslice = slc.subslices[0]
    coords = slice_coordinates(slc)
    axis_names = tuple(dim for dim in ("x", "y", "z") if dim != args.plane_dir)
    axis1_name, axis2_name = axis_names
    axis1 = coords[axis1_name]
    axis2 = coords[axis2_name]
    requested_end = None
    if args.time is None:
        requested_end = args.t1
    elif len(args.time) == 1:
        requested_end = args.time[0]
    else:
        requested_end = args.time[1]

    if args.time is not None and (args.t0 is not None or args.t1 is not None):
        raise RuntimeError("Use either --time T/--time T0 T1, or legacy --t0/--t1, not both.")

    if args.time is None:
        mode = "timeavg"
        load_mode = "legacy_timeavg"
    elif len(args.time) == 1:
        mode = "instant"
        load_mode = "instant"
    else:
        mode = "timeavg"
        load_mode = "timeavg"

    values, selected_times, frame_index, times, direct_info = load_slice_values(
        slc, load_mode, args, requested_end
    )

    clean_info = None
    if slc.quantity.name.upper() == "SOOT DENSITY" and not args.no_soot_clean:
        values, clean_info = clean_soot_density(values)

    plot_values = values
    obst_info = None
    if args.obst:
        obst_mask = None
        if args.obst_source in ("auto", "geometry"):
            obst_mask, obst_fds_path = geometry_obst_mask(
                args.smv, args.plane_dir, plane_value, axis1_name, axis1, axis2_name, axis2
            )
            if obst_mask is not None:
                obst_info = {
                    "source": f"geometry {obst_fds_path.name}",
                    "masked": int(np.count_nonzero(obst_mask)),
                    "eps": None,
                }
        if obst_mask is None:
            if args.obst_source == "geometry":
                raise RuntimeError("Cannot build OBST mask from input-file geometry.")
            if slc.quantity.name.upper() == args.obst_quantity.upper():
                obst_values = values
                obst_slice_name = slc.quantity.name
            else:
                obst_slc, obst_plane_value = select_slice(sim, args.obst_quantity, args.plane_dir, args.plane_value)
                if not math.isclose(obst_plane_value, plane_value, rel_tol=1.0e-9, abs_tol=1.0e-9):
                    raise RuntimeError(
                        f"OBST mask plane mismatch: {args.obst_quantity} at {obst_plane_value:g} m, "
                        f"current slice at {plane_value:g} m"
                    )
                obst_values, _, _, _, _ = load_slice_values(
                    obst_slc, load_mode, args, requested_end
                )
                obst_slice_name = obst_slc.quantity.name
            if obst_values.shape != values.shape:
                raise RuntimeError(
                    f"OBST mask shape {obst_values.shape} does not match plotted field shape {values.shape}"
                )
            obst_mask = np.abs(obst_values) <= args.obst_eps
            obst_info = {
                "source": f"{obst_slice_name} <= {args.obst_eps:.8g}",
                "masked": int(np.count_nonzero(obst_mask)),
                "eps": args.obst_eps,
            }
        plot_values = np.ma.masked_where(obst_mask, values)

    normalization_info = None
    if args.normalize:
        plot_values, normalization_info = normalize_field(plot_values)

    vector_overlay = None
    vector_info = None
    vector_normalization_scale = None
    csv_extra_fields = []
    if args.vector:
        vec1, vec2, vector_direct_infos = load_vector_overlay(
            sim, args, plane_value, load_mode, requested_end, values.shape,
            subslice.filename, base_slice=slc
        )
        csv_vec1 = np.array(vec1, dtype=float, copy=True)
        csv_vec2 = np.array(vec2, dtype=float, copy=True)
        if args.normalize:
            vec1, vec2, vector_normalization_scale = normalize_vector_components(vec1, vec2)
        if args.vector_normalize:
            vec1, vec2 = normalize_vector_directions(vec1, vec2)
        component_names = vector_component_names(args.plane_dir)
        csv_extra_fields.extend(((component_names[0], csv_vec1), (component_names[1], csv_vec2)))
        vector_mask = np.zeros(values.shape, dtype=bool)
        if np.ma.isMaskedArray(plot_values):
            vector_mask = np.ma.getmaskarray(plot_values)
        vector_overlay = (
            vec1,
            vec2,
            args.skip,
            args.scale,
            args.arrowsize,
            args.vector_color,
            args.vector_edgecolor,
            args.vector_linewidth,
            "width" if args.vector_normalize else "xy",
            vector_mask,
        )
        vector_info = {
            "components": component_names,
            "direct_infos": vector_direct_infos,
        }

    safe_quantity = args.quantity.strip().lower().replace(" ", "_")
    safe_value = f"{plane_value:g}".replace("-", "m").replace(".", "p")
    safe_plane = f"Plane{args.plane_dir}{safe_value}"
    if mode == "instant":
        safe_time = fmt_time(args.time[0]).replace("-", "m").replace(".", "p")
        default_prefix = f"{sim.chid}_{safe_quantity}_{safe_plane}_t{safe_time}"
    elif args.time is not None and len(args.time) == 2:
        safe_t0 = fmt_time(args.time[0]).replace("-", "m").replace(".", "p")
        safe_t1 = fmt_time(args.time[1]).replace("-", "m").replace(".", "p")
        default_prefix = f"{sim.chid}_{safe_quantity}_{safe_plane}_t{safe_t0}_{safe_t1}_avg"
    elif args.t0 is not None or args.t1 is not None:
        safe_t0 = fmt_time(args.t0) if args.t0 is not None else fmt_time(selected_times[0])
        safe_t1 = fmt_time(args.t1) if args.t1 is not None else fmt_time(selected_times[-1])
        safe_t0 = safe_t0.replace("-", "m").replace(".", "p")
        safe_t1 = safe_t1.replace("-", "m").replace(".", "p")
        default_prefix = f"{sim.chid}_{safe_quantity}_{safe_plane}_t{safe_t0}_{safe_t1}_avg"
    else:
        default_prefix = f"{sim.chid}_{safe_quantity}_{safe_plane}_timeavg"
    if args.normalize and args.output_prefix is None:
        default_prefix = f"{default_prefix}_normalized"
    prefix = args.output_prefix or default_prefix
    png_path = Path(f"{prefix}.png")
    csv_path = Path(f"{prefix}.csv")

    unit = slc.quantity.unit
    # Keep the figure title compact. Detailed plane/time information is
    # reported in the terminal and encoded in the output filename.
    title = slc.quantity.name.upper()
    # Keep all figure labels concise; normalization details stay in the log.
    colorbar_label = slc.quantity.name.upper()

    plot_axis1 = normalize_axis_coordinates(axis1) if args.normalize_coordinates else axis1
    plot_axis2 = normalize_axis_coordinates(axis2) if args.normalize_coordinates else axis2
    plot_axis1_name = f"{axis1_name}*" if args.normalize_coordinates else axis1_name
    plot_axis2_name = f"{axis2_name}*" if args.normalize_coordinates else axis2_name

    write_csv(csv_path, plot_axis1_name, plot_axis1, plot_axis2_name, plot_axis2, plot_values, csv_extra_fields)
    color_min, color_max = write_plot(
        png_path,
        plot_axis1_name,
        plot_axis1,
        plot_axis2_name,
        plot_axis2,
        plot_values,
        title,
        colorbar_label,
        args.cmap,
        args.plot_orientation,
        tuple(args.plot_range) if args.plot_range is not None else None,
        args.smooth,
        args.dpi,
        vector_overlay,
        normalization_info["mode"] if normalization_info is not None else None,
        args.show_axis,
    )

    print(f"SMV: {args.smv}")
    if len(slc.subslices) == 1:
        print(f"Slice file: {subslice.filename}")
    else:
        print(f"Slice files: merged {len(slc.subslices)} mesh subslices")
    print(f"Quantity: {slc.quantity.name} [{unit}]")
    print(f"Plane: {args.plane_dir} = {plane_value:g} m")
    if mode == "instant":
        print("Mode: instantaneous")
        print(f"Requested time: {args.time[0]:.8g} s")
        print(f"Nearest frame: {frame_index} at {selected_times[0]:.8g} s")
    else:
        print("Mode: time average")
        print(f"Frames used: {selected_times.size} / {times.size}")
        print(f"Time range used: {selected_times[0]:.8g} s to {selected_times[-1]:.8g} s")
    if clean_info is not None:
        print(
            "SOOT DENSITY clean: "
            f"reset {clean_info['reset']} values above {clean_info['cutoff']:.8g} to 0 "
            f"(raw max {clean_info['vmax']:.8g})"
        )
    if normalization_info is not None:
        print(
            "Field normalization: on, "
            f"{normalization_info['mode']}, scale={normalization_info['scale']:.8g}, "
            f"raw range=[{normalization_info['raw_min']:.8g}, "
            f"{normalization_info['raw_max']:.8g}]"
        )
    else:
        print("Field normalization: off")
    print(
        "Coordinate normalization: "
        + ("on, plotted axes in [0,1]" if args.normalize_coordinates else "off")
    )
    if direct_info is not None:
        print(f"SLCF read: {direct_info}")
    if vector_info is not None:
        print(
            "Vector: on, "
            f"{vector_info['components'][0]}, {vector_info['components'][1]}, "
            f"skip={args.skip}, scale={args.scale:g}, arrowsize={args.arrowsize:g}, "
            f"linewidth={args.vector_linewidth:g}, direction-normalized="
            f"{'yes' if args.vector_normalize else 'no'}"
        )
        for info in vector_info["direct_infos"]:
            print(f"Vector SLCF read: {info}")
        if vector_normalization_scale is not None:
            print(f"Vector normalization: signed [-1,1], scale={vector_normalization_scale:.8g}")
    else:
        print("Vector: off")
    if obst_info is not None:
        print(
            "OBST mask: "
            f"masked {obst_info['masked']} cells using {obst_info['source']}"
        )
    print(f"Grid: {axis1_name}={axis1.size}, {axis2_name}={axis2.size}")
    value_array = plot_values.compressed() if np.ma.isMaskedArray(plot_values) else np.asarray(plot_values).ravel()
    if value_array.size:
        print(f"Value min/max: {np.nanmin(value_array):.8g}, {np.nanmax(value_array):.8g}")
    range_mode = "fixed" if args.plot_range is not None else "auto"
    print(f"Color range: {color_min:g} to {color_max:g} ({range_mode})")
    print("Smooth: " + ("on, bilinear interpolation" if args.smooth else "off"))
    print(f"DPI: {args.dpi}")
    print(f"PNG: {png_path}")
    print(f"CSV: {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
