# GFDS

Version: 1.0

GFDS is a CUDA Fortran branch of the Fire Dynamics Simulator (FDS) for NVIDIA GPUs. It calculates three-dimensional low-Mach-number flow, multi-species transport, combustion, thermal radiation, and wall heat transfer. It also supports multi-mesh calculations on one or more GPUs.

Input files and configuration examples are provided in [`case/`](case/).

## Build

GFDS requires the NVIDIA HPC SDK with `nvfortran`, CUDA, an OpenMPI installation compatible with NVHPC, CMake 3.20 or newer, and GNU Make.

Double precision is the default:

```bash
cmake -S . -B build
cmake --build build -j8
```

To specify the OpenMPI installation directory explicitly:

```bash
cmake -S . -B build -DNVHPC_MPI_ROOT=/path/to/openmpi
cmake --build build -j8
```

For single precision, use a separate build directory:

```bash
cmake -S . -B build-single -DSINGLE_PRECISION=ON
cmake --build build-single -j8
```

The resulting executable is `build/fds`.

## Run

Run a case in GPU mode:

```bash
mpirun -np 1 build/fds "case/case-3-Billaud/RA6.fds" --g
```

The execution modes are:

| Option | Mode |
| --- | --- |
| `--g` or `--gpu` | GPU mode and the default execution mode |
| `--gc` | GFDS call path using the corresponding CPU numerical functions for comparison and module timing |
| `--c` or `--cpu` | Original FDS CPU time-step path |

A one-to-one mapping between MPI processes and GPUs is recommended. For example, run a four-mesh case with four MPI processes and four GPUs:

```bash
mpirun -np 4 build/fds "case/case-7-Speedup/Arup Tunnel/4MESH/T_100_4.fds" --g
```

Fewer MPI processes can also run a multi-mesh case when fewer GPUs are available. For example, one MPI process can run all four meshes on one GPU:

```bash
mpirun -np 1 build/fds "case/case-7-Speedup/Arup Tunnel/4MESH/T_100_4.fds" --g
```

## Python Plotting

[`GUtilities/slice_contour.py`](GUtilities/slice_contour.py) extracts FDS slice fields and generates contour or vector figures together with the plotted CSV data. Installation and example commands are provided in the [Python plotting guide](GUtilities/USAGE.md). The complete command reference is available from:

```bash
python GUtilities/slice_contour.py --help
```

## Contact

- [xiaoxm8@mail2.sysu.edu.cn](mailto:xiaoxm8@mail2.sysu.edu.cn)
- [mingpj@mail.sysu.edu.cn](mailto:mingpj@mail.sysu.edu.cn)

## Citation and License

Citation metadata are provided in [CITATION.cff](CITATION.cff). The original FDS/NIST license notice is retained in [LICENSE.md](LICENSE.md), and the GFDS modification notice is provided in [NOTICE.md](NOTICE.md).

GFDS makes no warranty of any kind. Users assume all risks arising from its use.

## Acknowledgment

Thanks to Marcos Vanella for feedback from testing GFDS on different GPU platforms.
