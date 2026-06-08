# AdSamp — Advanced Sampling for Rare Events

AdSamp is a C++ software package for running Forward Flux Sampling (FFS) simulations on high performance computing clusters. It supports both LAMMPS and GROMACS as molecular dynamics engines and is designed to run scalably on scavenge partitions, requeueing automatically after interruptions and running until completion without further user input.

It has been used to study protein folding (Trp-Cage, Hen Egg White Lysozyme), ion transport (first passage through a potassium channel), and conformational change (Adenylate Kinase compaction).

## How It Works

FFS estimates the rate of rare transitions by breaking them into a sequence of more achievable milestones. A scalar order parameter λ tracks progress along the event. The overall rate is:

$$T = \frac{\sum^{N}_{i=0} T_i/C_i}{\prod^{M}_{i=0} P(\lambda_{i+1}|\lambda_i)}$$

A full simulation consists of three phases run in sequence:

1. **Exploration** — relax input structures and characterize order parameter behavior
2. **Basin** — measure the flux of trajectories crossing from λ_A to λ_0
3. **Probabilistic steps** — repeated milestoning from basin outputs until a full transition is achieved

See [Algorithms](docs/algorithms.md) for a detailed description of each phase.


## Dependencies

AdSamp can be compiled in either **LAMMPS mode** or **No LAMMPS mode**.

### LAMMPS mode
Requires:

- LAMMPS patched with PLUMED ([LAMMPS Installation](LAMMPSinstallation.md))
- PLUMED
- OpenMPI
- FFTW

### No LAMMPS mode
Requires:

- Custom GROMACS patched with PLUMED [Custom Gromacs Installation](ModifiedGromacs/README.md) 
- OpenMPI
- FFTW

> No LAMMPS mode does NOT require LAMMPS to be installed. This can only run GROMACS simulations

> LAMMPS mode does NOT require GROMACS to be installed. This can run either GROMACS (if GROMACS is also installed, package is not used in compilation) or LAMMPS simulations. 

> AdSamp requires a modified build of GROMACS that extends its internal interruption mechanism to save checkpoints immediately upon any interruption signal. This is required for accurate FFS saves and provides a significant speedup over frequent periodic saves. See [MD Engines](docs/MDengines.md) for details.

> See [Custom Gromacs Installation](ModifiedGromacs/README.md) for installation information.

## Installation

AdSamp can be compiled in two different modes depending on the MD engine being used.

### GROMACS build

For GROMACS simulations, compile with:

```bash
mpicxx -std=c++17 AdSamp.cpp orderparameterpicker.cpp -o AdSamp -I.
```

### LAMMPS build

For LAMMPS simulations, compile with:

```bash
mpicxx -std=c++17 -DUSE_LAMMPS AdSamp.cpp orderparameterpicker.cpp -o AdSamp -I. -I/pathtolammps/include -L/pathtolammps/lib64 -llammps -Wl,-rpath,/pathtolammps/lib64
```


All further interactions with the software use:

```plaintext
./AdSamp
```



## Quick Start

AdSamp is controlled through a single binary that operates in different modes depending on how it is invoked. The typical workflow for one simulation phase is:

1. Run `./ManuelsAdvancedSampling` with no arguments to enter the login flow
2. Provide run parameters when prompted — inputs are validated before any cluster resources are used
3. Submit the generated `AdSampControl.sh` to your scheduler
4. Monitor progress via `Runinfo.txt` in the working directory, or wait for the email notification
5. Once complete, re-run `./ManuelsAdvancedSampling` and select postprocessing
6. Submit `PostProcessing.sh`, then repeat from step 1 for the next phase

See [Running a Full Simulation](docs/MultiStepWorkflow.md) for the complete multi-step workflow.


## Input Files

Before running the login flow, ensure the following are in place:

- Structure files in a single directory, named with integers and a shared delimiter (e.g. `1.gro`, `2.gro`)
- `ffs.mdp` in the working directory with all required fields for your MD engine
- `plumed.dat` in the working directory with a `COMMITTOR`, `PRINT`, and `FLUSH` section

See [Input Files](docs/Inputfiles.md) for full requirements and examples.


## Tutorials

There are three included tutorials AdSamp, each of these provides examples for running a full set of simulations [Multistep](docs/multistep.md) (involving exploration, basin, probabilistic, postprocessing, and binary export examples [Algorithms](docs/algorithms.md). The tutorials are as follows:

- [GROMACS](/tutorials/GROMACS/)
- [LAMMPSwithPLUMED](/tutorials/LAMMPSwithPLUMED)
- [LAMMPSnoPLUMED](/tutorials/LAMMPSnoPLUMED)


## Documentation

| Page | Contents |
|---|---|
| [Algorithms](docs/algorithms.md) | Simulation types that can be specified when interacting with the script: exploration, basin, probabilistic, postprocessing, binary export |
| [Input Files](docs/Inputfiles.md) | Required input files and naming conventions |
| [Outputs](docs/outputfiles.md) | Output directory structure, bookkeeping files, postprocessing |
| [MD Engines](docs/mdengines.md) | LAMMPS vs GROMACS differences, command reference, modified GROMACS |
| [HPC Usage](docs/hpc.md) | Scavenge partitions, requeue behavior, live monitoring, tmux |
| [Running a Full Simulation](docs/multistep.md) | End-to-end FFS workflow across all phases |
| [Software Internals](docs/InternalSoftwareOperation.md) | Login mode, controller, simulation run mode, binary export |
| [Custom OP](docs/customorderparameter.md) | Creating and using a custom order parameter script for LAMMPS (no PLUMED) |


## Repository Structure

```plaintext
AdSamp/
-AdSamp.cpp
-OrderParameter.hpp
-orderparameterpicker.cpp / .h
-tutorials/
-templates/
  -
  - 
  -
-docs/
  - algorithms.md
  - Inputfiles.md
  - outputfiles.md
  - MDengines.md
  - hpc.md
  - multistep.md
  - InternalSoftwareOperation.md
  - customorderparameter.md
- assets
-README.md
```


## Authors

Manuel Alejandro López, Betül Uralcan, Amir Haji-Akbari
