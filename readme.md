# Quantum tomography
This repository contains the matlab functions and scripts used to generate the figures in *Quantum process tomography via completely positive and trace-preserving projection*, [Physical Review A 98, 062336 (2018)](https://doi.org/10.1103/PhysRevA.98.062336).

### Dependencies
Many of the solvers are called via the excellent [`cvx` package](https://www.cvxr.com). In order to run them, install `cvx` in the parent directory. A cvx license file is necessary to run certain solvers, such as `mosek`.

## Reproducing plots
`do_paper_plots.m` is a script that sets certain global parameters:

- `ensemble`- either 'qp' (for quasi-pure) or 'fr' (for full rank);
- `drange` - a list of Hilbert space dimesions to iterate over;
- `LIswitch`- 0 to skip linear inversion, 1 to include linear inversion;
- `ensemble_size` - number of ground-truth processes to generate, for each Hilbert space dimension;
- `Npows` - a list of powers of ten to iterate over (representing the number of multinomial trials performed). As `Npow` increases, shot noise is decreased,

before generating a simulated dataset, running each of a number of tomography solvers on that dataset, and producing plots that summarise the runtime and accuracy of those solvers. The datasets and statistics are saved to disk so that (for example) the plots can be altered without having to re-run all of the code, or the solvers can be tweaked and run on the same dataset (if desired).

## CPTP projection
The main contribution of this work is the CPTP projection algorithm, which accepts a `d^4 x 1`, vectorised Choi matrix (in the column stacking convention), and returns a vectorised *projected* Choi matrix, which represents the closest completely positive, trace-preserving process to the input Choi matrix. For performance reasons, the function also accepts some precomputed helper matrices which are used in the subroutine that projects onto trace-preserving matrices. The function signature is

```matlab
function [ projected_choi_vec ] = CPTP_project( choi_vec, MdagM, Mdagb  )
```
The user may alter the exit criterion defined in `while` loop: 

```matlab
while GAP(end) >= 1e-4
```
in order to alter the tradeoff between runtime and accuracy.