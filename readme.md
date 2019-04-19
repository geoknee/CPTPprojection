# Quantum tomography
This repository contains the matlab functions and scripts used to generate the figures in *Quantum process tomography via completely positive and trace-preserving projection*, [Physical Review A 98, 062336 (2018)](https://doi.org/10.1103/PhysRevA.98.062336).

### Dependencies
Some of the benchmark solvers are called via the [`cvx` package](https://www.cvxr.com). In order to run them, you will need install `cvx` in the parent directory. A cvx license file is necessary to run certain solvers, such as `mosek`.

This code has been tested with MATLAB R2017a and CVX v2.1. 

## Reproducing plots


### Figures 1c, 3, 4, 7
These are generated using `fig1c.m`, `ROBUSTNESS_investigation.m`, `LI_typical_distance_to_CP.m`, and `how_often_cvx_fails.m`respectively.

### Figures 2, 5, 8

These are generated using `do_paper_plots.m`. It is a script that sets certain global parameters:

- `ensemble`- either 'qp' (for quasi-pure) or 'fr' (for full rank);
- `drange` - a list of Hilbert space dimesions to iterate over;
- `LIswitch`- 0 to skip linear inversion, 1 to include linear inversion;
- `ensemble_size` - number of ground-truth processes to generate, for each Hilbert space dimension;
- `Npows` - a list of powers of ten to iterate over (representing the number of multinomial trials performed). As `Npow` increases, shot noise is decreased,

before generating a simulated dataset, running each of a number of tomography solvers on that dataset, and producing plots that summarise the runtime and accuracy of those solvers. The datasets and statistics are saved to disk so that (for example) the plots can be altered without having to re-run all of the code, or the solvers can be tweaked and run on the same dataset (if desired).

NB: these figures take quite some time to generate. I recommend that you set the parameters above to 'easy' values initially, to make sure the code runs on your machine. When you are happy that things are functioning as you expect, you can then move on to the higher dimensions, larger ensemble sizes and so on.

The various solvers accept two arguments: 

- `A` a "prepare and measure" matrix describing the preparations and measurements used to generate the clicks;
- `n` list of clicks generated through a multinomial random number generator, with probabilities p = A * vec(C),

and return a vectorised Choi matrix.

Please see the paper for full definitions. 

*Please note that `pgd`
 is currently referred to as `gdap` in this repository.*

## CPTP projection
The main contribution of this work is the CPTP projection algorithm, which accepts a `d^4 x 1`, vectorised Choi matrix, and returns a vectorised *projected* Choi matrix, which represents the closest completely positive, trace-preserving process to the input Choi matrix. For performance reasons, the function also accepts some precomputed helper matrices which are used in the subroutine that projects onto trace-preserving matrices. The function signature is

```matlab
function [ projected_choi_vec ] = CPTP_project( choi_vec, MdagM, Mdagb  )
```
The user may alter the exit criterion defined in `while` loop: 

```matlab
while GAP(end) >= 1e-4
```
in order to alter the tradeoff between runtime and accuracy.

The CPTP projection is used in both the pgdB (projected gradient descent with backtracking) and LIFP (Liner Inversion with a Final Projection) solvers, for quantum process tomography. However, it is anticipated that the projection may be of wider use and interest.


### Vectorization convention
Vectorization assumes the following matlab convention: 

```matlab
vectorized_Choi = reshape(Choi,[],1);
```