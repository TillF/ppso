R-package ppso: Particle Swarm Optimization and Dynamically Dimensioned Search, optionally using parallel computing based on Rmpi
======================================================================

**This is a replacement of the repository formerly maintained at http://rforge.net/ppso/ (version 0.9994)**

## Description
(Optionally parallelized) optimization using PSO (Particle Swarm Optimzation) or DDS (Dynamically Dimensioned Search) algorithms, which excell for multidimensional (multi-parameter) functions with many local extrema and a restricted number of function evaluations.
The parallelized version builds on Rmpi and is intended for highly computationally intensive objective functions (>20 s evaluation time).
Another focus of this package is the possibility to resume interrupted optimization runs from the intermediate project files.
These features make this package useful for the automatic calibration of complex numerical models (e.g. hydrological models). 

For bug fixes, comments or further development please open an issue here at github.

## Installation

1. Start R
2. install package `devtools` and load it: `library(devtools)`
3. Install package `ppso` via `devtools`: 
`devtools::install_github("TillF/ppso")`
4. load package: 
`library(ppso)`;

## Dependencies
basic functionality: no further dependencies; optionally R-packages `rgl`, `lhs`
### for installation only
* R-packages: `devtools`
* system libraries for `devtools`
in debian install using: 
`apt-get update && apt-get install libssl-dev`;
in Windows may require RTools (I don't quite remember)
### in parallel mode (for computationally-demanding functions, e.g. model calibration)
* R-packages: `Rmpi` (may inquire further system installations, please check at http://www.stats.uwo.ca/faculty/yu/Rmpi/



