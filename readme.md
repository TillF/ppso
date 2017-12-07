Particle Swarm Optimization and Dynamically Dimensioned Search, optionally using parallel computing based on Rmpi
======================================================================

**This is a fork of ppso (version 0.9994) from https://www.rforge.net/ppso/**

**Originally written and maintained by Till Franke**

## Description
(Optionally parallelized) optimization using PSO (Particle Swarm Optimzation) or DDS (Dynamically Dimensioned Search) algorithms, which excell for multidimensional (multi-parameter) functions with many local extrema and a restricted number of function evaluations.
The parallelized version builds on Rmpi and is intended for highly computationally intensive objective functions (>20 s evaluation time).
Another focus of this package is the possibility to resume interrupted optimization runs from the intermediate project files.
These features make this package useful for the automatic calibration of complex numerical models (e.g. hydrological models). 

For bug fixes, comments or further development please open an issue or contact: franke@uni-potsdam.de or mreich@posteo.de.

## Installation

1. Start R
2. Install package via devtools: 
`devtools::install_github("marcianito/ppso")`

3. load packages: 
`library(ppso)`;

## Dependencies

### Computationally
* r-base version 3.3.1
* following R-packages: devtools,  ??
* system libraries for devtools

in debian install using: 
`apt-get update && apt-get install libssl-dev`

