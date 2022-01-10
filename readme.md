R-package ppso: Particle Swarm Optimization and Dynamically Dimensioned Search, optionally using parallel computing based on Rmpi
======================================================================

**This is a replacement of the repository formerly maintained at http://rforge.net/ppso/ (version 0.9994)**

## Description
(Optionally parallelized) optimization using PSO (Particle Swarm Optimzation) or DDS (Dynamically Dimensioned Search) algorithms, which excel for multidimensional (i.e. multi-parameter) functions with many local extrema and a restricted number of function evaluations.
The parallelized version builds on Rmpi and is intended for highly computationally intensive objective functions (>20 s evaluation time).
Another focus of this package is the possibility to resume interrupted optimization runs from the intermediate project files.
These features make this package useful for the automatic calibration of complex numerical models (e.g. hydrological models). 

For bug fixes, comments or further development please [open an issue here at github](https://github.com/TillF/ppso/issues).

## Installation
1. start R
2. run these commands:
```
install.packages("devtools")
devtools::install_github("TillF/ppso")
```
 If you are unable to get RTools and devtools running, you may also download the files (green button "Code" at top right, "Download ZIP") and proceed following [these instructions](https://stackoverflow.com/questions/30989027/how-to-install-a-package-from-a-download-zip-file/30989367).
 
## Dependencies
basic functionality: no further dependencies; optionally R-packages `rgl` (visualization)), `lhs` (improved initialization of starting values)
### for installation only
* R-packages: `devtools`
* system libraries for `devtools`
in debian, install using: 
`apt-get update && apt-get install libssl-dev`;

In Windows, this requires RTools, which will be installed along.
### in parallel mode (for computationally-demanding functions, e.g. model calibration)
* R-package: `Rmpi` (may inquire further system installations, please check [the Rmpi pages](http://www.stats.uwo.ca/faculty/yu/Rmpi/).



