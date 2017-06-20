# RandEM: Random Encounter Modelling in R

A package for estimating population densities and relevant parameters using Random Encounter Models (REMs).

Random Encounter Models allow the estimation of animal population densities and other ecological measures from camera trap data, for a given survey area. The REM can be applied equally to species which exhibit individually identifiable markings and those which do not. The method has been found to perform well when compared to other density estimaton techniques such as Distance sampling and, due to the nature of the data collection medium (i.e. camera traps) may be less demanding in terms of survey effort and allow for greater spatio-temporal analyses of population dynamics.

Here, we aim to produce a comprehensive toolkit for REM analyses in R. 

We encourage anyone interested in contributing to get in touch.


# Installation


```r
devtools::install_github("arcaravaggi/RandEM")
```


# Road map

Includes

* Effective Detection Distance (EDD)
* REMs
* Survey design tools
* Web app to calculate densities via REM


## Contributions

Contributors are welcome to fork the package and suggest additions or improvements. If you would like to be included as a contributor, please let us know.

#### I found a bug

Please report it to the [issue tracker][issues]. Please provide specific details, allowing the error to be reproduced and investigated. Always note the version of R you are using, along with any other relevant software (e.g. RStudio).  

[issues]: https://github.com/arcaravaggi/RandEM/issues

# License

This package is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

A copy of the GNU General Public License, version 3, is available at https://www.r-project.org/Licenses/GPL-
