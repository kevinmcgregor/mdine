# mdine
(**Note: This package is currently in development**) **M**icrobiome **Di**fferential **N**etwork **E**stimation (**mdine**) allows the estimation of OTU co-occurrence networks within two separate groups, where the networks are defined through precision matrices.  The difference between the two precision matrices is also estimated, along with corresponding interval estimates.

## Installation
The first step to installing **mdine** is to install [rstan](http://mc-stan.org/users/interfaces/rstan) along with the appropriate compiler.  The steps on how to do this are different based on your OS:

### Installing **mdine** for Mac/Linux
* Install **rstan** by following the [Mac/Linux instructions](https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Mac-or-Linux)
  * Make sure you have a C++ toolchain and configuration installed as described in the instructions.

Once **rstan** has been successfully installed, run the following code to install **mdine**:
* Run the following code
```r
if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}
Sys.setenv(USE_CXX14=1)
install_github("kevinmcgregor/mdine", dependencies=TRUE)
```


### Installing **mdine** Windows
* Install **rstan** by following the [Windows instructions](https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Windows)
  * Make sure to install **rtools** as described in the instructions.

Once **rstan** has been successfully installed, run the following code to install **mdine**:
* Run the following code
```r
if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}
install_github("kevinmcgregor/mdine", dependencies=TRUE)
```
