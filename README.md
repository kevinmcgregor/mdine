# mdine

## Installation
The first step to installing **mdine** is to first install [rstan](http://mc-stan.org/users/interfaces/rstan) along with the appropriate compiler.  The steps on how to do this are different based on your OS:

### Mac/Linux
* Install **rstan** by following these [Mac/Linux instructions](https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Mac-or-Linux)

### Windows
* Install **rstan** by following these [Windows instructions](https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Windows)
  * Make sure to **rtools** as described in the instructions.

Once **rstan** has been successfully installed, run the following code to install **mdine**:
* Run the following code
```r
if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}
install_github("kevinmcgregor/mdine")
```
