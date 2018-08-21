# mdine
### Installation
In order to install **mdine** you'll first need to make sure to [rstan](http://mc-stan.org/users/interfaces/rstan) along with the appropriate compiler.  Instructions based on OS are given below:

### Mac/Linux
  * Install **rstan** by following these [instructions](https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Mac-or-Linux)
  * Run the following code
  ```r
  if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}
install_github("kevinmcgregor/mdine")
  ```

### Windows

