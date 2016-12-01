# Diversity and Macro-Evolution Analysis

This package consists on a set of tools for Diversification process analysis. 

### Instalation. 

For instalation please 

 1. download the .zip file and save the containing folder. 
 2. open 'dmea.Rproj' with Rstudio
 3. Run the code above to be sure you have all dependences needed 
 
 ```{r}
 ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("ape","subplex","latex2exp","gridExtra","foreach","doParallel")
ipak(packages)
```

 3. under build menu, click on build and reload
 4. you should have the dmea package on your library now, type ```library('dmea')``` to check it.
 
 
 
 
### Some travis check (not using yet)
 
  [![Build Status](https://travis-ci.org/richelbilderbeek/dmea.svg?branch=master)](https://travis-ci.org/richelbilderbeek/dmea)
 [![codecov.io](https://codecov.io/github/richelbilderbeek/dmea/coverage.svg?branch=master)](https://codecov.io/github/richelbilderbeek/dmea?branch=master)
