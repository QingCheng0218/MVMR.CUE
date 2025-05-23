MVMR.CUE
=======
  
  **MVMR.CUE** is a package for Robust multivariant Mendelian Randomization method to account for Correlated and Uncorrelated Pleiotropy.

Installation
============
  Install the development version of **MVMR.CUE** by use of the 'devtools' package. Note that **MVMR.CUE** depends on the 'Rcpp' and 'RcppArmadillo' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.
```
library(devtools)
install_github("QingCheng0218/MVMR.CUE@main")
```

If you have errors installing this package on Linux, try the following commands in R:
  ```
Sys.setlocale(category = "LC_ALL", locale = "English_United States.1252") 
Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

library(devtools)
install_github("QingCheng0218/MVMR.CUE@main")
```

Examples
=========
```
library(MVMR.CUE);
# Run MVMR.CUE using independent GWAS samples.
mvmr_cue = MVMRCUEIndepSample(gamma1, Gamma2, se1, se2);

beta1hat = mvmr_cue$beta.hats[1];
beta1hat.pvalue = mvmr_cue$beta.p.values[1];
beta2hat = mvmr_cue$beta.hats[2];
beta2hat.pvalue = mvmr_cue$beta.p.values[2];
```


Development
===========
  
  This package is developed and maintained by Qing Cheng (qingcheng0218@gmail.com). 
