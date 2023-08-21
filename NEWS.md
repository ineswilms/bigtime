# bigtime 0.2.3
* package documentation adjusted due to breaking change in roxygen2 7.0.0 

# bigtime 0.2.2
* verbose arguments added to functions sparseVAR and sparseVARX
* correction lambda grid to hvaralgorithms
* removal SystemRequirements: C++11

# bigtime 0.2.1
* Dependency on Rcpp version >=1.0.7 added to the DESCRIPTION.

# bigtime 0.2.0

* Important starting note: Please update all R packages on which bigtime depends upon installing bigtime from CRAN. We have been informed on errors occurring when using older versions of the Rcpp package on which bigtime depends. In version 0.2.1 of bigtime (currently only available on github), we have therefore added the dependency on Rcpp version >=1.0.7
* `sparseVAR`, `sparseVARMA`, `sparseVARX` are all more efficient and faster now
* All estimation functions return corresponding S3 classes
* `fitted` and `residuals` functions were implemented for VAR, VARMA, VARX 
* `diagnostics_plot` was implemented for VAR, VARMA, VARX
* `plot_cv` was implemented; Allows to investigate the behavior in CV 
* `sparseVAR`, `sparseVARMA` and `sparseVARX` can now use information criteria to select the optimal penalization
* `simVAR` allow for the simulation of VAR models using different sparsity patterns; utility functions `summary` and `plot` were also implemented for simulated VARs
* Default selection procedure in `sparseVAR`, `sparseVARMA` and `sparseVARX` changed to "none". Will return a 3D array of estimations from this version on, unless a different selection procedure is chosen. **WARNING: This can break old code!**
* `recursiveforecast` was implemented for VAR models. 
* `is.stable` was implemented for VAR models
* in `plot_cv`, `directforecast`, and `lagmatrix` the model argument was removed. Functions know  automatically what model was used. **WARNING: This can break old code!** 
* All model functions will give warnings if the given data is not standardized
* added example data for VAR, VARMA, and VARX

# bigtime 0.1.0

* First release on CRAN
