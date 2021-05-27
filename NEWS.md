# BigTime 0.2.0

* `sparseVAR`, `sparseVARMA`, `sparseVARX` are all more efficient and faster now
* All estimation functions return corresponding S3 classes
* `fitted` and `residuals` functions were implemented for VAR, VARMA, VARX 
* `diagnostics_plot` was implemented for VAR, VARMA, VARX
* `plot_cv` was implemented; Allows to investigate the behaviour in CV 
* `sparseVAR` can now use information criteria to select the optimal penalisation
* `simVAR` allow for the simulation of VAR models using different sparsity patterns
* Default selection procedure in `sparseVAR` changed to "none". Will return a 3D array of estimations from this version on, unless a different selection procedure is chosen. **WARNING: This can break old code!**; utitlity functions `summary` and `plot` were also implemented for simulated VARs
* `recursiveforecast` was implemented for VAR models. 
* `is.stable` was implemented for VAR models
* in `plot_cv`, `directforecast`, and `lagmatrix` the model argument was removed. Functions know  automatically what model was used. **WARNING: This can break old code!** 
