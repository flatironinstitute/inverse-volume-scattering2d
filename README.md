This repository contains beta MATLAB codes for data generation of 
forward medium problems using a hybrid HPS-HBS solver, 
and an inverse medium solver based on the recursive linearization
approach. 

The code depends on the MATLAB installation of
[FINUFFT](https://github.com/flatironinstitute/finufft.git), and
[fmm2d](https://github.com/flatironinstitute/fmm2d.git). 
Please ensure that the matlab folder of finufft, and fmm2d containing the matlab `.mex*' 
file is included in the MATLAB path.

Forward medium solver
----------------------------
To run the solver in forward mode, run `generating_data.m` in
the examples directory. You can specify the domain by setting the number
of coefficients in the sine series `N` on line 33, and then setting the 
`(N,N) coefs` matrix between lines 42-44. 

The data is currently stored in `test_data.mat` in the examples folder.
This can be changed via changing the variable `filename` on line 73. 

The file stores the following quantities:
* The list of frequencies used for generating the data `khv`
* Number of incident directions at each frequency denoted by
  `ntheta(nkh)`
* Number of sensors correspdoning to each incident frequency
  `npoints(nkh)`
* Radius of the disc on which the sensors are located `radius`
* `u_meas(nkh)` a struct containing a field called `field` which
  stores the scattered data 

Here `nkh` is the number of frequencies or length(khv)

The parameters `khv, ntheta, npoints`, and radius can be set by editing
the block of lines 25-40 in `parameters_rla.m`


Inverse medium solver
----------------------------
To run the solver in inverse mode, run `driver_volume.m` in the examples
directory. The input to this script is the `test_data.mat` generated
by the `generating_data.m` script, and relies on the fields
`khv, u_meas, ntheta, npoints, radius` in the `.mat` file

The script runs recursive linearization to construct the inverse
medium properties, and stores the results in `test_result.mat`
containing an array of structs called `solution` of size `nkh`
which has the following fields

* `q_newton`: point values of the medium at a grid of points `XP, YP`
  set via `parameters_volume.m`
* `coefs`: the sine series coefficients corresponding to the solution
* `rhs`: the residual vector from the measurement data
* `rel_rhs`: relative norm of the residual vector
* `it`: number of newton iterations required at the frequency
* `stop`: The stopping criterion which terminated the optimization loop
* `lsqr`: max number of lsqr iterations used at any iterate


Codes by: Carlos Borges, Adrianna Gillman
