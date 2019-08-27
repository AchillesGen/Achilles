# Notes on interfacing Fortran and Python

## Basic steps

1. Create pyf files that define the interface for python and fortran, by running:
```
f2py input_files -m module_name -h output.pyf
```
Where the input files are all the files that contain code that needs to be called from the python code
2. After the file is generated, go into the output.pyf and update the function signatures defining what is input and what is output through the use of "intent(in/out)"
3. Compile the fortran code into a library that python can then use by running:
```
f2py -c output.pyf files.f90
```
Where the files that are listed are all the files needed to get the code to run. In other words, at this point you are creating a library, so any function that is called internally needs to be included.

## Example for Noemi's code:
```
f2py xsec_fact_new.f90 currents_opt_v1.f90 -m xsec_fact -h xsec_fact.pyf
```

To see the modifications, run the above command with something else instead of xsec\_fact.pyf and compare the files.

```
f2py -c xsec_fact.pyf xsec_fact_new.f90 currents_opt_v1.f90 nform.f90
```

This would then create a library called something like: xsec\_fact.cpython-36m-x86\_64-linux-gnu.so (this depends on the architecture of your machine)
Included is a file xsec.py which initializes all the needed variables, and evaluates the function at a given point. This has been compared to the fortran output, and is consistent to an accuracy of 1e-11.
