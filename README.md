# MCRS with promiscous activities

We have created this code to test the effects of catalytic promiscuity at the origin of life. For more deatails about MCRS models see: https://doi.org/10.1016/j.jtbi.2015.06.002
Please note, that this repository does not contains all the simulations that has been published. The necessary modifications are commented out in the source code.
Also note, that comments and variable names are hungarian.

## Prerequisites

The code has been written on Debian. You will need:

  * build-essential
  * pkg-config
  * libgsl23 (https://www.gnu.org/software/gsl/)
 
 ## Installing
 
At the codes library: 
 
```
make
```

## Building the parameter file

The parameter file is needed to run the simulations in the declared parameter space. To set the parameters, number of repeats edit *generalas* file. To build parameter file press

```
./generalas
```

## Starting the simulations

```
./start.sh
```

## Evaluating the output

The results of the simulations are stored in binary files. To evaluate them, you should generate the reader by

```
make olvaso
```

For further help see the help of the reader

```
./olvaso --help
```
