# swOpenFOAM
The accelerated version of OpenFOAM on the Sunway TaihuLight supercomputer 


## 1. How to compile
  * 1.1 edit variable `foamInstall` in `OpenFOAM-3.0.0/etc/bashrc` and `OpenFOAM-3.0.0/etc/bashrc_run`
    * `OpenFOAM-3.0.0/etc/bashrc foamInstall`: using content of `echo $PWD`
    * `OpenFOAM-3.0.0/etc/bashrc_run foamInstall`: using content of `pwd`
  * 1.2 source `OpenFOAM-3.0.0/basrc_compile`: compiler settings and parallel compiling env are loaded      

## 2. How to run
  * 2.1 source `OpenFOAM-3.0.0/basrc_run`
