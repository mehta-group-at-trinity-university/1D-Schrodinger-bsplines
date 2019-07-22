rm *.o *.x
gfortran -c TwoBodydata.f90
gfortran -c modules_qd.f90
make -f CalcBoundState.mak
