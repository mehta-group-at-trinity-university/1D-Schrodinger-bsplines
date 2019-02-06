rm *.o *.x
g95 -c TwoBodydata.f90
g95 -c modules_qd.f90
make -f CalcBoundState.mak
