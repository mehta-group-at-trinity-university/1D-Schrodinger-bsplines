CalcBoundState.x:	CalcBoundState.f90 Bsplines.f matrix_stuff.f CalcBoundState.o Bsplines.o matrix_stuff.o nrmat.o TwoBodydata.o ~/bin/minpack.o TwoBodydata.mod
	g95    -ftrace=full -freal-loops CalcBoundState.o matrix_stuff.o Bsplines.o nrmat.o ~/bin/minpack.o TwoBodydata.o -framework vecLib -L/opt/local/lib/ -I/opt/local/lib/ -L/users/mehtan/bin -larpack -o CalcBoundState.x

CalcBoundState.o:	CalcBoundState.f90
	g95    -ffixed-line-length-132 -freal-loops -ftrace=full -c CalcBoundState.f90

matrix_stuff.o:	matrix_stuff.f
	g95    -ffixed-line-length-132 -ftrace=full -c matrix_stuff.f

Bsplines.o:	Bsplines.f
	g95    -ffixed-line-length-132 -ftrace=full -c Bsplines.f	

nrmat.o:	nrmat.f
	g95    -ffixed-line-length-132 -ftrace=full -c nrmat.f

TwoBodydata.o:   TwoBodydata.f90
	g95    -ffixed-line-length-132 -ftrace=full -c TwoBodydata.f90 



       