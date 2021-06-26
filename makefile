# This is an commentary line in a makefile
# Start of the makefile
primo: nrtype.o nr.o nrutil.o bessel0.o levin.o
	gfortran -o primo *.o -llapack -lblas
	
nrtype.mod: nrtype.o nrtype.f90
	gfortran -c nrtype.f90 -llapack -lblas
nrtype.o: nrtype.f90
	gfortran -c nrtype.f90 -llapack -lblas
	
nr.mod: nr.o nr.f90
	gfortran -c nr.f90 -llapack -lblas
nr.o: nr.f90
	gfortran -c nr.f90 -llapack -lblas	
	
nrutil.mod: nrutil.o nrutil.f90
	gfortran -c nrutil.f90 -llapack -lblas
nrutil.o: nrutil.f90
	gfortran -c nrutil.f90 -llapack -lblas	
	
bessel0.mod: bessel0.o bessel0.f90
	gfortran -c bessel0.f90 -llapack -lblas
bessel0.o: bessel0.f90
	gfortran -c bessel0.f90 -llapack -lblas

	

	

	
solve.mod: solve.o solve.f90
	gfortran -c solve.f90 -llapack -lblas
solve.o: solve.f90
	gfortran -c solve.f90 -llapack -lblas

	

levin.o : nrtype.mod nr.mod nrutil.mod bessel0.mod solve.mod levin.f90
	gfortran -c levin.f90 -llapack -lblas

cleanall:
	rm -r *.o primo *.mod *.pdf *.png *.dat *.log *.txt *.TXT .fuse* ../data ../plot ../errorlog

clean:
	rm -r *.o primo *.mod *.pdf *.png *.dat *.log *.txt *.TXT .fuse* ../data

cleanSys:
	rm *.o primo *.mod .fuse*

cleanfuse:
	rm -rf ../outputdata/* .fuse*

cleandata:
	 rm -r ../outputdata/ #~/.local/share/Trash
	#sudo rm -fr ~/.local/share/Trash/*

cleanplot:
	rmdir plot

