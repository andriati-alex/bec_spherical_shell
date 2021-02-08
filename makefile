
# ***************************** MAKEFILE ***************************** #

# $ make main_program



obj = memoryHandling.o 	   \
	  arrayOperations.o    \
	  tridiagonalSystems.o \
	  calculus.o           \
	  observables.o        \
	  inout.o              \
	  imagtimeIntegrator.o



   # ------------------------------------------------------------------ #

                         ###     EXECUTABLES     ###

   # ------------------------------------------------------------------ #



time_evolution : libgp.a exe/time_evolution.c
	icc -o ./bin/time_evolution exe/time_evolution.c -lm -mkl -qopenmp \
		-L./lib -I./include -lgp -O3




# Libraries to be linked
# ----------------------

libgp.a : $(obj)
	ar rcs libgp.a $(obj)
	mv libgp.a lib
	mv $(obj) build





# Object files to the library
# ---------------------------

inout.o : src/inout.c
	icc -c -O3 -I./include src/inout.c

memoryHandling.o : src/memoryHandling.c
	icc -c -O3 -I./include src/memoryHandling.c

arrayOperations.o : src/arrayOperations.c
	icc -c -O3 -qopenmp -I./include src/arrayOperations.c

tridiagonalSystems.o : src/tridiagonalSystems.c
	icc -c -O3 -qopenmp -I./include src/tridiagonalSystems.c

calculus.o : src/calculus.c
	icc -c -O3 -qopenmp -I./include src/calculus.c

observables.o : src/observables.c
	icc -c -O3 -qopenmp -I./include src/observables.c

newtoncg.o : src/newtoncg.c
	icc -c -O3 -qopenmp -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -I./include src/newtoncg.c

imagtimeIntegrator.o : src/imagtimeIntegrator.c
	icc -c -O3 -qopenmp -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -I./include src/imagtimeIntegrator.c

clean :
	-rm build/*.o
	-rm lib/lib*
	-rm bin/*
