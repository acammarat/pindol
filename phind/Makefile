FC = gfortran

FLAGS = -I. -ffree-line-length-none -fopenmp

# debug
#FLAGS+= -Wall -Wextra -Wconversion -pedantic -O -fcheck=all -g -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid 

# opt
FLAGS+= -O3

EXE = phind

OBJECTS = main.o \
          init.o show_logo.o read_eigen.o \
          calc_delta.o calc_phiun.o reorder_phi.o write_map.o deallocate_all.o \
          credits.o

all: $(OBJECTS)
	$(FC) -o $(EXE) $(FLAGS) $(OBJECTS)

clean:
	-rm -f *__genmod.f90 *.o *.mod $(EXE)

%.o: %.f90
	$(FC) $(FLAGS) -c $<
