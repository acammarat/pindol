FC = gfortran

FLAGS = -I. -ffree-line-length-none 

# debug
#FLAGS+= -Wall -Wextra -Wconversion -pedantic -O -fcheck=all -g -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid 

# opt
FLAGS+= -O3

EXE = nd2xyz

OBJECTS = main.o \
	init.o read_eigen.o normal_to_cartesian.o show_logo.o \
        write_xyz.o write_xdot.o write_xddot.o write_last.o \
        deallocate_all.o credits.o

all: $(OBJECTS)
	$(FC) -o $(EXE) $(FLAGS) $(OBJECTS)

clean:
	-rm -f *__genmod.f90 *.o *.mod $(EXE)

%.o: %.f90
	$(FC) $(FLAGS) -c $<

