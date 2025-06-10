FC = gfortran

FLAGS = -I. -ffree-line-length-none 

# debug
#FLAGS+= -Wall -Wextra -Wconversion -pedantic -O -fcheck=all -g -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid 

# opt
FLAGS+= -O3

EXE = q4phind

OBJECTS = main.o init.o check_qpt.o write_qordered.o write_qopt.o check_qopt.o\
          show_logo.o credits.o

all: $(OBJECTS)
	$(FC) -o $(EXE) $(FLAGS) $(OBJECTS)

clean:
	-rm -f *__genmod.f90 *.o *.mod $(EXE)

%.o: %.f90
	$(FC) $(FLAGS) -c $<
