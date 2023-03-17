# makefile for compiling the code
# Contributors:
#   - Paolo Nicolini (Czech Technical University in Prague), paolo.nicolini22@gmail.com

FC = gfortran

# flags for optimization
FLAGS = -I. -ffree-line-length-none -fopenmp -O3
# flags for debug
#FLAGS = -I. -ffree-line-length-none -fopenmp -Wall -Wextra -Wconversion -pedantic -O -fcheck=all -g -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid

EXE = pindol

OBJECTS = modules.o main.o read_input.o read_refconf.o read_initconf.o read_restart.o read_eigen.o read_phi.o \
	setup_supercell.o initialize_velocities.o generate_distortions.o cartesian_to_normal.o normal_to_cartesian.o \
	calc_double_sum_phi.o calc_kinetic_energy.o calc_normal_accelerations.o calc_potential_energy.o run_nd.o \
	initialize_dvode.o dvode.o f.o write_finalconf.o write_restart.o deallocate_all.o credits.o

all: $(OBJECTS)
	$(FC) -o $(EXE) $(FLAGS) $(OBJECTS)

clean:
	-rm -f *__genmod.f90 *.o *.mod $(EXE)

dvode.f:
	wget https://computing.llnl.gov/sites/default/files/dvode.f

dvode.o: dvode.f
	$(FC) $(FLAGS) -c $<

%.o: %.f90
	$(FC) $(FLAGS) -c $<
