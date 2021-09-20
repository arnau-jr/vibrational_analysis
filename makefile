objects = molecule.o force_field.o vibration.o
comp = gfortran
OPT = -fbounds-check -O3 


molecule.o: molecule.f90
	$(comp) -c $(OPT) molecule.f90

force_field.o: force_field.f90
	$(comp) -c $(OPT) force_field.f90

vibration.o: vibration.f90 molecule.o force_field.o
	$(comp) -c $(OPT) vibration.f90 molecule.f90 force_field.f90 -llapack

.PHONY: clean_o all

all : molecule.o force_field.o vibration.o
	touch molecule.o
	touch force_field.o
	touch vibration.o

clean_o:
	rm -f $(objects)
