objects = molecule.o force_field.o
comp = gfortran
OPT = -fbounds-check -O3 


molecule.o: molecule.f90
	$(comp) -c $(OPT) molecule.f90

force_field.o: force_field.f90
	$(comp) -c $(OPT) force_field.f90

.PHONY: clean_o all

all : molecule.o force_field.o

clean_o:
	rm -f $(objects)
