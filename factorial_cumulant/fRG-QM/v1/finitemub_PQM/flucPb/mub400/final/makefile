# Start of the makefile
# Defining variables
objects=chebft.o QM.o 
# Makefile
exe: $(objects)
	gfortran -o exe $(objects) 
%.o: %.f90
	gfortran -c -O2 -g $< 
# Cleaning everything
clean:
	rm $(objects)
# End of the makefile
