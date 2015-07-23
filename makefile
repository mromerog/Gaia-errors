# Variables

FC = gfortran

FLAGS = -O 


####################################################################


#computes UVW and proper motions for the model and with errors
main_Gaiaerrors: main_Gaiaerrors.o coordrou.o Gaia-errors.o gaussdev.o
	$(FC) $(FLAGS) -o $@ main_Gaiaerrors.o coordrou.o Gaia-errors.o gaussdev.o


clean:  
	rm -f *.o
