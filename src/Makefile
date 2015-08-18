EXEA =  ../calc_nbar
#EXEZ =  ../creat-HR3-mock-RSDMocks_for_Test

OBJS = ap_tools.o ap_cosmo_funs.o mainA_calc_nbar_wfkp.o
F90C = mpif90 -heap-arrays
F90FLAGS = 

default: $(EXEA)

$(EXEA): $(OBJS)
	$(F90C) -o $(EXEA)  ap_tools.o ap_cosmo_funs.o mainA_calc_nbar_wfkp.o
#$(EXEZ): $(OBJS)
#	$(F90C) -o $(EXEZ)  ap_tools.o ap_cosmo_funs.o mainZ_RSDMocks_for_Test.o        

ap_cosmo_funs.o: ap_tools.o
mainA_calc_nbar_wfkp.o:  ap_cosmo_funs.o

%.o: %.f90
	$(F90C) -c -heap-arrays $*.f90 $(F90FLAGS)

clean :
	rm -rf *.o *.mod $(EXEA) $(EXEB) $(EXEC) $(EXED) 
