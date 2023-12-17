F90     = mpifort
F90OPT  = -Mpreprocess -fast -acc -gpu=cc80,managed -Minfo=acc
LOPT    = -Mpreprocess -fast -acc -gpu=cc80,managed -Minfo=acc
#F90OPT  = -Mpreprocess -fast -acc -gpu=managed,ccn90 -Minfo=acc  # for GH100
#LOPT    = -Mpreprocess -fast -acc -gpu=managed,ccn90 -Minfo=acc  # for GH100
TARGET  = ../../../run/ACC_U_s_c2

SRC     =         \
hpcmw_util.o hpcmw_fem_cntl.o hpcmw_fem_util.o hpcmw_solver_cntl.o\
hpcmw_fem_mesh.o hpcmw_solver_matrix.o hpcmw_all.o\
input_grid.o\
solver_SR_3.o solver_CG_3_SMP_novec.o \
solver33.o test1.o util.o sRCM.o sMC.o sCMRCM.o \
hpcmw_init.o input_cntl.o\
mat_con0.o mat_con1.o mat_ass_main.o mat_ass_bc.o hpcmw_finalize.o

OBJ     = ${SRC:.f=.o}

.f.o:
	$(F90) -c -loglist $(F90OPT) $<

$(TARGET):$(OBJ)
	$(F90) -o $@ $(LOPT) $(OBJ) $(LIB)

clean    :
	rm -f *.o *.mod *.log $(TARGET) *~
# Makefile End
