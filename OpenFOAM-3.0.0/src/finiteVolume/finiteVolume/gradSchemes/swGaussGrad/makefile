include ${WM_PROJECT_DIR}/swAccPart/makefile.env

LIB_TARGET=${LIB_ROOT}/libfiniteVolume.a

OBJS = swGaussGrad.o swGaussGradSub.o swGgradfRoutines.o swGgradfRoutines_slave.o

default:${OBJS}
	swar cur ${LIB_TARGET} $^

swGaussGrad.o:swGaussGrad.cpp
	${CPP} ${CPP_FLAGS} -o $@ -c $(INCLUDE) $<
swGaussGradSub.o:swGaussGradSub.cpp
	${CPP} ${CPP_FLAGS} -o $@ -c $(INCLUDE) $<
swGgradfRoutines.o:swGgradfRoutines.c
	${CC} $(HOST_FLAGS) -c -o $@ $(C_INCLUDE) $^
swGgradfRoutines_slave.o:swGgradfRoutines_slave.c
	${CC} $(SLAVE_FLAGS) -CG:pjump_all -msimd -c -o $@ $(C_INCLUDE) $^
