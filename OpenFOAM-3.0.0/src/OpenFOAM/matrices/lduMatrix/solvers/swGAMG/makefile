include ${WM_PROJECT_DIR}/swAccPart/makefile.env

LIB_TARGET=${LIB_ROOT}/libOpenFOAM.a

OBJS = swGAMGSolver.o swVcycle.o swScalingFactor.o scalingFactor_host.o scalingFactor_slave.o

default:${OBJS}
	swar cur ${LIB_TARGET} $^

swGAMGSolver.o:swGAMGSolver.cpp
	${CPP} ${CPP_FLAGS} -o $@ -c $(INCLUDE) $<
swScalingFactor.o:swScalingFactor.cpp
	${CPP} ${CPP_FLAGS} -o $@ -c $(INCLUDE) $<
swVcycle.o:swVcycle.cpp
	${CPP} ${CPP_FLAGS} -o $@ -c $(INCLUDE) $<

scalingFactor_host.o:scalingFactor_host.c
	${CC} $(HOST_FLAGS) -c -o $@ $(C_INCLUDE) $<
scalingFactor_slave.o:scalingFactor_slave.c
	${CC} $(SLAVE_FLAGS) -msimd -c -o $@ $(C_INCLUDE) $<
