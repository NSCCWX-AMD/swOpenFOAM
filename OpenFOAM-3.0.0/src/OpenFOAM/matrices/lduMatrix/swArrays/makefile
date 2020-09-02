include ${WM_PROJECT_DIR}/swAccPart/makefile.env

LIB_TARGET=${LIB_ROOT}/libOpenFOAM.a

OBJS = vectorOps_host.o vectorOps_slave.o userFunc_SAXPY.o vectorCopy_host.o vectorCopy_slave.o gSum_host.o gSum_slave.o residualNormFactor_host.o residualNormFactor_slave.o

default:${OBJS}
	swar cur ${LIB_TARGET} $^

vectorOps_host.o:vectorOps_host.c
	${CC} $(HOST_FLAGS) -c -o $@ $(C_INCLUDE) $^
vectorOps_slave.o:vectorOps_slave.c
	${CC} $(SLAVE_FLAGS) -c -o $@ $(C_INCLUDE) $^
vectorCopy_host.o:vectorCopy_host.c
	${CC} $(HOST_FLAGS) -c -o $@ $(C_INCLUDE) $^
vectorCopy_slave.o:vectorCopy_slave.c
	${CC} $(SLAVE_FLAGS) -c -o $@ $(C_INCLUDE) $^
userFunc_SAXPY.o:userFunc_SAXPY.c
	${CC} $(SLAVE_FLAGS) -msimd -c -o $@ $(C_INCLUDE) $^
gSum_host.o:gSum_host.c
	${CC} $(HOST_FLAGS) -c -o $@ $(C_INCLUDE) $<
gSum_slave.o:gSum_slave.c
	${CC} $(SLAVE_FLAGS) -msimd -c -o $@ $(C_INCLUDE) $<
residualNormFactor_host.o:residualNormFactor_host.c
	${CC} $(HOST_FLAGS) -c -o $@ $(C_INCLUDE) $^
residualNormFactor_slave.o:residualNormFactor_slave.c
	${CC} $(SLAVE_FLAGS) -msimd -c -o $@ $(C_INCLUDE) $^
