include ${WM_PROJECT_DIR}/swAccPart/makefile.env

LIB_TARGET=${LIB_ROOT}/libOpenFOAM.a

OBJS = restrictData_slave.o restrictData_host.o interpolateData_slave.o interpolateData_host.o

default:${OBJS}
	swar cur ${LIB_TARGET} $^

restrictData_host.o:restrictData_host.c
	${CC} $(HOST_FLAGS) -c -o $@ $(C_INCLUDE) $<
restrictData_slave.o:restrictData_slave.c
	${CC} $(SLAVE_FLAGS) -msimd -c -o $@ $(C_INCLUDE) $<

interpolateData_host.o:interpolateData_host.c
	${CC} $(HOST_FLAGS) -c -o $@ $(C_INCLUDE) $<
interpolateData_slave.o:interpolateData_slave.c
	${CC} $(SLAVE_FLAGS) -msimd -c -o $@ $(C_INCLUDE) $<
