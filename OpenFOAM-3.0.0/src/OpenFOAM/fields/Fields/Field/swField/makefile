include ${WM_PROJECT_DIR}/swAccPart/makefile.env

LIB_TARGET=${LIB_ROOT}/libOpenFOAM.a

OBJS = swField.o swField_host.o swField_slave.o

default:${OBJS}
	swar cur ${LIB_TARGET} $^

swField.o: swField.cpp
	${CPP} $(CPP_FLAGS) -o $@ -c $(INCLUDE) $<

swField_host.o: swField_host.c
	${CC} $(HOST_FLAGS) -c -o $@ $(C_INCLUDE) $<

swField_slave.o: swField_slave.c
	${CC} $(SLAVE_FLAGS) -msimd -c -o $@ $(C_INCLUDE) $<
