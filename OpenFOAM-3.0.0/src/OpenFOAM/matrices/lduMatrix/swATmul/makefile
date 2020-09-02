include ${WM_PROJECT_DIR}/swAccPart/makefile.env

LIB_TARGET=${LIB_ROOT}/libOpenFOAM.a

OBJS = swAmul_slave.o swAmul.o swAmulTranslate.o swAmul_host.o util.o swAmulRedesigned.o

default:${OBJS}
	swar cur ${LIB_TARGET} $^

swAmul_slave.o:swAmul_slave.c
	${CC} $(SLAVE_FLAGS) -msimd -c -o $@ $(C_INCLUDE) $<
swAmul.o:swAmul.cpp
	${CPP} $(CPP_FLAGS) -o $@ -c $(INCLUDE) $<
swAmulTranslate.o:swAmulTranslate.cpp
	${CPP} $(CPP_FLAGS) -o $@ -c $(INCLUDE) $<
swAmul_host.o:swAmul_host.c
	${CC} $(HOST_FLAGS) -c -o $@ $(C_INCLUDE) $<
util.o:util.c
	${CC} $(HOST_FLAGS) -c -o $@ $(C_INCLUDE) $<
swScaling_factor.o:swScaling_factor.c
	${CC} $(SLAVE_FLAGS) -msimd -c -o $@ $(C_INCLUDE) $<
swAmulRedesigned.o:swAmulRedesigned.cpp
	${CPP} $(CPP_FLAGS) -o $@ -c $(INCLUDE) -I$(FV_ROOT)/lnInclude $<
