include ${WM_PROJECT_DIR}/swAccPart/makefile.env

LIB_TARGET=${LIB_ROOT}/libfiniteVolume.a

OBJS = swFvcSurfaceIntegrate.o surfaceIntegrate_host.o surfaceIntegrate_link.o surfaceIntegrate_slave.o

default:${OBJS}
	swar cur ${LIB_TARGET} $^

swFvcSurfaceIntegrate.o:swFvcSurfaceIntegrate.cpp
	${CPP} ${CPP_FLAGS} -o $@ -c $(INCLUDE) $<
surfaceIntegrate_host.o:surfaceIntegrate_host.cpp
	${CPP} ${CPP_FLAGS} -o $@ -c $(C_INCLUDE) $<
surfaceIntegrate_link.o:surfaceIntegrate_link.c
	${CC} $(HOST_FLAGS) -c -o $@ $(C_INCLUDE) $<
surfaceIntegrate_slave.o:surfaceIntegrate_slave.c
	${CC} $(SLAVE_FLAGS) -msimd -c -o $@ $(C_INCLUDE) $<
