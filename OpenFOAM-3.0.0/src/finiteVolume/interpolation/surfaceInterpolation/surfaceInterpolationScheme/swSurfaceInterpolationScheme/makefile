include ${WM_PROJECT_DIR}/swAccPart/makefile.env

LIB_TARGET=${LIB_ROOT}/libfiniteVolume.a

OBJS = interpolate1_host.o interpolate1_slave.o swSurfaceInterpolationScheme_C.o swSurfaceInterpolationScheme.o

default:${OBJS}
	swar cur ${LIB_TARGET} $^

swSurfaceInterpolationScheme.o:swSurfaceInterpolationScheme.cpp
	${CPP} $(CPP_FLAGS) -o $@ -c $(INCLUDE) $<
swSurfaceInterpolationScheme_C.o:swSurfaceInterpolationScheme_C.cpp
	${CPP} $(CPP_FLAGS) -o $@ -c $(C_INCLUDE) $<
interpolate1_host.o:interpolate1_host.c
	${CC} $(HOST_FLAGS) -c -o $@ $(C_INCLUDE) $<
interpolate1_slave.o:interpolate1_slave.c
	${CC} $(SLAVE_FLAGS) -msimd -c -o $@ $(C_INCLUDE) $<
