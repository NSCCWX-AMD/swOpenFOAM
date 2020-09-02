include ${WM_PROJECT_DIR}/swAccPart/makefile.env

LIB_TARGET=${LIB_ROOT}/libOpenFOAM.a

OBJS = swScalarFieldOps.o

default:${OBJS}
	swar cur ${LIB_TARGET} $^

swScalarFieldOps.o:swScalarFieldOps.cpp
	${CPP} $(CPP_FLAGS) -o $@ -c $(INCLUDE) $<
