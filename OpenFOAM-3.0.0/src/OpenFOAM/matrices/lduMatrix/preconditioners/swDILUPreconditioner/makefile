include ${WM_PROJECT_DIR}/swAccPart/makefile.env

LIB_TARGET=${LIB_ROOT}/libOpenFOAM.a

OBJS = swDILUPreconditioner.o

default:${OBJS}
	swar cur ${LIB_TARGET} $^

swDILUPreconditioner.o:swDILUPreconditioner.cpp
	${CPP} $(CPP_FLAGS) -o  $@  -c $(INCLUDE) $<
