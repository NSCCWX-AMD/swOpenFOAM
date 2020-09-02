include ${WM_PROJECT_DIR}/swAccPart/makefile.env

LIB_TARGET=${LIB_ROOT}/libOpenFOAM.a

OBJS = swDiagonalPreconditioner.o

default:${OBJS}
	swar cur ${LIB_TARGET} $^

swDiagonalPreconditioner.o:swDiagonalPreconditioner.cpp
	${CPP} $(CPP_FLAGS) -o  $@  -c $(INCLUDE) $<
