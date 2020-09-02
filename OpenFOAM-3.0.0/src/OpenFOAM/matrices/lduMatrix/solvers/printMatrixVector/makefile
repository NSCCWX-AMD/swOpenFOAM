include ${WM_PROJECT_DIR}/swAccPart/makefile.env

LIB_TARGET=${LIB_ROOT}/libOpenFOAM.a

OBJS = printMatrixVector.o

default:${OBJS}
	swar cur ${LIB_TARGET} $^

printMatrixVector.o:printMatrixVector.cpp
	${CPP} $(CPP_FLAGS) -o $@ -c $(INCLUDE) $<
