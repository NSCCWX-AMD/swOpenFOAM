include ${WM_PROJECT_DIR}/swAccPart/makefile.env

LIB_TARGET=${LIB_ROOT}/libOpenFOAM.a

OBJS = srPCG.o

default:${OBJS}
	swar cur ${LIB_TARGET} $^

srPCG.o:srPCG.cpp
	${CPP} $(CPP_FLAGS) -o $@ -c $(INCLUDE) $<
