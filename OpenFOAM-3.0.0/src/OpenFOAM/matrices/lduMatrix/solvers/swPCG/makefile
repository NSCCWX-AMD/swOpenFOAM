include ${WM_PROJECT_DIR}/swAccPart/makefile.env

LIB_TARGET=${LIB_ROOT}/libOpenFOAM.a

OBJS = swPCG.o

default:${OBJS}
	swar cur ${LIB_TARGET} $^

swPCG.o:swPCG.C
	${CPP} $(CPP_FLAGS) -o $@ -c $(INCLUDE) $<