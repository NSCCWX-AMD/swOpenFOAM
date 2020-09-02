include ${WM_PROJECT_DIR}/swAccPart/makefile.env

LIB_TARGET=${LIB_ROOT}/libOpenFOAM.a

OBJS = swPBiCGStab.o

default:${OBJS}
	swar cur ${LIB_TARGET} $^

swPBiCGStab.o:swPBiCGStab.C
	${CPP} $(CPP_FLAGS) -o $@ -c $(INCLUDE) $<
