include ${WM_PROJECT_DIR}/swAccPart/makefile.env

LIB_TARGET=${LIB_ROOT}/libOpenFOAM.a

OBJS = swChebyIfpack2srPCG.o

default:${OBJS}
	swar cur ${LIB_TARGET} $^

swChebyIfpack2srPCG.o:swChebyIfpack2srPCG.cpp
	${CPP} $(CPP_FLAGS) -o  $@  -c $(INCLUDE) $<
