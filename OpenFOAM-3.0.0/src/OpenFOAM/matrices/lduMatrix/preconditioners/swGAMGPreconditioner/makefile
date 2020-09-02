include ${WM_PROJECT_DIR}/swAccPart/makefile.env

LIB_TARGET=${LIB_ROOT}/libOpenFOAM.a

OBJS = swGAMGPreconditioner.o

default:${OBJS}
	swar cur ${LIB_TARGET} $^

swGAMGPreconditioner.o:swGAMGPreconditioner.C
	${CPP} $(CPP_FLAGS) -o  $@  -c $(INCLUDE) $<
