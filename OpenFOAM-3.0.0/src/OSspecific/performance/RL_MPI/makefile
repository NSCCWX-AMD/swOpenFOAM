include ${WM_PROJECT_DIR}/swAccPart/makefile.env

LIB_TARGET=${LIB_ROOT}/libOSspecific.a

OBJS = RlmpiInitializer.o rlmpi.o

default:${OBJS}
	swar cr ${LIB_TARGET} $^

RlmpiInitializer.o:RlmpiInitializer.cxx
	sw5CC ${HOST_FLAGS} -w -DSUNWAY -o $@ -c $<
#	sw5CC ${HOST_FLAGS} -w -DSUNWAY -o $@ -c $(INCLUDE) $<
rlmpi.o:rlmpi.c
	${CC} -slave -w -DSUNWAY -msimd -o $@ -c $(INCLUDE) $<
