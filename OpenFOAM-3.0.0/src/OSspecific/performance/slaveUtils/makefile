include ${WM_PROJECT_DIR}/swAccPart/makefile.env

SLAVEFLAGS=-slave -msimd -OPT:IEEE_arith=1
LIB_TARGET=${LIB_ROOT}/libOSspecific.a

OBJS = slaveUtils.o

default:${OBJS}
	swar cru ${LIB_TARGET} $^

slaveUtils.o: slaveUtils.c
	${CC} ${SLAVEFLAGS} -g -c $^
