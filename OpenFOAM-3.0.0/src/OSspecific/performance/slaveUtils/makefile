include ${WM_PROJECT_DIR}/swAccPart/makefile.env

SLAVEFLAGS=-mslave -msimd -mieee -D_SW_COMPILER_VERSION
LIB_TARGET=${LIB_ROOT}/libOSspecific.a

OBJS = slaveUtils.o

default:${OBJS}
	swar cru ${LIB_TARGET} $^

slaveUtils.o: slaveUtils.c
	${CC} ${SLAVEFLAGS} -g -c $^
