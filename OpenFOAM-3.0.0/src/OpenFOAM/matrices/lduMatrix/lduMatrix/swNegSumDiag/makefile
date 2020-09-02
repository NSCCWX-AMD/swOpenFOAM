include ${WM_PROJECT_DIR}/swAccPart/makefile.env

LIB_TARGET=${LIB_ROOT}/libOpenFOAM.a

OBJS = swNegSumDiag.o swNegSumDiag_reg.o swNegSumDiag_host.o swNegSumDiag_slave.o lduMatrixOperations.o

default:${OBJS}
	swar cr ${LIB_TARGET} $^

swNegSumDiag_reg.o:swNegSumDiag_reg.cpp
#	${CC} ${HOST_FLAGS} -w -o $@ -c $(INCLUDE) $<
	sw5CC -host -O1 -o $@ -c $(C_INCLUDE) $<
swNegSumDiag_host.o:swNegSumDiag_host.c
	${CC} ${HOST_FLAGS} -o $@ -c $(INCLUDE) $<
swNegSumDiag.o:swNegSumDiag.c
	${CC} ${HOST_FLAGS} -o $@ -c $(INCLUDE) $<
swNegSumDiag_slave.o:swNegSumDiag_slave.c
	${CC} ${SLAVE_FLAGS} -msimd -O1 -o $@ -c $(C_INCLUDE) $<
lduMatrixOperations.o:../lduMatrixOperations.C
	${CPP} ${CPP_FLAGS} -o $@ -c $(INCLUDE) $<

rollback:
	swar cr ${LIB_ROOT}/libOpenFOAM.a ${LIB_ROOT}/../src/OpenFOAM/matrices/lduMatrix/lduMatrix/lduMatrixOperations.o
