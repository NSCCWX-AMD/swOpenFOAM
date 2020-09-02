include ${WM_PROJECT_DIR}/swAccPart/makefile.env

LIB_TARGET=${LIB_ROOT}/libOSspecific.a

OBJS = BlockOrdering.o BOrderUtils.o

default:${OBJS}
	swar cr ${LIB_TARGET} $^

BlockOrdering.o:BlockOrdering.c
	${CC} ${HOST_FLAGS} -o $@ -c $(INCLUDE) $<
BOrderUtils.o:BOrderUtils.c
	${CC} ${HOST_FLAGS} -o $@ -c $(INCLUDE) $<
