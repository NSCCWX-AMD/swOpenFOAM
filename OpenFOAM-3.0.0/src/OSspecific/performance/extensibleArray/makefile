include ${WM_PROJECT_DIR}/swAccPart/makefile.env

LIB_TARGET=${LIB_ROOT}/libOSspecific.a

OBJS = extensibleLabelArray.o extensibleScalarArray.o

default:${OBJS}
	swar cr ${LIB_TARGET} $^

extensibleLabelArray.o:extensibleLabelArray.c
	${CC} ${HOST_FLAGS} -o $@ -c $(INCLUDE) $<
extensibleScalarArray.o:extensibleScalarArray.c
	${CC} ${HOST_FLAGS} -o $@ -c $(INCLUDE) $<
