include ${WM_PROJECT_DIR}/swAccPart/makefile.env

LIB_TARGET=${LIB_ROOT}/libOSspecific.a

OBJS = rowSubsections.o

default:${OBJS}
	swar cr ${LIB_TARGET} $^

rowSubsections.o:rowSubsections.cpp
	${CPP} ${CPP_FLAGS} -o $@ -c $(INCLUDE) $<
