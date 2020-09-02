include ${WM_PROJECT_DIR}/swAccPart/makefile.env

LIB_TARGET=${LIB_ROOT}/libOSspecific.a

OBJS = Timer.o

default:${OBJS}
	swar cr ${LIB_TARGET} $^

Timer.o:Timer.cpp
	${CPP} ${CPP_FLAGS} -o $@ -c $(INCLUDE) $<