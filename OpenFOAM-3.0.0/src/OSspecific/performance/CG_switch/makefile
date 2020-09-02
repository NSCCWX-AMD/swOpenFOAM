include ${WM_PROJECT_DIR}/swAccPart/makefile.env

LIB_TARGET=${LIB_ROOT}/libOSspecific.a

OBJS = CG_switch.o athread_switch.o

default: ${OBJS}
	swar cr ${LIB_TARGET} $^

CG_switch.o:CG_switch.cpp
	${CPP} -o $@ -c $(INCLUDE) $<

athread_switch.o:athread_switch.c
	${CC} $(HOST_FLAGS) -c -o $@ $^
