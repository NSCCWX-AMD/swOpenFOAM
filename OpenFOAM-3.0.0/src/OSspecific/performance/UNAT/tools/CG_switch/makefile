defult: CG_switch.o athread_switch.o
	swar cr $(LIB_ROOT)/libOSspecific.a $^
#cp ./*.h ../../lnInclude/

CG_switch.o:CG_switch.cpp
	${CPP} -o $@ -c $(INCLUDE) $<

athread_switch.o:athread_switch.c
	${CC} $(HOST_FLAGS) -c -o $@ $^

clean:
	rm *.o
	swar d $(LIB_ROOT)/libOSspecific.a CG_switch.o athread_switch.o

