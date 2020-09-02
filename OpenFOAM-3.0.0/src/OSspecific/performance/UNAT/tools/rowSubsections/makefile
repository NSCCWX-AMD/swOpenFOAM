all:rowSubsections.o
	swar cr $(LIB_ROOT)/libOSspecific.a $^

rowSubsections.o:rowSubsections.cpp
	${CPP} ${CPP_FLAGS} -o $@ -c $(INCLUDE) $<
clean:
	rm *.o
	swar d $(LIB_ROOT)/libOSspecific.a rowSubsections.o

