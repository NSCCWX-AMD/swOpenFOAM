#############################################################
# path setting
#############################################################
PROJECT=${PWD}
LIBPATH=${PROJECT}/../lib
INCLUDE=${PROJECT}/include
OBJPATH=${PROJECT}/objects

export PROJECT
export LIBPATH
export INCLUDE
export OBJPATH
export LIBROOT

#############################################################
# build tool set setting
#############################################################
CC="sw5cc -host -OPT:IEEE_arithmetic=2 -fPIC"
#CXX="mpiCC -OPT:IEEE_arithmetic=2"
CXX="mpiswg++ -mieee -fPIC"
SLAVECC="sw5cc -slave -msimd -OPT:IEEE_arithmetic=2 -fPIC"
AR=swar cru
RANLIB=swranlib
LD=${CXX}

#############################################################
# source directory
#############################################################
SRCPATH:=${PROJECT}/tools \
		 ${PROJECT}/topology \
		 ${PROJECT}/iterator \
		 ${PROJECT}/SW_Interface \
		 ${PROJECT}/wrappedInterface 

TESTPATH=${PROJECT}/test

#############################################################
# make
#############################################################
.PHONY:all libobjs makepath

${LIBPATH}/libUNAT.a: libUNAT.a
	@cp $< $@
	@echo -e "\033[32mBuild completed! \033[0m"

libUNAT.a: libobjs
	${AR} $@ ${OBJPATH}/*.o
	${RANLIB} $@ 
	
test: libUNAT.a
	(cd ${path}; make CC=${CC} CXX=${CXX} OBJPATH=${OBJPATH} \
	INLCUDE=${INLCUDE} SLAVECC=${SLAVECC} RANLIB=${RANLIB} LD=${LD})

libobjs: makepath
	@for path in $(SRCPATH); do \
		echo -e "\033[32mEntering path $$path \033[0m"; \
		(cd $$path; \
		make -f ../Makefile.rule CC=${CC} CXX=${CXX} OBJPATH=${OBJPATH} \
		INLCUDE=${INLCUDE} SLAVECC=${SLAVECC} RANLIB=${RANLIB}); \
	done

makepath:
	@if [ ! -d ${OBJPATH} ]; then \
		echo -e "\033[32mMaking directory ${OBJPATH} \033[0m"; \
		mkdir ${OBJPATH}; \
	fi
	@if [ ! -d ${INCLUDE} ]; then \
		echo -e "\033[32mMaking directory ${INCLUDE} \033[0m"; \
		mkdir ${INCLUDE}; \
	fi
	@if [ ! -d ${LIBPATH} ]; then \
		echo -e "\033[32mMaking directory ${LIBPATH} \033[0m"; \
		mkdir ${LIBPATH}; \
	fi

.PHONY:clean
clean:
	rm -rf ${LIBPATH}/libUNAT.a ${OBJPATH}/* ${INCLUDE}/*
	@for path in $(SRCPATH); do \
		(cd $$path; \
		make -f ../Makefile.rule clean); \
	done
