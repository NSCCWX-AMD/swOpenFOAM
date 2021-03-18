 CC  = sw5cc
 CPP = swg++453
 HOST_FLAGS  = -host -OPT:IEEE_arith=2 -DWM_DP
 SLAVE_FLAGS = -slave -OPT:IEEE_arith=2 -DWM_DP
 CPP_FLAGS = -DSWTIMER=1
 ROOT = ${HOME}/lhb/OpenFOAM_RIGHT/OpenFOAM-3.0.0
 LIB_ROOT = ${ROOT}/platforms/linux64swgccDPInt32Opt/lib
 FOAM_ROOT = $(ROOT)/src/OpenFOAM
 OS_ROOT = ${ROOT}/src/OSspecific/POSIX
 FV_ROOT = ${ROOT}/src/finiteVolume
 MESH_ROOT = ${ROOT}/src/fvMesh
 PERF_ROOT = ${ROOT}/src/OSspecific
 INCLUDE = -fPIC -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=64 -Wall -Wextra -Wno-old-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -O3 -DNoRepository -ftemplate-depth-100 -DSW_SLAVE -I/usr/sw-mpp/mpi2/include -I/usr/sw-mpp/swcc/swgcc453-binary/cross-tools/include -I${FOAM_ROOT}/lnInclude -I${OS_ROOT}/lnInclude -I${FV_ROOT}/lnInclude -I${PERF_ROOT}/lnInclude
# C_INCLUDE=-I$(CORE_ROOT)/lnInclude/

all:gaussLaplacianSchemes.o
	swar cur $(LIB_ROOT)/libfiniteVolume.a $^

gaussLaplacianSchemes.o:gaussLaplacianSchemes.C
	${CPP} ${CPP_FLAGS} -o $@ -c $(INCLUDE) gaussLaplacianSchemes.C
#swGgradfRoutines.o:swGgradfRoutines.c
#	${CC} $(HOST_FLAGS) -c -o $@ $(C_INCLUDE) $^
#swGgradfRoutines_slave.o:swGgradfRoutines_slave.c
#	${CC} $(SLAVE_FLAGS) -CG:pjump_all -msimd -c -o $@ $(C_INCLUDE) $^

clean:
	rm -f *.o
