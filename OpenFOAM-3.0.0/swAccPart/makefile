include ./makefile.env

default:
	make all |&tee log.make

all:
	@cd performance/slaveUtils && make default
	@cd performance/CG_switch &&  make default
	@cd performance/timers &&  make default
#	////@cd performance/rowSubsections && make default
#	////@cd performance/BlockOrdering && make default
#	////@cd performance/extensibleArray && make default
#	////@cd swNegSumDiag &&  make default
#   ////@cd lduMatrix &&  make default
	@cd swArrays &&  make default
	@cd swATmul && make default
	#@cd swGAMG && make default
	@cd swPBiCGSTAB && make default
	#@cd swCHEBY && make default
	#@cd swDiagonalPreconditioner && make default
	#@cd swDILUPreconditioner && make default
	#@cd srPCG && make default
	#@cd swScalarFieldOps && make default
	#@cd swField && make default
#   ////@cd swGaussGrad && make default
##	////@cd swSurfaceIntegrate &&  make default
#   ////@cd swSurfaceInterpolationScheme &&  make default
	#@cd swGAMGAgglomeration && make default
#   ////@cd swGAMGPreconditioner && make default
	#@cd printMatrixVector && make default

cleanAll:
	@cd performance/slaveUtils && make cleanObjs
	@cd performance/CG_switch &&  make cleanObjs
	@cd performance/timers &&  make cleanObjs
	@cd performance/rowSubsections && make cleanObjs
	@cd performance/BlockOrdering && make cleanObjs
	@cd performance/extensibleArray && make cleanObjs
	@cd swNegSumDiag &&  make cleanObjs
	@cd swArrays &&  make cleanObjs
	@cd swATmul && make cleanObjs
	@cd swGAMG && make cleanObjs
	@cd swPBiCGSTAB && make cleanObjs
	@cd swCHEBY && make cleanObjs
	@cd swDiagonalPreconditioner && make cleanObjs
	@cd swDILUPreconditioner && make cleanObjs
	@cd srPCG && make cleanObjs
	@cd swScalarFieldOps && make cleanObjs
	@cd swField && make cleanObjs
	@cd swGaussGrad && make cleanObjs
	@cd swSurfaceIntegrate &&  make cleanObjs
	@cd swSurfaceInterpolationScheme &&  make cleanObjs
	@cd swGAMGAgglomeration && make cleanObjs
	@cd swGAMGPreconditioner && make cleanObjs
	@cd printMatrixVector && make cleanObjs
