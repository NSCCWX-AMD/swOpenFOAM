#ifndef SWGAMGAgglomeration_H
#define SWGAMGAgglomeration_H

#include "GAMGAgglomeration.H"
#include "swRestInterStruct.h"
#include <vector>

namespace Foam
{
class GAMGAgglomeration;

class swRestInterMap
{
private:
	const GAMGAgglomeration& aggl_;
	// int32 copy of PtrList<labelField> restrictAddressing_;
	std::vector<std::vector<int> > restrictAddressing_int32;

public:
	swRestInterMap(const GAMGAgglomeration& aggl);
	
	friend class GAMGAgglomeration;

	static label bandSize_;
	static label minCellsUsingSW_;

	//- restrict
	static bool* restFirstUse_;
	static restStruct* restStructLevels_;

	//- face restrict
	static bool* faceRestFirstUse_;
	static restStruct* faceRestStructLevels_;

	//- interpolate
	static bool* interFirstUse_;
	static interStruct* interStructLevels_;

	//- agglomerate matrix upper
	static bool* aggMatrixUpperFirstUse_;
	static aggMatrixUpperStruct* aggMatrixUpperStructLevels_;

	void initRestInterMap();

	void initRestStruct
	(
		Field<scalar>& cf,
		const Field<scalar>& ff,
    	const label fineLevelIndex
	);

	void initInterStruct
	(
		Field<scalar>& ff,
	    const Field<scalar>& cf,
	    const label levelIndex
	);

	void agglomerateMatrixUpper
	(
		scalarField& coarseUpper,
		scalarField& coarseDiag,
		const scalarField& fineUpper,
		const label fineLevelIndex
	);

	void initAggMatrixUpperStruct
	(
		scalarField& coarseUpper,
		scalarField& coarseDiag,
		const scalarField& fineUpper,
		const label& fineLevelIndex
	);
};

} //- end namespace Foam

#endif
