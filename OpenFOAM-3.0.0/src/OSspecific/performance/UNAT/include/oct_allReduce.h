#ifndef OCT_ALLREDUCE_H
#define OCT_ALLREDUCE_H

//ARRAY_PUTR(sData, size, destC); \
// Row All Reduce
#define oct_rowAllReduce( buffer, type, size, opt ) \
{ \
	if(  sizeof(oct_pack) < sizeof(type) ||\
				sizeof(oct_pack) % sizeof(type) != 0 )\
	{\
		if(_MYID == 0)\
			LOG("size of oct_pack is not integeral times of the size of type!\n");\
		exit(-1);\
	}\
\
	swInt rowID = ROW(_MYID);\
	swInt colID = COL(_MYID);\
\
	swInt i;\
	for(i = 0; i < 3; i++)\
	{\
		if( ( colID >> (2-i) ) % 2 == 0 )\
		{\
			swInt destC = colID + (swInt) ( 8 >> (i+1) ); \
			type *sData = (type*)buffer; \
			ARRAY_PUTR(sData, size, destC); \
		} \
		else\
		{\
			swInt packNum = (size) / PACKSIZE;\
			swInt restCount = \
				( (size) % PACKSIZE );\
			type* realData = (type*) buffer;\
			swInt j; \
			swInt typeCount = 0;\
			oct_pack pack;\
			for(j = 0; j < packNum; j++)\
			{\
				REG_GETR( *(int256*) &pack );\
				swInt k;\
				for(k = 0; k < pack.valSize; k++)\
					opt(realData[typeCount+k], pack.val[k]);\
				typeCount += pack.valSize;\
			}\
			if(restCount != 0)\
			{\
				REG_GETR( *(int256*) &pack );\
				swInt k;\
				for(k = 0; k < pack.valSize; k++)\
					opt(realData[typeCount+k], pack.val[k]);\
				typeCount += pack.valSize;\
			}\
		}\
		\
		ROWSYN;\
\
	}\
	for(i = 0; i < 3; i++)\
	{\
		if( ( colID >> (2-i) ) % 2 != 0 )\
		{\
			swInt destC = colID - (swInt) ( 8 >> (i+1) ); \
			type *sData = (type*)buffer; \
			ARRAY_PUTR(sData, size, destC); \
		}\
		else\
		{\
			swInt packNum = (size) / PACKSIZE;\
			swInt restCount = \
				( (size) % PACKSIZE );\
			type* rData = (type*) buffer;\
			swInt j; \
			swInt typeCount = 0;\
			oct_pack pack;\
			for(j = 0; j < packNum; j++)\
			{\
				REG_GETR( *(int256*) &pack );\
				UNPACK_DATA(pack, rData); \
			}\
			if(restCount != 0)\
			{\
				REG_GETR( *(int256*) &pack );\
				UNPACK_DATA(pack, rData); \
			}\
		}\
		\
		ROWSYN;\
\
	}\
\
	ALLSYN;\
\
}

// Column All Reduce
#define oct_colAllReduce( buffer, type, size, opt ) \
{ \
	if(  sizeof(oct_pack) < sizeof(type) ||\
				sizeof(oct_pack) % sizeof(type) != 0 )\
	{\
		if(_MYID == 0)\
			LOG("size of oct_pack is not integeral times of the size of type!\n");\
		exit(-1);\
	}\
\
	swInt rowID = ROW(_MYID);\
	swInt colID = COL(_MYID);\
\
	swInt i;\
	for(i = 0; i < 3; i++)\
	{\
		if( ( rowID >> (2-i) ) % 2 == 0 )\
		{\
			swInt destC = rowID + (swInt) ( 8 >> (i+1) ); \
			destC = destC * 8 + colID; \
			type *sData = (type*)buffer; \
			ARRAY_PUTC(sData, size, destC); \
		}\
		else\
		{\
			swInt packNum = (size) / PACKSIZE;\
			swInt restCount = \
				( (size) % PACKSIZE );\
			type* realData = (type*) buffer;\
			swInt j; \
			swInt typeCount = 0;\
			oct_pack pack;\
			for(j = 0; j < packNum; j++)\
			{\
				REG_GETC( *(int256*) &pack );\
				swInt k;\
				for(k = 0; k < pack.valSize; k++)\
					opt(realData[typeCount+k], pack.val[k]);\
				typeCount += pack.valSize;\
			}\
			if(restCount != 0)\
			{\
				REG_GETC( *(int256*) &pack );\
				swInt k;\
				for(k = 0; k < pack.valSize; k++)\
					opt(realData[typeCount+k], pack.val[k]);\
				typeCount += pack.valSize;\
			}\
		}\
		\
		COLSYN;\
\
	}\
	for(i = 0; i < 3; i++)\
	{\
		if( ( rowID >> (2-i) ) % 2 != 0 )\
		{\
			swInt destC = rowID - (swInt) ( 8 >> (i+1) ); \
			destC = destC * 8 + colID; \
			type *sData = (type*)buffer; \
			ARRAY_PUTC(sData, size, destC); \
		}\
		else\
		{\
			swInt packNum = (size) / PACKSIZE;\
			swInt restCount = \
				( (size) % PACKSIZE );\
			type* rData = (type*) buffer;\
			swInt j; \
			swInt typeCount = 0;\
			oct_pack pack;\
			for(j = 0; j < packNum; j++)\
			{\
				REG_GETC( *(int256*) &pack );\
				UNPACK_DATA(pack, rData); \
			}\
			if(restCount != 0)\
			{\
				REG_GETC( *(int256*) &pack );\
				UNPACK_DATA(pack, rData); \
			}\
		}\
		\
		COLSYN;\
\
	}\
\
	ALLSYN;\
\
}

//// All Reduce
//#define oct_allReduce( buffer, type, size, opt ) \
//{ \
//	if(  sizeof(oct_pack) < sizeof(type) ||\
//				sizeof(oct_pack) % sizeof(type) != 0 )\
//	{\
//		if(_MYID == 0)\
//			LOG("size of oct_pack is not integeral times of the size of type!\n");\
//		exit(-1);\
//	}\
//\
//	swInt rowID = (swInt) (_MYID >> 3);\
//	swInt colID = (swInt) (_MYID%8);\
//\
//	swInt i;\
//	for(i = 0; i < 3; i++)\
//	{\
//		if( ( rowID >> (2-i) ) % 2 == 0 )\
//		{\
//			swInt destR = rowID + (swInt) ( 8 >> (i+1) ); \
//			oct_pack* packs = (oct_pack*) (buffer);\
//			swInt packSize = (size)*sizeof(type) / sizeof(oct_pack);\
//			swInt restCount = \
//				( ( (size)*sizeof(type) ) % sizeof(oct_pack) ) \
//				/ sizeof(type);\
//			swInt j;\
//			for(j = 0; j < packSize; j++)\
//				REG_PUTC(packs[j], destR);\
//			if( restCount != 0 )\
//			{\
//				oct_pack restPack;\
//				for(j = 0; j < restCount; j++)\
//					((type*) &restPack)[j] = \
//						((type*) &packs[packSize])[j];\
//				REG_PUTC(restPack, destR);\
//			}	\
//		}\
//		else\
//		{\
//			swInt srcR = rowID - (swInt) (8 >> (i+1)); \
//			swInt packSize = (size)*sizeof(type) / sizeof(oct_pack);\
//			swInt restCount = \
//				( ( (size)*sizeof(type) ) % sizeof(oct_pack) ) \
//				/ sizeof(type);\
//			type* realData = (type*) buffer;\
//			swInt numTypeInPack = sizeof(oct_pack) / sizeof(type);\
//			swInt j, typeCount = 0;\
//			oct_pack pack;\
//			for(j = 0; j < packSize; j++)\
//			{\
//				REG_GETC( pack );\
//				swInt k;\
//				for(k = 0; k < numTypeInPack; k++)\
//					opt(realData[typeCount+k], ((type*) &pack)[k]);\
//				typeCount += numTypeInPack;\
//			}\
//			if(restCount != 0)\
//			{\
//				REG_GETC( pack );\
//				swInt k;\
//				for(k = 0; k < restCount; k++)\
//					opt(realData[typeCount+k], ((type*) &pack)[k]);\
//				typeCount += restCount;\
//			}\
//			if(typeCount != size)\
//			{\
//				LOG("received wrong number of data!\n");\
//				printf("typeCount is %d , size is %d\n", typeCount, size);\
//				exit(-1);\
//			}\
//		}\
//		\
//		COLSYN;\
//\
//		if( ( rowID >> (2-i) ) % 2 != 0 )\
//		{\
//			swInt destR = rowID - (swInt) ( 8 >> (i+1) ); \
//			oct_pack* packs = (oct_pack*) (buffer);\
//			swInt packSize = (size)*sizeof(type) / sizeof(oct_pack);\
//			swInt restCount = \
//				( ( (size)*sizeof(type) ) % sizeof(oct_pack) ) \
//				/ sizeof(type);\
//			swInt j;\
//			for(j = 0; j < packSize; j++)\
//				REG_PUTC(packs[j], destR);\
//			if( restCount != 0 )\
//			{\
//				oct_pack restPack;\
//				for(j = 0; j < restCount; j++)\
//					((type*) &restPack)[j] = \
//						((type*) &packs[packSize])[j];\
//				REG_PUTC(restPack, destR);\
//			}	\
//		}\
//		else\
//		{\
//			swInt srcR = rowID + (swInt) (8 >> (i+1)); \
//			swInt packSize = (size)*sizeof(type) / sizeof(oct_pack);\
//			swInt restCount = \
//				( ( (size)*sizeof(type) ) % sizeof(oct_pack) ) \
//				/ sizeof(type);\
//			type* realData = (type*) buffer;\
//			swInt numTypeInPack = sizeof(oct_pack) / sizeof(type);\
//			swInt j, typeCount = 0;\
//			oct_pack pack;\
//			for(j = 0; j < packSize; j++)\
//			{\
//				REG_GETC( pack );\
//				swInt k;\
//				for(k = 0; k < numTypeInPack; k++)\
//					realData[typeCount+k] = ((type*) &pack)[k];\
//				typeCount += numTypeInPack;\
//			}\
//			if(restCount != 0)\
//			{\
//				REG_GETC( pack );\
//				swInt k;\
//				for(k = 0; k < restCount; k++)\
//					realData[typeCount+k] = ((type*) &pack)[k];\
//				typeCount += restCount;\
//			}\
//			if(typeCount != size)\
//			{\
//				LOG("received wrong number of data!\n");\
//				exit(-1);\
//			}\
//		}\
//		COLSYN;\
//	}\
//\
//	ALLSYN;\
//\
//	for(i = 0; i < 3; i++)\
//	{\
//		if( ( colID >> (2-i) ) % 2 == 0 )\
//		{\
//			swInt destC = colID + (swInt) ( 8 >> (i+1) ); \
//			oct_pack* packs = (oct_pack*) (buffer);\
//			swInt packSize = (size)*sizeof(type) / sizeof(oct_pack);\
//			swInt restCount = \
//				( ( (size)*sizeof(type) ) % sizeof(oct_pack) ) \
//				/ sizeof(type);\
//			swInt j;\
//			for(j = 0; j < packSize; j++)\
//				REG_PUTR(packs[j], destC);\
//			if( restCount != 0 )\
//			{\
//				oct_pack restPack;\
//				for(j = 0; j < restCount; j++)\
//				((type*) &restPack)[j] = \
//					((type*) &packs[packSize])[j];\
//				REG_PUTR(restPack, destC);\
//			}	\
//		}\
//		else\
//		{\
//			swInt srcC = colID - (swInt) (8 >> (i+1)); \
//			swInt packSize = (size)*sizeof(type) / sizeof(oct_pack);\
//			swInt restCount = \
//				( ( (size)*sizeof(type) ) % sizeof(oct_pack) ) \
//				/ sizeof(type);\
//			type* realData = (type*) buffer;\
//			swInt numTypeInPack = sizeof(oct_pack) / sizeof(type);\
//			swInt j, typeCount = 0;\
//			oct_pack pack;\
//			for(j = 0; j < packSize; j++)\
//			{\
//				REG_GETR( pack );\
//				swInt k;\
//				for(k = 0; k < numTypeInPack; k++)\
//					opt(realData[typeCount+k], ((type*) &pack)[k]);\
//				typeCount += numTypeInPack;\
//			}\
//			if(restCount != 0)\
//			{\
//				REG_GETR( pack );\
//				swInt k;\
//				for(k = 0; k < restCount; k++)\
//					opt(realData[typeCount+k], ((type*) &pack)[k]);\
//				typeCount += restCount;\
//			}\
//			if(typeCount != size)\
//			{\
//				LOG("received wrong number of data!\n");\
//				exit(-1);\
//			}\
//		}\
//		\
//		ROWSYN;\
//\
//		if( ( colID >> (2-i) ) % 2 != 0 )\
//		{\
//			swInt destC = colID - (swInt) ( 8 >> (i+1) ); \
//			oct_pack* packs = (oct_pack*) (buffer);\
//			swInt packSize = (size)*sizeof(type) / sizeof(oct_pack);\
//			swInt restCount = \
//				( ( (size)*sizeof(type) ) % sizeof(oct_pack) ) \
//				/ sizeof(type);\
//			swInt j;\
//			for(j = 0; j < packSize; j++)\
//				REG_PUTR(packs[j], destC);\
//			if( restCount != 0 )\
//			{\
//				oct_pack restPack;\
//				for(j = 0; j < restCount; j++)\
//				((type*) &restPack)[j] = \
//					((type*) &packs[packSize])[j];\
//				REG_PUTR(restPack, destC);\
//			}	\
//		}\
//		else\
//		{\
//			swInt srcC = colID + (swInt) (8 >> (i+1)); \
//			swInt packSize = (size)*sizeof(type) / sizeof(oct_pack);\
//			swInt restCount = \
//				( ( (size)*sizeof(type) ) % sizeof(oct_pack) ) \
//				/ sizeof(type);\
//			type* realData = (type*) buffer;\
//			swInt numTypeInPack = sizeof(oct_pack) / sizeof(type);\
//			swInt j, typeCount = 0;\
//			oct_pack pack;\
//			for(j = 0; j < packSize; j++)\
//			{\
//				REG_GETR( pack );\
//				swInt k;\
//				for(k = 0; k < numTypeInPack; k++)\
//					realData[typeCount+k] = ((type*) &pack)[k];\
//				typeCount += numTypeInPack;\
//			}\
//			if(restCount != 0)\
//			{\
//				REG_GETR( pack );\
//				swInt k;\
//				for(k = 0; k < restCount; k++)\
//					realData[typeCount+k] = ((type*) &pack)[k];\
//				typeCount += restCount;\
//			}\
//			if(typeCount != size)\
//			{\
//				LOG("received wrong number of data!\n");\
//				exit(-1);\
//			}\
//		}\
//		ROWSYN;\
//	}\
//\
//	ALLSYN;\
//}

#endif // OCT_ALLREDUCE_H
