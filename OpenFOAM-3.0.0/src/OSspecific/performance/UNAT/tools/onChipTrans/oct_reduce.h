#ifndef OCT_REDUCE_H
#define OCT_REDUCE_H

#define CYCLE_OFFSET_8( id, offset) \
{ \
	(id) = ( (id) - (offset) + 8) % 8; \
}

#define oct_reduce( buffer, type, size, opt, id) \
{ \
	if(  sizeof(oct_pack) < sizeof(type) ||\
				sizeof(oct_pack) % sizeof(type) != 0 )\
	{\
		if(_MYID == 0)\
			LOG("size of oct_pack is not integeral times of the size of type!\n");\
		exit(-1);\
	}\
	if( id > 63 ) \
	{\
		if(_MYID == 0)\
			LOG("wrong core id to do reduce!\n");\
		exit(-1);\
	}\
\
	swInt mid_row = (swInt) ( (id) >> 3);\
	swInt mid_col = (swInt) ( id%8);\
	swInt rowID = (swInt) (_MYID >> 3);\
	swInt colID = (swInt) (_MYID%8);\
	CYCLE_OFFSET_8(rowID, mid_row+1);\
	CYCLE_OFFSET_8(colID, mid_col+1);\
\
	swInt i;\
	for(i = 0; i < 3; i++)\
	{\
		if( rowID > 7 - ( 8 >> i ) )\
			if( ( rowID >> (2-i) ) % 2 == 0 )\
			{\
				swInt destR = rowID + (swInt) ( 8 >> (i+1) );\
				CYCLE_OFFSET_8(destR, -(mid_row+1));	\
				oct_pack* packs = (oct_pack*) (buffer);\
				swInt packSize = (size)*sizeof(type) / sizeof(oct_pack);\
				swInt restCount = \
					( ( (size)*sizeof(type) ) % sizeof(oct_pack) ) \
					/ sizeof(type);\
				swInt j;\
				for(j = 0; j < packSize; j++)\
					REG_PUTC(packs[j], destR);\
				if( restCount != 0 )\
				{\
					oct_pack restPack;\
					for(j = 0; j < restCount; j++)\
						((type*) &restPack)[j] = \
							((type*) &packs[packSize])[j];\
					REG_PUTC(restPack, destR);\
				}	\
			}\
			else\
			{\
				swInt srcR = rowID - (swInt) (8 >> (i+1)); \
				CYCLE_OFFSET_8(srcR, -(mid_row+1));	\
				swInt packSize = (size)*sizeof(type) / sizeof(oct_pack);\
				swInt restCount = \
					( ( (size)*sizeof(type) ) % sizeof(oct_pack) ) \
					/ sizeof(type);\
				type* realData = (type*) buffer;\
				swInt numTypeInPack = sizeof(oct_pack) / sizeof(type);\
				swInt j, typeCount = 0;\
				oct_pack pack;\
				for(j = 0; j < packSize; j++)\
				{\
					REG_GETC( pack );\
					swInt k;\
					for(k = 0; k < numTypeInPack; k++)\
						opt(realData[typeCount+k], ((type*) &pack)[k]);\
					typeCount += numTypeInPack;\
				}\
				if(restCount != 0)\
				{\
					REG_GETC( pack );\
					swInt k;\
					for(k = 0; k < restCount; k++)\
						opt(realData[typeCount+k], ((type*) &pack)[k]);\
					typeCount += restCount;\
				}\
				if(typeCount != size)\
				{\
					LOG("received wrong number of data!\n");\
					printf("typeCount is %d , size is %d\n", typeCount, size);\
					exit(-1);\
				}\
			}\
\
		COLSYN;\
	}\
\
	ALLSYN;\
\
	for(i = 0; i < 3; i++)\
	{\
		if( rowID = 7 && colID > 7 - (8 >> i) )\
			if( ( colID >> (2-i) ) % 2 == 0 )\
			{\
				swInt destC = colID + (swInt) ( 8 >> (i+1) ); \
				CYCLE_OFFSET_8(destC, -(mid_col+1));	\
				oct_pack* packs = (oct_pack*) (buffer);\
				swInt packSize = (size)*sizeof(type) / sizeof(oct_pack);\
				swInt restCount = \
					( ( (size)*sizeof(type) ) % sizeof(oct_pack) ) \
					/ sizeof(type);\
				swInt j;\
				for(j = 0; j < packSize; j++)\
					REG_PUTR(packs[j], destC);\
				if( restCount != 0 )\
				{\
					oct_pack restPack;\
					for(j = 0; j < restCount; j++)\
					((type*) &restPack)[j] = \
						((type*) &packs[packSize])[j];\
					REG_PUTR(restPack, destC);\
				}	\
			}\
			else\
			{\
				swInt srcC = colID - (swInt) (8 >> (i+1)); \
				CYCLE_OFFSET_8(srcC, -(mid_col+1));	\
				swInt packSize = (size)*sizeof(type) / sizeof(oct_pack);\
				swInt restCount = \
					( ( (size)*sizeof(type) ) % sizeof(oct_pack) ) \
					/ sizeof(type);\
				type* realData = (type*) buffer;\
				swInt numTypeInPack = sizeof(oct_pack) / sizeof(type);\
				swInt j, typeCount = 0;\
				oct_pack pack;\
				for(j = 0; j < packSize; j++)\
				{\
					REG_GETR( pack );\
					swInt k;\
					for(k = 0; k < numTypeInPack; k++)\
						opt(realData[typeCount+k], ((type*) &pack)[k]);\
					typeCount += numTypeInPack;\
				}\
				if(restCount != 0)\
				{\
					REG_GETR( pack );\
					swInt k;\
					for(k = 0; k < restCount; k++)\
						opt(realData[typeCount+k], ((type*) &pack)[k]);\
					typeCount += restCount;\
				}\
				if(typeCount != size)\
				{\
					LOG("received wrong number of data!\n");\
					exit(-1);\
				}\
			}\
		\
		ROWSYN;\
	}\
\
	ALLSYN;\
}

#endif //OCT_REDUCE_H
