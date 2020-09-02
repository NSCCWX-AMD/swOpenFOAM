#include "util.h"
#include <stdio.h>
//LABEL binarySearch(LABEL *arr, LABEL start, LABEL end, LABEL target) ;//supremum

void get_block_row_end_index(const LABEL *arr , LABEL *res , LABEL size , LABEL block_row_size,LABEL offset_size) {

    LABEL i;
    res[0] = 0 ;
    for(i = 0 ; i < block_row_size ; i++ ) {
         res[i+1] = binarySearch(arr,res[i],size - 1,(i+1)<<offset_size,offset_size);
    }
}

LABEL binarySearch(const LABEL *arr, LABEL start, LABEL end, LABEL target,LABEL offset_size) {//supremum
    LABEL mid ;
    if(!LESS_TARGET_OFFSET(arr[start],target,offset_size)) {
        mid = start ;
    }
    else {
        if(LESS_TARGET_OFFSET(arr[end],target,offset_size)) {
            mid = end + 1;
        }
        else {
            mid = (start + end) / 2;
            while(start < end ) {
          //      if (arr[mid] == target && arr[mid] != target) {
                if(LESS_TARGET_OFFSET(arr[mid],target,offset_size)){

                    if(LESS_TARGET_OFFSET(arr[mid+1],target,offset_size)) {
                       start = mid + 1;
                    }
                    else {
                //     prLABELf("has find %d ",mid);
                        break;

                    }
                }
                else {
                     end = mid - 1;
                }
                mid = (end + start) / 2;
            }
            mid++;
        }
    }
    return mid;
}


