#include "streamgraph.c"
#include "timeit.c"
#include "returncodes.h"
#include <stdbool.h>
#include <stdint.h>
#include <math.h>


int main(){
    bool arr1[1001];
    TIMEIT("bool array", 100, 
    for (int i = 0; i < 1001; i++){
        arr1[i] = false;
    }
    );
    
    bool *arr = malloc(1001*sizeof(bool));
    if (arr == NULL) return NULL_POINTER_ERROR;
    TIMEIT("bool array malloc", 100, 
    for (int i = 0; i < 1001; i++){
        arr[i] = false;
    }
    );
    



    const s = sizeof(uint8_t);
    int arrlen = (int)div(1001, s) + 1001 % s;
    uint8_t arr2[arrlen];
    // TIMEIT("uint8 bool array", 100, 
    // for (int i = 0; i < arrlen; i++){
    //     for (int j = 0; j < 8; j++){
    //         arr2[i] <<= 0;
    //     }
    // }
    // );
    return SUCCESS;
}