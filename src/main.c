#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

#include "timeit.h"
#include "array.h"
#include "arraymath.h"

#define M 255
#define N 255
#define REPS 30

int main(){
    const double dRAND_MAX = (double)RAND_MAX;

    clock_t t0, t1;

    TIMEIT("Allocating Array", REPS,
        Arr2D arr1 = empty_Arr2D(M, N);
        free_Arr2D(&arr1);
    );
    Arr2D arr1 = empty_Arr2D(M, N);

    TIMEIT("Filling Array with random values", REPS,
        for (int i = 0; i < M; i++){
            for (int j = 0; j < N; j++){
                double x = rand() / dRAND_MAX;
                set_Arr2D(&arr1, i, j, x);
            }
        }
    );
    for (int i = 0; i < M; i++){
            for (int j = 0; j < N; j++){
                double x = rand() / dRAND_MAX;
                set_Arr2D(&arr1, i, j, x);
            }
        }

    TIMEIT("Copying Array", REPS,
        Arr2D arr2 = copy_Arr2D(arr1);
        free_Arr2D(&arr2);
    );
    
    TIMEIT("Adding Arrays", REPS,
        Arr2D arr2 = copy_Arr2D(arr1);
        Arr2D result = empty_Arr2D(M, N);
        add_Arr2D(arr1, arr2, &result);
        free_Arr2D(&arr2);
        free_Arr2D(&result);
    );
    TIMEIT("Adding Arrays By Reference", REPS,
        Arr2D arr2 = copy_Arr2D(arr1);
        Arr2D result = empty_Arr2D(M, N);
        add_Arr2D_by_ref(&arr1, &arr2, &result);
        free_Arr2D(&arr2);
        free_Arr2D(&result);
    );


    // for (int i = 0; i < M; i++){
    //     printf("[");
    //     for (int j = 0; j < N; j++){
    //         printf("%.2f ", get_Arr2D(arr1, i, j));
    //     }
    //     printf("]\n");
    // }
    // printf("\n\n\n");

    free_Arr2D(&arr1);
    return 0;
}