#ifdef _OPENMP
#include <omp.h>
#endif

#include <ctime>

inline double get_time(void){
#ifdef _OPENMP
    return omp_get_wtime();
#else
    clock_t time = clock();
    return ((double) time / CLOCKS_PER_SEC);
#endif
}

inline void set_num_threads(int numthreads){
#ifdef _OPENMP
    omp_set_dynamic(0);
        omp_set_num_threads(numthreads);
#endif
}

inline int get_thread_id(void){
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}
