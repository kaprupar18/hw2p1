#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "common.h"
#include "mkl.h"

int create_buckets(particle_t *particles, int particle_count, int bucket_count){
	for(int i = 0; i < bucket_count; i++){

	}
}
double split_bucket(std::vector < std::vector<particle_t> > &particles, std::vector < std::vector<particle_t> > &ghosts)
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}
