#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "common.cpp"
#include "omp.h"
//#include "bucket.cpp"
//
//  benchmarking program
//
int main(int argc, char **argv) {
	int navg, nabsavg = 0, numthreads;
	double dmin, absmin = 1.0, davg, absavg = 0.0;

	if (find_option(argc, argv, "-h") >= 0) {
		printf("Options:\n");
		printf("-h to see this help\n");
		printf("-n <int> to set number of particles\n");
		printf("-o <filename> to specify the output file name\n");
		printf("-s <filename> to specify a summary file name\n");
		printf("-no turns off all correctness checks and particle output\n");
		return 0;
	}

	int n = read_int(argc, argv, "-n", 1000);
	char *savename = read_string(argc, argv, "-o", NULL);
	char *sumname = read_string(argc, argv, "-s", NULL);

	FILE *fsave = savename ? fopen(savename, "w") : NULL;
	FILE *fsum = sumname ? fopen(sumname, "a") : NULL;

	particle_t *particles = (particle_t*) malloc(n * sizeof(particle_t));
	double size = set_sizep(n);
	init_particles(n, particles); 

	//
	//  simulate a number of time steps
	//
	double simulation_time = read_timer();
    
    

#pragma omp parallel private(dmin) /*shared(particles,p,ghost)*/
	{
		numthreads = omp_get_num_threads();
		for (int step = 0; step < 100; step++) {
			navg = 0;
			davg = 0.0;
			dmin = 1.0;
#pragma omp for reduction (+:navg) reduction(+:davg) 
            		for (int i =0; i<numthreads; i++)
            		{
				//Intialize the Partitions 
				bucket *p = new bucket(); 
				bucket *ghost = new bucket();

                		//printf(" num of threads:%i num of particle:%i \n",numthreads,n); 
                		int r = i+1;
                		p->sx = (size/numthreads)*i; 
                		p->sy = (size/numthreads)*i;
                		p->ex = (size/numthreads)*(r) + .1; 
                		p->ey = (size/numthreads)*(r);
                		p->count = 0; 
                		insort(p, particles, n, ghost);
                		//printf(" pi count: %i \n, bucket name: %i ", p->count,i);
                		for(int k = 0; k < p->arr.size(); k++)
                		{
                    			p->arr[k]->ax = p->arr[k]->ay = 0;
                    			for (int j = 0; j < ghost->arr.size(); j++)
                    			{
                        
					    opapply_force(p->arr[k], ghost->arr[j], &dmin, &davg, &navg); 
                      
                    			}
                		}
				for (int q = 0; q < p->arr.size(); i++) {
                                	move(p->&(arr[q])); 
                        	}

				delete p; 
				delete ghost;    
            		}



			//
			//  compute all forces
			//
/*#pragma omp for reduction (+:navg) reduction(+:davg)
			for (int i = 0; i < n; i++) {
				particles[i].ax = particles[i].ay = 0;
				for (int j = 0; j < n; j++)
					apply_force(particles[i], particles[j], &dmin, &davg, &navg);
			}
*/
			//
			//  move particles
			//
/*#pragma omp for
			for (int i = 0; i < n; i++)
				move(particles[i]);
*/

			if (find_option(argc, argv, "-no") == -1) {
				//
				//  compute statistical data
				//
//#pragma omp master
				if (navg) {
					absavg += davg / navg;
					nabsavg++;
				}

//#pragma omp critical
				if (dmin < absmin)
					absmin = dmin;

				//
				//  save if necessary
				//
//#pragma omp master
				if (fsave && (step % SAVEFREQ) == 0)
					save(fsave, n, particles);
			}
		
		}
	}
	simulation_time = read_timer() - simulation_time;

	printf("n = %d,threads = %d, simulation time = %g seconds", n, numthreads,simulation_time);

	if (find_option(argc, argv, "-no") == -1) {
		if (nabsavg)
			absavg /= nabsavg;
		//
		//  -the minimum distance absmin between 2 particles during the run of the simulation
		//  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
		//  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
		//
		//  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
		//
		printf(", absmin = %lf, absavg = %lf", absmin, absavg);
		if (absmin < 0.4)
			printf(
					"\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
		if (absavg < 0.8)
			printf(
					"\nThe average distance is below 0.8 meaning that most particles are not interacting");
	}
	printf("\n");

	//
	// Printing summary data
	//
	if (fsum)
		fprintf(fsum, "%d %d %g\n", n, numthreads, simulation_time);

	//
	// Clearing space
	//
	if (fsum)
		fclose(fsum);

	free(particles);
	if (fsave)
		fclose(fsave);

	return 0;
}

