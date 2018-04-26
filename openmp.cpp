#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "do_serial.h"
#include "omp.h" 

int main(int argc, char **argv) {
	int navg, nabsavg = 0;
	double davg, dmin, absmin = 1.0, absavg = 0.0;

	if (find_option(argc, argv, "-h") >= 0) {
		printf("Options:\n");
		printf("-h to see this help\n");
		printf("-n <int> to set the number of particles\n");
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
	set_size(n);
	int threads = omp_get_num_threads();
	init_particles(n, particles);

	//
	//  simulate a number of time steps
	//
	double simulation_time = read_timer();

	const int BUCKET_COUNT = 100;
	std::vector <std::vector<particle_t> > particle_vec; 
	std::vector <std::vector<particle_t> > ghost_vec; 

	for (int step = 0; step < NSTEPS; step++) {
		navg = 0;
		davg = 0.0;
		dmin = 1.0;	
		create_buckets(particles, n, threads, particle_vec, ghost_vec);

		//printf("Running %d threads on step %d\n", threads, step);
		#pragma omp parallel for reduction (+:navg) reduction (+:davg)
		for (int i = 0; i < threads; i++){
			//printf("i = %d, particle_vec size = %d ghost_vec size = %d \n", i, particle_vec[i].size(), ghost_vec[i].size());
	 		//do_serial_process(particles, n, BUCKET_COUNT, dmin, davg, navg);
			do_serial_process(particle_vec[i].data(), particle_vec[i].size(), BUCKET_COUNT, dmin, davg, navg);

			//printf("davg = %f, navg = %d\n", davg, navg);

			//  compute ghost forces
			for (int ii = 0; ii < particle_vec[i].size(); ii++) {
				particle_vec[i].data()[ii].ax = particle_vec[i].data()[ii].ay = 0.0;
				// Apply from ghosts
				for (int j = 0; j < ghost_vec[i].size(); j++) {
					apply_force(particle_vec[i].data()[ii], ghost_vec[i].data()[j], &dmin, &davg, &navg);
				}
			}
			//  move particles
			for (int iii = 0; iii < particle_vec[i].size(); iii++) {
				acc_move(particle_vec[i].data()[iii]);
			}

			//printf("davg = %f, navg = %d\n", davg, navg);

		}

		//printf("davg = %f, navg = %d\n", davg, navg);

		int counter = 0; 
		for (int curr_bucket =0; curr_bucket < (int)particle_vec.size(); curr_bucket++){
			memcpy(&particles[counter], particle_vec[curr_bucket].data(), particle_vec[curr_bucket].size() *sizeof(particle_t));
			counter += particle_vec[curr_bucket].size(); 
		}
		

		if (find_option(argc, argv, "-no") == -1) {

			//
			// Computing statistical data
			//

			if (navg) {
				absavg += davg / navg;
				nabsavg++;
			}

			if (dmin < absmin)
				absmin = dmin;

			//
			//  save if necessary
			//
			if (fsave && (step % SAVEFREQ) == 0)
				save(fsave, n, particles);
		}
	}

	simulation_time = read_timer() - simulation_time;

	printf("n = %d, simulation time = %g seconds", n, simulation_time);

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
		if (absmin < 0.4) {
			printf(
					"\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
		}
		if (absavg < 0.8) {
			printf(
					"\nThe average distance is below 0.8 meaning that most particles are not interacting");
		}
	}
	printf("\n");

	//
	// Printing summary data
	//
	if (fsum)
		fprintf(fsum, "%d %g\n", n, simulation_time);

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
