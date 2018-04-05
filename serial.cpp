#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "bucket.cpp"

//
//  benchmarking program
//
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

	particle_t *local = (particle_t*) malloc(n * sizeof(particle_t));
	particle_t *local_ghost = (particle_t*) malloc(n * sizeof(particle_t));

	particle_t *particles = (particle_t*) malloc(n * sizeof(particle_t));
	set_size(n);
	init_particles(n, particles);

	//
	//  simulate a number of time steps
	//
	double simulation_time = read_timer();

	std::vector < std::vector<particle_t> > particle_vec;
	std::vector < std::vector<particle_t> > ghost_vec;

	const int BUCKET_COUNT = 50;
	for (int step = 0; step < NSTEPS; step++) {
		navg = 0;
		davg = 0.0;
		dmin = 1.0;
		
		// Create bucket
		create_buckets(particles, n, BUCKET_COUNT, particle_vec, ghost_vec);

		//
		//  compute forces
		//
		for(int curr_bucket = 0; curr_bucket < (int)particle_vec.size(); curr_bucket++){
			local = particle_vec[curr_bucket].data();
			local_ghost = ghost_vec[curr_bucket].data();
			for (int i = 0; i < (int)particle_vec[curr_bucket].size(); i++) {
				local[i].ax = local[i].ay = 0;
				// Apply from fellow local
				for (int j = 0; j < (int)particle_vec[curr_bucket].size(); j++) {
					apply_force(local[i], local[j], &dmin, &davg, &navg);
				}
				// Apply from ghosts
				for (int j = 0; j < (int)ghost_vec[curr_bucket].size(); j++) {
					apply_force(local[i], local_ghost[j], &dmin, &davg, &navg);
				}
			}
		}


		//
		// Receive all completed particles
		//
		int counter = 0;
		for(int curr_bucket = 0; curr_bucket < (int)particle_vec.size(); curr_bucket++){
			memcpy(&particles[counter], particle_vec[curr_bucket].data(), particle_vec[curr_bucket].size() * sizeof(particle_t));
			counter += particle_vec[curr_bucket].size();
			//printf("Receive %d particles from slave %d\n", (int)particle_vec[curr_slave - 1].size(), curr_slave);
		}


		//
		//  move particles
		//
		for (int i = 0; i < n; i++) {
			move(particles[i]);
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
		if (absmin < 0.4){
			printf(
					"\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
		}
		if (absavg < 0.8){
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
