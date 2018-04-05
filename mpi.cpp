#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
//#include <vector>
#include "bucket.cpp"

//
//  benchmarking program
//
int main(int argc, char **argv) {
	int navg = 0;
	int nabsavg = 0;
	double dmin = 1.0;
	double absmin = 1.0;
	double davg = 0.0; 
	double absavg = 0.0;
	double rdavg, rdmin;
	int rnavg = 0;

	//
	//  process command line parameters
	//
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

	//
	//  set up MPI
	//
	int n_proc, rank, num_slaves;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &n_proc); // Total number of processors available
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Current processor number
	num_slaves = n_proc - 1;

	//
	//  allocate generic resources
	//
	FILE *fsave = savename && rank == 0 ? fopen(savename, "w") : NULL;
	FILE *fsum = sumname && rank == 0 ? fopen(sumname, "a") : NULL;

	MPI_Datatype PARTICLE;
	MPI_Type_contiguous(6, MPI_DOUBLE, &PARTICLE);
	MPI_Type_commit(&PARTICLE);
	double simulation_time = 0.0;

	//****************  MASTER  ****************
	if (rank == 0) {
		//printf("Starting master\n");
		std::vector < std::vector<particle_t> > particle_vec;
		std::vector < std::vector<particle_t> > ghost_vec;

		//  initialize and distribute the particles (that's fine to leave it unoptimized)
		particle_t *particles = (particle_t*) malloc(n * sizeof(particle_t));
		set_size(n);
		init_particles(n, particles);

		//
		//  simulate a number of time steps
		//
		simulation_time = read_timer();
		for (int step = 0; step < NSTEPS; step++) {

			if (find_option(argc, argv, "-no") == -1)
				if (fsave && (step % SAVEFREQ) == 0)
					save(fsave, n, particles);

			// STEP - 1
			create_buckets(particles, n, num_slaves, particle_vec, ghost_vec);
			if(num_slaves != 1){
				//printf("Created %d buckets for %d slaves\n", (int)particle_vec.size(), num_slaves);
			}

			// STEP - 2
			// Send all vectors to their respective processes
			for (int curr_slave = 1; curr_slave < n_proc; curr_slave++) {
				// Set particle vec
				if(num_slaves != 1){
					//printf("Sending %d particles to slave %d\n", (int)particle_vec[curr_slave - 1].size(), curr_slave);
				}
				MPI_Send(particle_vec[curr_slave - 1].data(),
						particle_vec[curr_slave - 1].size(), PARTICLE,
						curr_slave, 1, MPI_COMM_WORLD);
			}

			// STEP - 3
			// Send all ghost vectors
			for (int curr_slave = 1; curr_slave < n_proc; curr_slave++) {
				// Set particle vec
				MPI_Send(ghost_vec[curr_slave - 1].data(),
						ghost_vec[curr_slave - 1].size(), PARTICLE, curr_slave,
						2, MPI_COMM_WORLD);
			}

			// STEP - 4
			// Receive all completed particles
			int counter = 0;
			for (int curr_slave = 1; curr_slave < n_proc; curr_slave++) {
				// Receive new particles
				MPI_Recv(&particles[counter],
						particle_vec[curr_slave - 1].size(), PARTICLE,
						curr_slave, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				counter += particle_vec[curr_slave - 1].size();
				//printf("Receive %d particles from slave %d\n", (int)particle_vec[curr_slave - 1].size(), curr_slave);
			}

		}

		simulation_time = read_timer() - simulation_time;

	} else		//****************  SLAVES  ****************
	{
		printf("Starting slave %d\n", rank);
		particle_t *local = (particle_t*) malloc(n * sizeof(particle_t));
		particle_t *local_ghost = (particle_t*) malloc(n * sizeof(particle_t));

		MPI_Status status;
		int local_count = 0;
		int local_ghost_count = 0;

		for (int step = 0; step < NSTEPS; step++) {
			navg = 0;
			dmin = 1.0;
			davg = 0.0;
//			printf("Slave %d, step %d\n", rank, step);

			// STEP - 1
			// Receive all the main particles
			MPI_Recv(local, n, PARTICLE, 0, 1, MPI_COMM_WORLD, &status);
//			printf("Slave %d, recieved local particles\n", rank);
			MPI_Get_count(&status, PARTICLE, &local_count);
//			printf("Slave %d, recieved local particles count = %d\n", rank, local_count);

			// Receive all the main particles
			MPI_Recv(local_ghost, n, PARTICLE, 0, 2, MPI_COMM_WORLD, &status);
			MPI_Get_count(&status, PARTICLE, &local_ghost_count);

			//
			//  compute all forces
			//
			for (int i = 0; i < local_count; i++) {
				local[i].ax = local[i].ay = 0;
				// Apply from fellow local
				for (int j = 0; j < local_count; j++) {
					apply_force(local[i], local[j], &dmin, &davg, &navg);
				}
				// Apply from ghosts
				for (int j = 0; j < local_ghost_count; j++) {
					apply_force(local[i], local_ghost[j], &dmin, &davg, &navg);
				}
			}

			//
			//  move particles
			//
			for (int i = 0; i < local_count; i++) {
				move(local[i]);
			}

			if (find_option(argc, argv, "-no") == -1) {

				//
				// Computing statistical data
				//
				if (rnavg) {
					absavg += rdavg / rnavg;
					nabsavg++;
				}
				if (rdmin < absmin)
					absmin = rdmin;

			}

			MPI_Send(local, local_count, PARTICLE, 0, 3, MPI_COMM_WORLD);
		}
		
		// Each slave is responsible for its own
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

	}

	if (rank == 0) {
		printf("n = %d, simulation time = %g seconds", n, simulation_time);


		//
		// Printing summary data
		//
		if (fsum)
			fprintf(fsum, "%d %d %g\n", n, n_proc, simulation_time);
	}

	//
	//  release resources
	//
	if (fsum)
		fclose(fsum);
	if (fsave)
		fclose(fsave);

	MPI_Finalize();

	return 0;
}
