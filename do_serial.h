#include "bucket.cpp"

//
//  benchmarking program
//
void do_serial_process(particle_t *particles_ptr, int n_particles,
		int BUCKET_COUNT, double &dmin, double &davg, int &navg) {

	std::vector < std::vector<particle_t> > particle_vec;
	std::vector < std::vector<particle_t> > ghost_vec;
	particle_t *local = (particle_t*) malloc(n_particles * sizeof(particle_t));
	particle_t *local_ghost = (particle_t*) malloc(n_particles * sizeof(particle_t));

	// Create bucket
	create_buckets(particles_ptr, n_particles, BUCKET_COUNT, particle_vec, ghost_vec);
	

	//
	//  compute forces
	//
	for (int curr_bucket = 0; curr_bucket < (int) particle_vec.size();
			curr_bucket++) {
		local = particle_vec[curr_bucket].data();
		local_ghost = ghost_vec[curr_bucket].data();
		for (int i = 0; i < (int) particle_vec[curr_bucket].size(); i++) {
			local[i].ax = local[i].ay = 0;
			// Apply from fellow local
			for (int j = 0; j < (int) particle_vec[curr_bucket].size(); j++) {
				apply_force(local[i], local[j], &dmin, &davg, &navg);
				//printf("davg = %f, navg = %d\n", davg, navg);
			}
			// Apply from ghosts
			for (int j = 0; j < (int) ghost_vec[curr_bucket].size(); j++) {
				apply_force(local[i], local_ghost[j], &dmin, &davg, &navg);
			}
		}
	}

	//
	// Receive all completed particles
	//
	int counter = 0;
	for (int curr_bucket = 0; curr_bucket < (int) particle_vec.size();
			curr_bucket++) {
		memcpy(&particles_ptr[counter], particle_vec[curr_bucket].data(),
				particle_vec[curr_bucket].size() * sizeof(particle_t));
		counter += particle_vec[curr_bucket].size();
		//printf("Receive %d particles from slave %d\n", (int)particle_vec[curr_slave - 1].size(), curr_slave);
	}
	//
	//  move particles
	//
	for (int i = 0; i < n_particles; i++) {
		move (particles_ptr[i]);
	}
}
