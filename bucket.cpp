//#include <stdlib.h>
//#include <stdio.h>
//#include <assert.h>
//#include <float.h>
//#include <string.h>
//#include <math.h>
//#include <time.h>
//#include <sys/time.h>
//#include <vector>
#include "common.cpp"
#include "mkl.h"

int split_bucket(particle_t *all_particles, int particle_count,
		std::vector<std::vector<particle_t> > &particles,
		std::vector<std::vector<particle_t> > &ghosts, int index_to_split) {
	double sample_percent = 0.02;

	int particle_subsamp = particles[index_to_split].size() * sample_percent;

	double min_x = 9.0e9;
	double max_x = -9.0e9;
	double min_y = 9.0e9;
	double max_y = -9.0e9;
	double mean_x = 0.0;
	double mean_y = 0.0;
	double std_x = 0.0;
	double std_y = 0.0;
	for (int i = 0; i < particle_subsamp; i++) {
		mean_x += particles[index_to_split][i].x;
		mean_y += particles[index_to_split][i].y;

	}
	mean_x = mean_x / particle_subsamp;
	mean_y = mean_y / particle_subsamp;

	for (int i = 0; i < particle_subsamp; i++) {
		std_x += fabs(mean_x - particles[index_to_split][i].x);
		std_y += fabs(mean_y - particles[index_to_split][i].y);
	}

	std_x = std_x / particle_subsamp;
	std_y = std_y / particle_subsamp;

	// Split actual particles, add two new vectors, we will be deleting the original element later
	std::vector<particle_t> new_bucket_1;
	std::vector<particle_t> new_bucket_2;
	particles.push_back(new_bucket_1);
	particles.push_back(new_bucket_2);

	// The ghosts are the obnoxious part, you gotta look through all the particles... there must be a better way
	std::vector<particle_t> new_ghosts_1;
	std::vector<particle_t> new_ghosts_2;
	ghosts.push_back(new_ghosts_1);
	ghosts.push_back(new_ghosts_2);

	if (std_x >= std_y) {

		for (std::vector<particle_t>::iterator curr_part =
				particles[index_to_split].begin();
				curr_part != particles[index_to_split].end(); ++curr_part) {
			if (curr_part->x < mean_x) {
				particles.rbegin()[1].push_back(*curr_part);
			} else {
				particles.rbegin()[0].push_back(*curr_part);
			}
			//printf("minx=%lf,curr_part->x=%lf | ", min_x, curr_part->x);
			min_x = dmin(min_x, curr_part->x);
			min_y = dmin(min_y, curr_part->y);
			max_x = dmax(max_x, curr_part->x);
			max_y = dmax(max_y, curr_part->y);
		}
		printf("Bucket Split Stats on X : minx=%lf, maxx=%lf, miny=%lf, maxy=%lf, meanx=%lf, meany=%lf, stdx=%lf, stdy=%lf \n",
							min_x, max_x, min_y, max_y, mean_x, mean_y, std_x, std_y);

		for (int i = 0; i < particle_count; i++) {
			if ((all_particles[i].x > min_x - cutoff
					&& all_particles[i].y > min_y - cutoff)
					&& (all_particles[i].x < mean_x + cutoff
							&& all_particles[i].y < max_y + cutoff)
					&& (all_particles[i].x < min_x
							|| all_particles[i].x > mean_x
							|| all_particles[i].y < min_y
							|| all_particles[i].y > max_y)) {

				ghosts.rbegin()[1].push_back(all_particles[i]);

			} else if ((all_particles[i].x > mean_x - cutoff
					&& all_particles[i].y > min_y - cutoff)
					&& (all_particles[i].x < max_x + cutoff
							&& all_particles[i].y < max_y + cutoff)
					&& (all_particles[i].x < mean_x
							|| all_particles[i].x > max_x
							|| all_particles[i].y < min_y
							|| all_particles[i].y > max_y)) {

				ghosts.rbegin()[0].push_back(all_particles[i]);
			}
		}
	} else {
		// In this scenario there is higher variance in the y direction so we split it that way
		for (std::vector<particle_t>::iterator curr_part =
				particles[index_to_split].begin();
				curr_part != particles[index_to_split].end(); ++curr_part) {
			if (curr_part->y < mean_y) {
				particles.rbegin()[1].push_back(*curr_part);
			} else {
				particles.rbegin()[0].push_back(*curr_part);
			}

			min_x = dmin(min_x, curr_part->x);
			min_y = dmin(min_y, curr_part->y);
			max_x = dmax(max_x, curr_part->x);
			max_y = dmax(max_y, curr_part->y);
		}
		printf("Bucket Split Stats on Y : minx=%lf, maxx=%lf, miny=%lf, maxy=%lf, meanx=%lf, meany=%lf, stdx=%lf, stdy=%lf \n",
							min_x, max_x, min_y, max_y, mean_x, mean_y, std_x, std_y);

		for (int i = 0; i < particle_count; i++) {
			if ((all_particles[i].x > min_x - cutoff
					&& all_particles[i].y > min_y - cutoff)
					&& (all_particles[i].x < max_x + cutoff
							&& all_particles[i].y < mean_y + cutoff)
					&& (all_particles[i].x < min_x || all_particles[i].x > max_x
							|| all_particles[i].y < min_y
							|| all_particles[i].y > mean_y)) {

				ghosts.rbegin()[1].push_back(all_particles[i]);

			} else if ((all_particles[i].x > min_x - cutoff
					&& all_particles[i].y > min_y - cutoff)
					&& (all_particles[i].x < max_x + cutoff
							&& all_particles[i].y < max_y + cutoff)
					&& (all_particles[i].x < min_x || all_particles[i].x > max_x
							|| all_particles[i].y < min_y
							|| all_particles[i].y > max_y)) {

				ghosts.rbegin()[0].push_back(all_particles[i]);
			}
		}
	}

	// Get rid of the old complete particles and ghosts
	particles.erase(particles.begin() + index_to_split);
	ghosts.erase(ghosts.begin() + index_to_split);

	return 0;

}

int create_buckets(particle_t *particles, int particle_count, int bucket_count,
		std::vector<std::vector<particle_t> > &particle_vec,
		std::vector<std::vector<particle_t> > &ghost_vec) {

	// Initialize vectors
	particle_vec.resize(1);
	ghost_vec.resize(1);

	particle_vec[0].resize(particle_count);
	ghost_vec[0].resize(0);

	// Copy particle pointer to first bucket
	//memcpy(&x[0], source, particle_count*sizeof(particle_t));
	particle_vec[0].assign(particles, particles + particle_count);

	for (int i = 1; i < bucket_count; i++) {
		// Find the biggest bucket
		int max_bucket = 0;
		int max_bucket_count = 0;
		for (int j = 0; j < (int)particle_vec.size(); j++) {
			if ((int) particle_vec[j].size() > max_bucket_count) {
				max_bucket = j;
				max_bucket_count = particle_vec[j].size();
			}
		}

		printf("Splitting bucket %d of %d with %d particles\n", max_bucket, (int)particle_vec.size(),
				max_bucket_count);
		split_bucket(particles, particle_count, particle_vec, ghost_vec,
				max_bucket);

	}

	return 0;
}

