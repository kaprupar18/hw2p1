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

double size;

//
//  tuned constants
//
#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005

//
//  timer
//
double read_timer() {
	static bool initialized = false;
	static struct timeval start;
	struct timeval end;
	if (!initialized) {
		gettimeofday(&start, NULL);
		initialized = true;
	}
	gettimeofday(&end, NULL);
	return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//  keep density constant
//
void set_size(int n) {
	size = sqrt( density * n);
}

//
//  Initialize the particle positions and velocities
//
void init_particles(int n, particle_t *p) {
	srand48 (time(NULL) );

int	sx = (int)ceil(sqrt((double)n));
	int sy = (n+sx-1)/sx;

	int *shuffle = (int*)malloc( n * sizeof(int) );
	for( int i = 0; i < n; i++ )
	shuffle[i] = i;

	for( int i = 0; i < n; i++ )
	{
		//
		//  make sure particles are not spatially sorted
		//
		int j = lrand48()%(n-i);
		int k = shuffle[j];
		shuffle[j] = shuffle[n-i-1];

		//
		//  distribute particles evenly to ensure proper spacing
		//
		p[i].x = size*(1.+(k%sx))/(1+sx);
		p[i].y = size*(1.+(k/sx))/(1+sy);

		//
		//  assign random velocities within a bound
		//
		p[i].vx = drand48()*2-1;
		p[i].vy = drand48()*2-1;
	}
	free( shuffle );
}

//
//  interact two particles
//
void apply_force(particle_t &particle, particle_t &neighbor, double *dmin,
		double *davg, int *navg) {

	double dx = neighbor.x - particle.x;
	double dy = neighbor.y - particle.y;
	double r2 = dx * dx + dy * dy;
	if (r2 > cutoff * cutoff) {
		return;
	}

	// Since the algorithm will comparit a particle with itself we need to take this into
	// account by having the particle be returned uneffected
	if (r2 == 0.0) {
		return;
	}

	if (r2 / (cutoff * cutoff) < (*dmin) * (*dmin)) {
		*dmin = sqrt(r2) / cutoff;
	}
	(*davg) += sqrt(r2) / cutoff;
	(*navg)++;
	//printf("Applied force a %dth time with davg = %lf \n", (*navg), (*davg));

	// Enable a mimimum radius so we dont have crazy high forces causing instability
	r2 = fmax(r2, min_r * min_r);
	double r = sqrt(r2);

	//
	//  very simple short-range repulsive force
	//
	double coef = (1 - cutoff / r) / r2 / mass;
	particle.ax += coef * dx;
	particle.ay += coef * dy;
}
void opapply_force(particle_t *particle, particle_t *neighbor, double *dmin,
		double *davg, int *navg) {

	double dx = neighbor->x - particle->x;
	double dy = neighbor->y - particle->y;
	// printf("Decimals: %f , \n",dy);
	double r2 = dx * dx + dy * dy;

	if (r2 > cutoff * cutoff)
		return;

	if (r2 != 0) {
		if (r2 / (cutoff * cutoff) < *dmin * (*dmin))
			*dmin = sqrt(r2) / cutoff;
		(*davg) += sqrt(r2) / cutoff;
		(*navg)++;
	}

	r2 = fmax(r2, min_r * min_r);
	double r = sqrt(r2);

	//
	//  very simple short-range repulsive force
	//
	double coef = (1 - cutoff / r) / r2 / mass;
	particle->ax += coef * dx;
	particle->ay += coef * dy;

}
//
//  integrate the ODE
//
void move(particle_t &p) {
	//
	//  slightly simplified Velocity Verlet integration
	//  conserves energy better than explicit Euler method
	//
	p.vx += p.ax * dt;
	p.vy += p.ay * dt;
	p.x += p.vx * dt;
	p.y += p.vy * dt;

	//
	//  bounce from walls
	//
	if (p.x < 0 || p.x > size) {
		p.x = p.x < 0 ? -p.x : 2 * size - p.x;
		p.vx = -p.vx;
	}
	if (p.y < 0 || p.y > size) {
		p.y = p.y < 0 ? -p.y : 2 * size - p.y;
		p.vy = -p.vy;
	}
}
//void move( particle_t &p_in, particle_t &p_out )
//{
//    //
//    //  slightly simplified Velocity Verlet integration
//    //  conserves energy better than explicit Euler method
//    //
//	p_out.vx += p_in.ax * dt;
//	p_out.vy += p_in.ay * dt;
//	p_out.x  += p_in.vx * dt;
//	p_out.y  += p_in.vy * dt;
//
//    //
//    //  bounce from walls
//    //
//    if( p_out.x < 0 || p_out.x > size )
//    {
//    	p_out.x  = p_out.x < 0 ? -p_out.x : 2*size-p_out.x;
//    	p_out.vx = -p_out.vx;
//    }
//    if( p_out.y < 0 || p_out.y > size )
//    {
//    	p_out.y  = p_out.y < 0 ? -p_out.y : 2*size-p_out.y;
//        p_out.vy = -p_out.vy;
//    }
//}

//
//	state update
//
//void stateUpdate( double *p_in, double *p_out ){
//	cblas_dgemm('C', CblasNoTrans, CblasNoTrans,
//	                m, n, k, alpha, A, k, B, n, beta, C, n);
//}

//
//  I/O routines
//
void save(FILE *f, int n, particle_t *p) {
	static bool first = true;
	if (first) {
		fprintf(f, "%d %g\n", n, size);
		first = false;
	}
	for (int i = 0; i < n; i++)
		fprintf(f, "%g %g\n", p[i].x, p[i].y);
}

//
//  command line option processing
//
int find_option(int argc, char **argv, const char *option) {
	for (int i = 1; i < argc; i++)
		if (strcmp(argv[i], option) == 0)
			return i;
	return -1;
}

int read_int(int argc, char **argv, const char *option, int default_value) {
	int iplace = find_option(argc, argv, option);
	if (iplace >= 0 && iplace < argc - 1)
		return atoi(argv[iplace + 1]);
	return default_value;
}

char *read_string(int argc, char **argv, const char *option,
		char *default_value) {
	int iplace = find_option(argc, argv, option);
	if (iplace >= 0 && iplace < argc - 1)
		return argv[iplace + 1];
	return default_value;
}
// Bucket
void setbounds(bucket *p, double sx, double sy, double ex, double ey) //  sets boundaries
		{
	p->sx = sx;
	p->sy = sy;
	p->ex = ex;
	p->ey = ey;

}
;
void addparticle(bucket *p, particle_t *k) {
	p->arr.push_back(k);
	p->count = p->count + 1;
}
;
void deleteparticle(bucket *p, particle_t *k) {
	for (unsigned int i = 0; i < (p->arr.size()); i++) {
		if (k == p->arr[i]) {
			p->arr.erase(p->arr.begin() + i);
			p->count = p->count - 1;
			break;
		}
	}
}
;

double sumx(bucket *p) {
	double allx = 0.0; // av[1] = x average
	double x;
	for (unsigned int i = 0; i < p->arr.size(); i++) {
		x = p->arr[i]->x;
		allx += x;
	}
	return allx;
}
;
double sumy(bucket *p) {
	double ally = 0.0; // av[1] = x average
	double y;
	for (unsigned int i = 0; i < p->arr.size(); i++) {
		y = p->arr[i]->y;
		ally += y;
	}
	return ally;
}
;
void insort(bucket *p, particle_t *t, int size) {
	double sx = p->sx;
	double sy = p->sy;
	double ex = p->ex;
	double ey = p->ey;
	for (int i = 0; i < size; i++) {
		if (t[i].x >= sx && t[i].x <= ex && t[i].y >= sy && t[i].y <= ey) {
			addparticle(p, &t[i]);
		}
	}

}
;

