#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

#include <vector>
//#include "common.cpp"
inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

inline double dmin( double a, double b ) { return a < b ? a : b; }
inline double dmax( double a, double b ) { return a > b ? a : b; }

//
//  saving parameters
//
const int NSTEPS = 1000;
const int SAVEFREQ = 10;

//
// particle data structure
//
#pragma pack() // Pack so we know the order of the data
typedef struct 
{
  double x;
  double y;
  double vx;
  double vy;
  double ax;
  double ay;
} particle_t;

typedef struct
{
double sx; //starting of x
double sy; // starting of y
double ex; //end of x
double ey; // end of y
int count;
std::vector<particle_t*> arr; // vector which points to the addresses of the particles
}bucket;

typedef struct
{
	particle_t *particles;
	int numel;
}bucket_t;

//
//
// Buckets routine 
//
void setbounds(bucket *p,double sx,double sy, double ex, double ey); // sets boundaries
void addparticle(bucket *p, particle_t *k);
void deleteparticle(bucket *p, particle_t *k);
double sumx(bucket *p);
double sumy(bucket *p); 
void insort(bucket *p, particle_t *t, int size, bucket *ghost);

//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//
void set_size( int n );
double set_sizep(int n);
void init_particles( int n, particle_t *p );
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg);
void opapply_force( particle_t *particle, particle_t *neighbor , double *dmin, double *davg, int *navg);
void move( particle_t &p );


//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

#endif
