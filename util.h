#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <fstream>
#include <cstring>
#include <sys/stat.h>

using namespace std;

#ifndef PI
#define PI (3.14159265358979324)
#endif

#ifndef hbarc
#define hbarc (0.1973)
#endif

#ifndef yes
#define yes 1
#endif

#ifndef no 
#define no 0
#endif

#ifndef true 
#define true 1
#endif

#ifndef false 
#define false 0
#endif

#ifndef default_tol
#define default_tol (1.0e-8)
#endif

#define absol(a) ((a) > 0 ? (a) : (-(a)))
#define maxi(a, b) ((a) > (b) ? (a) : (b))
#define mini(a, b) ((a) < (b) ? (a) : (b))
#define sgn(x) ((x) < 0.0 ? (-1.0) : (1.0))
#define theta(x) ((x) < 0.0 ? (0.0) : (1.0))
#define SQR(a) ((a)*(a))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define DMAX(a,b) ((a) > (b) ? (a) : (b))
#define DMIN(a,b) ((a) < (b) ? (a) : (b))
#define SIGN(a,b) ((b) >= 0.0 ? absol(a): -absol(a)) 


class Util{
 public:
double ***cube_malloc(int , int , int ); 
double **mtx_malloc(int , int );
double *vector_malloc(int );
char **char_mtx_malloc(int , int );
int **int_mtx_malloc(int , int );
char *char_malloc(int );
void char_free(char *);
int *int_malloc(int );

void cube_free(double ***, int, int, int);
void mtx_free(double **, int, int);
void vector_free(double *);

double Power(double, int);
double sech(double x);
void itoa(int , char *);
int is_yes_no(char *);
void reverse( char *);
void prterr(char *);

int integer(double );
int count_lines(char *);
int IsFile(string );

double Dot(double *, double *, int);
double m4p(double *, double *);
double LinearPara(double , double , double , double *);

double DFind(string file_name, const char *st);
string StringFind(string file_name, const char *st);
char *StringFind2(char *file_name, const char *st);
string StringFind3(string file_name, const char *st);
string StringFind4(string file_name, string str_in);
int IFind(string file_name, const char *st);

void FileCopy(const char *in_file, const char *out_file);
void FileCat(const char *in_file, const char *out_file);

void ReWrite(char *file_name, char *st, double x);
void ReWriteString(char *file_name, const char *st, char *x);

int IFindXInVx(double x, double *Vx, int ymax);

void Shout(void);

double Solve 
(double nu, double (*func)(double), 
 double s_down, double s_up, double tol, int *count);

int NumCol(char *infile);
int binning(double x, double xi, double xf, double *bins, int num_bins);
void InBinLocal(double x, double xdown, double xup, double inc, double *num);
double d_read_in(char *s);
int i_read_in(char *s);
char *ch_read_in(char *s);
void CountTime(int i, int inc1, int inc2, int last);
int *heapSort(double *numbers, int array_size);
double Theta(double x, double a);
double lin_int(double x1,double x2,double f1,double f2,double x);

/* nrutil.h */

void nrerror(const char error_text[]);
float *vector(int , int );
int *ivector(int , int );
double *dvector(int , int );
float **matrix(int , int , int , int );
double **dmatrix(int , int , int , int );
int **imatrix(int , int , int , int );
float **submatrix(float **, int , int , int , int , int , int );
float **convert_matrix(float *, int , int , int , int );
void free_vector(float *, int , int );
void free_ivector(int *, int , int );
void free_dvector(double *, int , int );
void free_matrix(float **, int , int , int , int );
void free_dmatrix(double **, int , int , int , int );
void free_imatrix(int **, int , int , int , int );
void free_submatrix(float **, int , int , int , int );
void free_convert_matrix(float **, int , int , int , int );
int siftDown(double *numbers, int root, int bottom, int *re_arrange);
int CheckMono(double *Vx, int ymax, int *mono_ind);
bool fileExists(const std::string& filename);
double four_dimension_linear_interpolation(double* lattice_spacing, double** fraction, double**** cube);
double three_dimension_linear_interpolation(double* lattice_spacing, double** fraction, double*** cube);

};

#endif
