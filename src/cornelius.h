#ifndef CORNELIUS_H
#define CORNELIUS_H

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>

using namespace std;

/**
 *
 * General class for elements. All other element classes are inherited
 * from this.
 *
 */
class GeneralElement
{
  protected:
    static const int DIM = 4;
    double *centroid;
    double *normal;
    int normal_calculated;
    int centroid_calculated;
    virtual void calculate_centroid() {};
    virtual void calculate_normal() {};
    void check_normal_direction(double *normal, double *out);
  public:
    GeneralElement();
    ~GeneralElement();
    double *get_centroid();
    double *get_normal();
};

/**
 *
 * This class is used for line elements gotten from the squares. Can
 * calculate the centroid and normal of the line.
 *
 */
class Line : public GeneralElement
{
  private:
    static const int LINE_CORNERS = 2;
    static const int LINE_DIM = 2;
    int x1,x2;
    int start_point;
    int end_point;
    double **corners;
    double *out;
    int *const_i;
    void calculate_centroid();
    void calculate_normal();
  public:
    Line();
    ~Line();
    void init(double**,double*,int*);
    void flip_start_end();
    double *get_start();
    double *get_end();
    double *get_out();
};

/**
 *
 * A class for storing polygons given by Cornelius. Can calculate
 * the centroid and normal of the polygon.
 *
 */
class Polygon : public GeneralElement
{
  private:
    static const int MAX_LINES = 24;
    static const int POLYGON_DIM = 3;
    Line **lines;
    int Nlines;
    int x1,x2,x3;
    int const_i;
    void calculate_centroid();
    void calculate_normal();
  public:
    Polygon();
    ~Polygon();
    void init(int);
    bool add_line(Line*,int);
    int get_Nlines();
    Line** get_lines();
    void print(ofstream &file,double*);
};

/**
 *
 * A class for storing polyhedrons given by Cornelius. Can calculate
 * the normal and centroid of the polyhedron.
 *
 */
class Polyhedron : public GeneralElement
{
  private:
    static const int MAX_POLYGONS = 24;
    Polygon **polygons;
    int Npolygons;
    int Ntetrahedra;
    int x1,x2,x3,x4;
    bool lines_equal(Line*,Line*);
    void tetravolume(double*,double*,double*,double*);
    void calculate_centroid();
    void calculate_normal();
  public:
    Polyhedron();
    ~Polyhedron();
    void init();
    bool add_polygon(Polygon*,int);
};

/**
 *
 * This class handles the squares which are splitted from the cube.
 * Finds the edges of the surface and also the poinst which point
 * the outward direction.
 *
 * 13.10.2011 Hannu Holopainen
 *
 */
class Square
{
  private:
    static const int DIM = 4;
    static const int SQUARE_DIM = 2;
    static const int MAX_POINTS = 4;
    static const int MAX_LINES = 2;
    double **points;
    double **cuts;
    double **out;
    double **points_temp;
    double *out_temp;
    int *const_i;
    double *const_value;
    int x1, x2;
    double *dx;
    int Ncuts;
    int Nlines;
    Line *lines;
    int ambiguous;
    void ends_of_edge(double);
    void find_outside(double);
  public:
    Square();
    ~Square();
    void init(double**,int*,double*,double*);
    void construct_lines(double);
    int is_ambiguous();
    int get_Nlines();
    Line* get_lines();
};

/**
 *
 * This class handles 3d-cubes. Splits them into squares and collects
 * polygons from the lines given from squares.
 *
 * 13.10.2011 Hannu Holopainen
 *
 */
class Cube
{
  private:
    static const int DIM = 4;
    static const int CUBE_DIM = 4;
    static const int MAX_POLY = 8;
    static const int NSQUARES = 6;
    static const int STEPS = 2;
    double ***cube;
    Line **lines;
    Polygon *polygons;
    Square *squares;
    int Nlines;
    int Npolygons;
    int ambiguous;
    int const_i;
    double const_value;
    int x1,x2,x3;
    double *dx;
    void split_to_squares();
    void check_ambiguous(int);
  public:
    Cube();
    ~Cube();
    void init(double***&,int,double,double*&);
    void construct_polygons(double);
    int get_Nlines();
    int get_Npolygons();
    int is_ambiguous();
    Polygon* get_polygons();
};


/**
 *
 * This class handles 4d-cubes. Splits them into squares and collects
 * polygons from the lines given from squares.
 *
 * 13.10.2011 Hannu Holopainen
 *
 */
class Hypercube
{
  private:
    static const int DIM = 4;
    static const int MAX_POLY = 10;
    static const int NCUBES = 8;
    static const int STEPS = 2;
    double ****hcube;
    Polyhedron *polyhedrons;
    Polygon **polygons;
    Cube *cubes;
    int Npolyhedrons;
    int ambiguous;
    int x1,x2,x3,x4;
    double *dx;
    void split_to_cubes();
    void check_ambiguous(double);
  public:
    Hypercube();
    ~Hypercube();
    void init(double****,double*);
    void construct_polyhedrons(double);
    int get_Npolyhedrons();
    Polyhedron* get_polyhedrons();
};

/**
 *
 * A class for finding a constant value surface from a 2-4 dimensional
 * cube. The values at corners and length of the side are needed as
 * an input.
 *
 * Algorithm by Pasi Huovinen. This code is based on the original FORTRAN
 * code by Pasi Huovinen.
 *
 * 23.04.2012 Hannu Holopainen
 *
 */
class Cornelius
{
  private:
    static const int STEPS = 2;
    static const int DIM = 4;
    static const int MAX_ELEMENTS = 10;
    int Nelements;
    double **normals;
    double **centroids;
    int cube_dim;
    int initialized;
    int print_initialized;
    double value0;
    double *dx;
    ofstream output_print;
    void surface_3d(double***,double*,int);
    Square cu2d;
    Cube cu3d;
    Hypercube cu4d;
  public:
    Cornelius();
    ~Cornelius();
    void init(int,double,double*);
    void init_print(string);
    void find_surface_2d(double**);
    void find_surface_3d(double***);
    void find_surface_3d_print(double***,double*);
    void find_surface_4d(double****);
    int get_Nelements();
    double **get_normals();
    double **get_centroids();
    double **get_normals_4d();
    double **get_centroids_4d();
    double get_centroid_elem(int,int);
    double get_normal_elem(int,int);
};

#endif /* CORNELIUS_H */
