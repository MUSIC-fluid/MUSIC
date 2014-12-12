/**
 *
 * cornelius++ version 1.3: Copyright 2012, Pasi Huovinen and Hannu Holopainen
 *
 * This subroutine is aimed to be used as a part of the fluid dynamical models
 * of the heavy-ion physics community. Permission to use it for any purpose
 * except for any commercial purpose is granted, provided that any publication
 * cites the paper describing the algorithm:
 *   P. Huovinen and H. Petersen, arXiv:1206.3371 [nucl-th]
 *
 * Permission to distribute this subroutine is granted, provided that no fee is
 * charged, and that this copyright and permission notice appear in all the
 * copies. Permission to modify this subroutine is granted provided that the
 * modified subroutine is made publicly available latest when any results obtained
 * using the modified subroutine are published, the modified subroutine is
 * distributed under terms similar to this notice, and the modified code carries
 * both the original copyright notice and notices stating that you modified it,
 * and a relevant date/year.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of FITNESS FOR A PARTICULAR PURPOSE. 
 *
 * Last update 03.08.2012 Hannu Holopainen
 *
 */

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
 * Constructor allocates memory for centroid and normal. 
 *
 */
GeneralElement::GeneralElement()
{
  normal_calculated = 0;
  centroid_calculated = 0;
  centroid = new double[DIM];
  normal = new double[DIM];
}

/**
 *
 * Destructor frees memory.
 *
 */
GeneralElement::~GeneralElement()
{
  delete[] centroid;
  delete[] normal;
}

/**
 *
 * Checks if the normal is pointing into right direction. If the normal
 * is pointing in the wrong direction, it is turned to the right direction. 
 *
 * @param [in/out] normal The normal vector
 * @param [in]     out    Vector pointing into outward direction.
 *
 */
void GeneralElement::check_normal_direction(double *normal, double *out)
{
  //We calculate the dot product, if less than 4 dimensions the elements,
  //which are not used, are just zero and do not effect this
  double dot_product = 0;
  for (int i=0; i < DIM; i++) {
    dot_product += out[i]*normal[i];
  }
  //If the dot product is negative, the normal is pointing in the
  //wrong direction and we have to flip it's direction
  if ( dot_product < 0 ) {
    for (int i=0; i < DIM; i++) {
      normal[i] = -normal[i];
    }
  }
}

/**
 *
 * Returns the centroid as a (double*) always in 4d.
 *
 * @return Table containing centroid in 4D
 *
 */
double *GeneralElement::get_centroid()
{
  if ( !centroid_calculated )
    calculate_centroid();
  return centroid;
}

/**
 *
 * Returns the normal as a (double*) always in 4d.
 *
 * @return Table containing normal in 4D
 *
 */
double *GeneralElement::get_normal()
{
  if ( !normal_calculated )
    calculate_normal();
  return normal;
}


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
 * Constructor allocates memory for end points, point which is always
 * outside and indices of the elements which are constant.
 *
 */
Line::Line()
{
  GeneralElement();
  corners = new double*[LINE_CORNERS];
  for (int i=0; i < LINE_CORNERS; i++) {
    corners[i] = new double[DIM];
  }
  out = new double[DIM];
  const_i = new int[DIM-LINE_DIM]; 
}

/**
 *
 * Frees memory allocated in constructor.
 *
 */
Line::~Line()
{
  for (int i=0; i < LINE_CORNERS; i++) {
    delete[] corners[i];
  }
  delete[] corners;
  delete[] out;
  delete[] const_i;
}

/**
 *
 * Initializes line element. Same element can be initialized several
 * times with this function.
 *
 * @param [in] p Table containing end points
 * @param [in] o Table containing point which is outside surface
 * @param [in] o Table containing constant indices
 *
 */
void Line::init(double **p, double *o, int *c)
{
  //We copy the values of the end points and the outside point
  for (int i=0; i < LINE_CORNERS; i++) {
    for (int j=0; j < DIM; j++) {
      corners[i][j] = p[i][j];
      if ( i==0 ) {
        out[j] = o[j];
      }
    }
  }
  //Let's set indices for start and end points, so we can conveniently
  //flip the start and end point without changing the actual values
  start_point = 0;
  end_point = 1;
  //Here we copy the indices which are constant at this line
  for (int i=0; i < DIM-LINE_DIM; i++) {
    const_i[i] = c[i];
  }
  //Here we fix the non-zero indices in such a way that x1 is
  //always smaller
  if ( c[0] == 0 ) {
    if ( c[1] == 1 ) {
      x1 = 2; x2 = 3;
    } else if ( c[1] == 2 ) {
      x1 = 1; x2 = 3;
    } else {
      x1 = 1; x2 = 2;
    }
  } else if ( c[0] == 1 ) {
    if ( c[1] == 2 ) {
      x1 = 0; x2 = 3;
    } else {
      x1 = 0; x2 = 2;
    }
  } else {
    x1 = 0; x2 = 1;
  }
  //We need to set these to zero so when this function is called
  //again we know that normal and centroid are not yet calculated
  //with new values
  normal_calculated = 0;
  centroid_calculated = 0;
}

/**
 *
 * Determines the centroid for a line. Thus this a very trivial
 * operation.
 *
 */
void Line::calculate_centroid()
{
  //Centroid is just the average of the points
  for (int i=0; i < DIM; i++) {
    centroid[i] = (corners[0][i] + corners[1][i])/2.0;
  }
  centroid_calculated = 1;
}

/**
 *
 * Calculates the normal vector for a line.
 *
 */
void Line::calculate_normal()
{
  //Centroid must be calculated before we can calculate normal
  if ( !centroid_calculated )
    calculate_centroid();
  //Normal is just (-dy,dx)
  normal[x1] = -(corners[1][x2] - corners[0][x2]);
  normal[x2] = corners[1][x1] - corners[0][x1];
  normal[const_i[0]] = 0;
  normal[const_i[1]] = 0;
  //Now we check if the normal is in the correct direction
  double *Vout = new double[DIM];
  for (int j=0; j < DIM; j++) {
    Vout[j] = out[j] - centroid[j];
  }
  check_normal_direction(normal,Vout);
  delete[] Vout;
  normal_calculated = 1;
}

/**
 *
 * Flips the start and end point of the line.
 *
 */
void Line::flip_start_end()
{
  int temp = start_point;
  start_point = end_point;
  end_point = temp;
}

/**
 *
 * Returns the start point in 4d as double*.
 *
 * @return Start point as a table
 *
 */
double *Line::get_start()
{
  return corners[start_point];
}

/**
 *
 * Returns the end point in 4d as double*.
 *
 * @return End point as a table
 *
 */
double *Line::get_end()
{
  return corners[end_point];
}

/**
 *
 * Returns the point, which is always outside, in 4d as double*.
 *
 * @return Point which is always outside as a table
 *
 */
double *Line::get_out()
{
  return out;
}

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
 * Allocates memory for the list of pointers to the lines which
 * form the polygon.
 *
 */
Polygon::Polygon()
{
  GeneralElement();
  lines = new Line*[MAX_LINES];
}

/**
 *
 * Frees memory allocated at the constructor.
 *
 */
Polygon::~Polygon()
{
  delete[] lines;
}

/**
 *
 * This initializes the polygon. Can be used several times.
 * 
 * @param [in] c Index which is constant in this polygon 
 *
 */
void Polygon::init(int c)
{
  //Let's copy the constant index
  const_i = c;
  //Here we fix the indices which are not constant
  x1 = -1;
  x2 = -1;
  x3 = -1;
  for (int i=0; i < DIM; i++) {
    if ( i != c ) {
      if ( x1 < 0 ) {
        x1 = i;
      } else if ( x2 < 0 ) {
        x2 = i;
      } else {
        x3 = i;
      }
    }
  }
  //We need to set these to zero so when this function is called
  //again we know that polygon has no lines added and normal and
  //centroid are not yet calculated with new values
  Nlines = 0;
  normal_calculated = 0;
  centroid_calculated = 0;
}

/**
 *
 * Adds a line to the polygon. This can also check if the line
 * is connected to the ones already in the polygon and the line
 * is added only when the line is connected. 
 * 
 * @param [in] l          Line which is tried to add to polygon 
 * @param [in] donotcheck 1 if no check is needed, 0 otherwise
 *
 * @return True if line was added, false if not added
 *
 */
bool Polygon::add_line(Line *l, int donotcheck)
{
  double eps = 1e-10;
  //If this is the first line or we do not want to check, line is added
  //automatically
  if ( donotcheck || Nlines == 0 ) {
    lines[Nlines] = l;
    Nlines++;
    return true; 
  } else {
    //Here we check if the new line is connected with the last line in
    //polygon. This since lines are ordered.
    double *p1 =  l->get_start();
    double *p2 =  l->get_end();
    double *p3 =  lines[Nlines-1]->get_end();
    double diff1 = 0;
    double diff2 = 0;
    for (int j=0; j < DIM; j++) {
      diff1 += fabs(p1[j]-p3[j]);
      diff2 += fabs(p2[j]-p3[j]);
    }
    //If start or end point is the same as the end point of previous line, line
    //is connected to polygon and it is added
    if ( diff1 < eps || diff2 < eps ) {
      //We always want that start point is connected to end so we flip the
      //start and end in the new line if needed
      if ( diff2 < eps )
        l->flip_start_end();
      lines[Nlines] = l;
      Nlines++;
      return true;
    }
    //If we are here, the line was not connected to the current polygon
    return false;
  }
}

/**
 *
 * This determines the centroid (center of mass) of the polygon.
 *
 */
void Polygon::calculate_centroid()
{
  //We need a vector for the mean of the corners.
  double *mean = new double[DIM];
  for (int i=0; i < DIM; i++ ) {
    mean[i] = 0;
  }
  //First we determine the mean of the corner points. Every point
  //is here twice since lines are not ordered
  for (int i=0; i < Nlines; i++ ) {
    double *p1 = lines[i]->get_start(); 
    double *p2 = lines[i]->get_end(); 
    for (int j=0; j < DIM; j++ ) {
      mean[j] += p1[j] + p2[j];
    }
  }
  for (int j=0; j < DIM; j++ ) {
    mean[j] = mean[j]/double(2.0*Nlines);
  }
  //If only 3 lines, i.e. 3 corners, the polygon is always on a plane
  //and centroid is the mean
  if ( Nlines == 3 ) {
    for (int i=0; i < DIM; i++) {
      centroid[i] = mean[i];
    }
    centroid_calculated = 1;
    delete[] mean;
    return;
  }
  //If more than 3 corners, calculation of the centroid is more
  //complicated
  //Here we from triangles from the lines and the mean point
  double *sum_up = new double[DIM]; //areas of the single triangles 
  double sum_down = 0; //area of all the triangles
  for (int i=0; i < DIM; i++) {
    sum_up[i] = 0;
  }
  //a and b are vectors which from the triangle
  double *a = new double[DIM];
  double *b = new double[DIM];
  //centroid of the triangle (this is always on a plane)
  double *cm_i = new double[DIM];
  for (int i=0; i < Nlines; i++) {
    double *p1 = lines[i]->get_start();
    double *p2 = lines[i]->get_end();
    for (int j=0; j < DIM; j++) {
      cm_i[j] = (p1[j] + p2[j] + mean[j])/3.0;
    } 
    //Triangle is defined by these two vectors
    for (int j=0; j < DIM; j++) {
      a[j] = p1[j] - mean[j];
      b[j] = p2[j] - mean[j];
    }
    //Area is calculated from cross product of these vectors
    double A_i = 0.5*sqrt( pow(a[x2]*b[x3]-a[x3]*b[x2],2.0) + pow(a[x1]*b[x3]-a[x3]*b[x1],2.0) + pow(a[x2]*b[x1]-a[x1]*b[x2],2.0) );
    //Then we store the area and update the total area
    for (int j=0; j < DIM; j++) {
      sum_up[j] += cm_i[j]*A_i;
    }
    sum_down += A_i;
  }
  //Now we can determine centroid as a weighted average of the centroids
  //of the triangles
  for (int i=0; i < DIM; i++) {
    centroid[i] = sum_up[i]/sum_down;
  }
  //centroid is now calculated and memory is freed
  centroid_calculated = 1;
  delete[] sum_up;
  delete[] mean;
  delete[] a;
  delete[] b;
  delete[] cm_i;
}

/**
 *
 * Determines the normal vector for polygon. Makes sure that
 * it points in the outward direction.
 *
 */
void Polygon::calculate_normal()
{
  //centroid must be calculated before normal can be determined
  if ( !centroid_calculated )
    calculate_centroid();
  //First we find the normal for each triangle formed from
  //one edge and centroid
  double **normals = new double*[Nlines]; //these are the normals
  double *Vout = new double[DIM]; //the point which is always outside
  for (int i=0; i < Nlines; i++) {
    normals[i] = new double[DIM];
    for (int j=0; j < DIM; j++) {
      normals[i][j] = 0;
    }
  }
  //Normal is defined by these two vectors
  double *a = new double[DIM];
  double *b = new double[DIM];
  //Loop over all triangles
  for (int i=0; i < Nlines; i++) {
    //First we calculate the vectors which form the triangle
    double *p1 = lines[i]->get_start();
    double *p2 = lines[i]->get_end();
    for (int j=0; j < DIM; j++) {
      a[j] = p1[j] - centroid[j];
      b[j] = p2[j] - centroid[j];
    }
    //Normal is calculated as a cross product of these vectors
    normals[i][x1] =  0.5*(a[x2]*b[x3]-a[x3]*b[x2]);
    normals[i][x2] = -0.5*(a[x1]*b[x3]-a[x3]*b[x1]);
    normals[i][x3] =  0.5*(a[x1]*b[x2]-a[x2]*b[x1]);
    normals[i][const_i] = 0;
    //Then we construct a vector which points out
    double *o = lines[i]->get_out();
    for (int j=0; j < DIM; j++) {
      Vout[j] = o[j] - centroid[j];
    }
    //then we check that normal is point in the correct direction
    check_normal_direction(normals[i],Vout);
  }
  //Finally the normal is a sum of the normals of the triangles
  for (int i=0; i < DIM; i++) {
    normal[i] = 0;
  }
  for (int i=0; i < Nlines; i++) {
    for (int j=0; j < DIM; j++) {
      normal[j] += normals[i][j];
    }
  }
  //Normal is now calculated and memory can be freed
  normal_calculated = 1;
  delete[] a;
  delete[] b;
  for (int i=0; i < Nlines; i++) {
    delete[] normals[i];
  }
  delete[] normals;
  delete[] Vout;
}

/**
 *
 * Return the number of lines in this polygon.
 *
 * @return Number of lines in this polygon
 *
 */
int Polygon::get_Nlines()
{
  return Nlines;
}

/**
 *
 * Returns an array of pointers of the lines in the polygon.
 *
 * @return Array of pointer to the lines in polygon
 *
 */
Line** Polygon::get_lines()
{
  return lines;
}

/**
 *
 * Prints the triangles formed from the polygon into a given file. Prints
 * the absolute points, so this file can be used to plot the surface.
 *
 * @param [in] file Output stream where the points are printed.
 * @param [in] file Absolute position of the cube point (0,0,0) 
 *
 */
void Polygon::print(ofstream &file, double *pos)
{
  for (int i=0; i < Nlines; i++) {
    double *p1 = lines[i]->get_start();
    double *p2 = lines[i]->get_end();
    file << pos[x1] + p1[x1] << " " << pos[x2] + p1[x2]  << " " << pos[x3] + p1[x3];
    file << " " << pos[x1] + p2[x1] << " " << pos[x2] + p2[x2]  << " " << pos[x3] + p2[x3];
    file << " " << pos[x1] + centroid[x1] << " " << pos[x2] + centroid[x2]  << " " << pos[x3] + centroid[x3] << endl;
  }
}

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
 * Allocates memory for a list pointers to the polygons which form
 * the polyhedron.
 *
 */
Polyhedron::Polyhedron()
{
  GeneralElement();
  polygons = new Polygon*[MAX_POLYGONS];
}

/**
 *
 * Frees memory allocated at the constructor.
 *
 */
Polyhedron::~Polyhedron()
{
  delete[] polygons;
}

/**
 *
 * This initializes the polygon. Can be used several times.
 * 
 */
void Polyhedron::init()
{
  //Here we fix the non-constant indices
  x1 = 0;
  x2 = 1;
  x3 = 2;
  x4 = 3;
  //We need to set these to zero so when this function is called
  //again we know that polyhedron has no polyfons added and normal
  //and centroid are not yet calculated with new values
  Npolygons = 0;
  Ntetrahedra = 0;
  normal_calculated = 0;
  centroid_calculated = 0;
}

/**
 *
 * Adds a polygon to the polyhedron. Also possible to check if the new
 * polygon is connected with previous polygons and only make the addition
 * if this is the case.
 * 
 * @param [in] p          Polygon which is tried to add to polyhedron
 * @param [in] donotcheck 1 if no check is needed, 0 otherwise
 *
 * @return True if polygon was added, false if not added
 *
 */
bool Polyhedron::add_polygon(Polygon *p, int donocheck)
{
  //If this is the first polygon or we want to add it in any case, we
  //just add it automatically
  if ( donocheck || Npolygons == 0 ) {
    polygons[Npolygons] = p;
    Npolygons++;
    Ntetrahedra += p->get_Nlines();
    return true; 
  } else {
    //Here we check if the new polygon is connected with the last polygon in
    //polyhedron
    for (int i=0; i < Npolygons; i++) {
      int Nlines1 = p->get_Nlines();
      int Nlines2 = polygons[i]->get_Nlines();
      Line **lines1 = p->get_lines();
      Line **lines2 = polygons[i]->get_lines();
      for (int j=0; j < Nlines1; j++) {
        for (int k=0; k < Nlines2; k++) {
          if ( lines_equal(lines1[j],lines2[k]) ) {
            polygons[Npolygons] = p;
            Npolygons++;
            Ntetrahedra += p->get_Nlines();
            return true; 
          }
        }
      }
    }
    //If we are here, the polygon was not connected to polyhedron
    return false;
  }
}

/**
 *
 * Checks if two lines are connected.
 * 
 * @param [in] l1 Line 1
 * @param [in] l2 Line 2
 *
 * @return true if lines are connected, false if not 
 *
 */
bool Polyhedron::lines_equal(Line *l1, Line *l2)
{
  double eps = 1e-10;
  double *p11 = l1->get_start();
  double *p12 = l1->get_end();
  double *p21 = l2->get_start();
  double diff1 = 0;
  double diff2 = 0;
  for (int i=0; i < DIM; i++) {
    diff1 += fabs(p11[i]-p21[i]);
    diff2 += fabs(p12[i]-p21[i]);
    if ( diff1 > eps && diff2 > eps )
      return false;
  }
  //If we are here, the lines are connected
  return true;
}

/**
 *
 * Calculates the normal vector for tetrahedra, which also tells
 * the volume of the tetrahedra.
 *
 * @param [in]  a Vector that spans the tetrahedra
 * @param [in]  b Vector that spans the tetrahedra
 * @param [in]  c Vector that spans the tetrahedra
 * @param [out] n Normal vector for tetrahedra
 *
 */
void Polyhedron::tetravolume(double *a, double *b, double *c, double *n)
{
  double bc01 = b[0]*c[1]-b[1]*c[0];
  double bc02 = b[0]*c[2]-b[2]*c[0];
  double bc03 = b[0]*c[3]-b[3]*c[0];
  double bc12 = b[1]*c[2]-b[2]*c[1];
  double bc13 = b[1]*c[3]-b[3]*c[1];
  double bc23 = b[2]*c[3]-b[3]*c[2];
  n[0] =  1.0/6.0*(a[1]*bc23 - a[2]*bc13 + a[3]*bc12);
  n[1] = -1.0/6.0*(a[0]*bc23 - a[2]*bc03 + a[3]*bc02);
  n[2] =  1.0/6.0*(a[0]*bc13 - a[1]*bc03 + a[3]*bc01);
  n[3] = -1.0/6.0*(a[0]*bc12 - a[1]*bc02 + a[2]*bc01);
  /*n[0] =  1.0/6.0*(a[1]*(b[2]*c[3]-b[3]*c[2]) - a[2]*(b[1]*c[3]-b[3]*c[1])
                  +a[3]*(b[1]*c[2]-b[2]*c[1]));
  n[1] = -1.0/6.0*(a[0]*(b[2]*c[3]-b[3]*c[2]) - a[2]*(b[0]*c[3]-b[3]*c[0])
                  +a[3]*(b[0]*c[2]-b[2]*c[0]));
  n[2] =  1.0/6.0*(a[0]*(b[1]*c[3]-b[3]*c[1]) - a[1]*(b[0]*c[3]-b[3]*c[0])
                  +a[3]*(b[0]*c[1]-b[1]*c[0]));
  n[3] = -1.0/6.0*(a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0])
                  +a[2]*(b[0]*c[1]-b[1]*c[0]));*/
}

/**
 *
 * This determines the centroid (center of mass) of the polyhedron.
 * 
 */
void Polyhedron::calculate_centroid()
{
  double *mean = new double[DIM];
  for (int i=0; i < DIM; i++ ) {
    mean[i] = 0;
  }
  //First we determine the mean of the corners
  for (int i=0; i < Npolygons; i++ ) {
    int Nlines = polygons[i]->get_Nlines();
    Line **lines = polygons[i]->get_lines();
    for (int j=0; j < Nlines; j++ ) {
      double *p1 = lines[j]->get_start();
      double *p2 = lines[j]->get_end();
      for (int k=0; k < DIM; k++ ) {
        mean[k] += p1[k] + p2[k];
      }
    }
  }
  for (int j=0; j < DIM; j++ ) {
    mean[j] = mean[j]/double(2.0*Ntetrahedra);
  }
  //Some memory allocated for temporary variables
  double *a = new double[DIM];
  double *b = new double[DIM];
  double *c = new double[DIM];
  double *n = new double[DIM];
  double *cm_i = new double[DIM];
  double *sum_up = new double[DIM];
  double sum_down = 0;
  for (int i=0; i < DIM; i++) {
    sum_up[i] = 0;
  }
  for (int i=0; i < Npolygons; i++ ) {
    int Nlines = polygons[i]->get_Nlines();
    Line **lines = polygons[i]->get_lines();
    double *cent = polygons[i]->get_centroid();
    for (int j=0; j < Nlines; j++) {
      double *p1 = lines[j]->get_start();
      double *p2 = lines[j]->get_end();
      //CM of the tetrahedra
      for (int k=0; k < DIM; k++) {
        cm_i[k] = (p1[k] + p2[k] + cent[k] + mean[k])/4.0;
      } 
      //Volume of a tetrahedron is determined next
      //Tetrahedron is defined by these three vectors
      for (int k=0; k < DIM; k++) {
        a[k] = p1[k] - mean[k];
        b[k] = p2[k] - mean[k];
        c[k] = cent[k] - mean[k];
      }
      //Then we calculate the volume from the normal n 
      tetravolume(a,b,c,n);
      double V_i = 0; 
      for (int k=0; k < DIM; k++) {
        V_i += n[k]*n[k];
      }
      V_i = sqrt(V_i);
      //Then we add contributions to sums
      for (int k=0; k < DIM; k++) {
        sum_up[k] += cm_i[k]*V_i;
      }
      sum_down += V_i;
    }
  }
  //Now the centroid is the volume weighted average of individual
  //tetrahedra
  for (int i=0; i < DIM; i++) {
    centroid[i] = sum_up[i]/sum_down;
  }
  //Centroid is now calculated and we can free memory
  centroid_calculated = 1;
  delete[] a;
  delete[] b;
  delete[] c;
  delete[] n;
  delete[] cm_i;
  delete[] sum_up;
  delete[] mean;
}

/**
 *
 * Determines the normal vector for polygon. Makes sure that
 * it points in the outward direction.
 *
 */
void Polyhedron::calculate_normal()
{
  //We need centroid in the following calculation and thus we
  //need to check that it is calculated
  if ( !centroid_calculated )
    calculate_centroid();
  //First we allocate memory for temporary variables
  double *Vout = new double[DIM];
  double *a = new double[DIM];
  double *b = new double[DIM];
  double *c = new double[DIM];
  double **normals = new double*[Ntetrahedra];
  for (int i=0; i < Ntetrahedra; i++) {
    normals[i] = new double[DIM];
  }
  int Ntetra = 0; //Index of tetrahedron we are dealing with
  for (int i=0; i < Npolygons; i++ ) {
    int Nlines = polygons[i]->get_Nlines();
    Line **lines = polygons[i]->get_lines();
    double *cent = polygons[i]->get_centroid();
    for (int j=0; j < Nlines; j++) {
      double *p1 = lines[j]->get_start();
      double *p2 = lines[j]->get_end();
      //Tetrahedron is defined by these three vectors
      for (int k=0; k < DIM; k++) {
        a[k] = p1[k] - centroid[k];
        b[k] = p2[k] - centroid[k];
        c[k] = cent[k] - centroid[k];
      }
      //Normal is calculated with the same function as volume
      tetravolume(a,b,c,normals[Ntetra]);
      //Then we determine the direction towards lower energy
      double *o = lines[j]->get_out();
      for (int k=0; k < DIM; k++) {
        Vout[k] = o[k] - centroid[k];
      }
      check_normal_direction(normals[Ntetra],Vout);
      Ntetra++;
    }
  }
  //Finally the element normal is a sum of the calculated normals
  for (int i=0; i < DIM; i++) {
    normal[i] = 0;
  }
  for (int i=0; i < Ntetrahedra; i++) {
    for (int j=0; j < DIM; j++) {
      normal[j] += normals[i][j];
    }
  }
  //Normal is now determined and we can free memory
  normal_calculated = 1;
  for (int i=0; i < Ntetrahedra; i++) {
    delete[] normals[i];
  }
  delete[] a;
  delete[] b;
  delete[] c;
  delete[] normals;
  delete[] Vout;
}

/**
 *
 * This class handles the squares. Finds the edges of the surface and also
 * the points which are always outside the surface so that we can determine
 * correct direction for normal vector.
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
 * Constructor allocates memory for the basic things needed in the
 * process.
 *
 */
Square::Square()
{
  cuts = new double*[MAX_POINTS];
  out = new double*[MAX_POINTS];
  points = new double*[SQUARE_DIM];
  const_i = new int[DIM-SQUARE_DIM];
  const_value = new double[DIM-SQUARE_DIM];
  lines = new Line[MAX_LINES];
  for (int i=0; i < SQUARE_DIM; i++) {
    points[i] = new double[SQUARE_DIM];
  }
  for (int i=0; i < MAX_POINTS; i++) {
    cuts[i] = new double[SQUARE_DIM];
    out[i] = new double[SQUARE_DIM];
  }
  int Npoints = 2;
  points_temp = new double*[SQUARE_DIM];
  out_temp = new double[DIM];
  for (int j=0; j < Npoints; j++) {
    points_temp[j] = new double[DIM];
  }
  ambiguous = 0;
}

/**
 *
 * Destructor frees memory allocated in the constructor.
 *
 */
Square::~Square()
{
  for (int i=0; i < DIM-SQUARE_DIM; i++) {
    delete[] points[i];
  }
  delete[] points;
  for (int i=0; i < MAX_POINTS; i++) {
    delete[] cuts[i];
    delete[] out[i];
  }
  delete[] cuts;
  delete[] out;
  delete[] const_i;
  delete[] const_value;
  delete[] lines;
  for (int i=0; i < SQUARE_DIM; i++) {
    delete[] points_temp[i];
  }
  delete[] points_temp;
  delete[] out_temp;
}

/**
 *
 * Initializes the square. Can be used several times to replace the old square
 * with a new one.
 *
 * @param [in] sq    Values at ends of the square.
 * @param [in] c_i   Indices that are constant on this square
 * @param [in] s_i   Indices that are x and y in this square
 * @param [in] c_val Values for the constant indices
 * @param [in] dex   Length of sides x and y
 *
 */
void Square::init(double **sq, int *c_i, double *c_v, double *dex)
{
  for (int i=0; i < SQUARE_DIM; i++) {
    for (int j=0; j < SQUARE_DIM; j++) {
      points[i][j] = sq[i][j];
    }
  }
  for (int i=0; i < DIM-SQUARE_DIM; i++) {
    const_i[i] = c_i[i];
    const_value[i] = c_v[i];
  }
  dx = dex;
  x1 = -1;
  x2 = -1;
  for (int i=0; i < DIM; i++) {
    if ( i != const_i[0] && i != const_i[1] ) {
      if ( x1 < 0 ) {
        x1 = i;
      } else {
        x2 = i;
      }
    }
  }
  //We need to set these to zero so when this function is called
  //again we know that nothing has been done for this square
  Ncuts = 0;
  Nlines = 0;
  ambiguous = 0;
}

/**
 *
 * Construct the lines which represent the surface on this square.
 *
 * @param [in] E0 The value which defines surface
 *
 */
void Square::construct_lines(double E0)
{
  //First we check the corner points and see if there are any lines
  int above=0;
  for (int i=0; i < DIM-SQUARE_DIM; i++) {
    for (int j=0; j < DIM-SQUARE_DIM; j++) {
      if ( points[i][j] >= E0 )
        above++; 
    }
  }
  //If all corners are above or below the surface value, no lines in this
  //square
  if ( above == 0 || above == 4 ) {
    Nlines = 0;
    return;
  }
  //First we find the cut points and the points which are always outside of the
  //surface. Also find_outside() arranges cuts so that first two cuts form a line
  //as defined in the algorithm (in case there are 4 cuts)
  ends_of_edge(E0);
  if ( Ncuts > 0 )
    find_outside(E0);
  //Then we go through the cut points and form the line elements 
  for (int i=0; i < Ncuts; i++) {
    points_temp[i%2][x1] = cuts[i][0];
    points_temp[i%2][x2] = cuts[i][1];
    points_temp[i%2][const_i[0]] = const_value[0];
    points_temp[i%2][const_i[1]] = const_value[1];
    //When we have inserted both endpoints we insert outside point
    //and we are ready to create line element.
    if ( i%2 == 1 ) {
      out_temp[x1] = out[i/2][0];
      out_temp[x2] = out[i/2][1];
      out_temp[const_i[0]] = const_value[0];
      out_temp[const_i[1]] = const_value[1];
      lines[i/2].init(points_temp,out_temp,const_i);
      Nlines++;
    }
  }
}

/**
 *
 * This finds the points from the edges where the surface crosses
 * the given square.
 *
 * The edges are gone through in the following order
 *      4
 *    -----
 *  1 |   | 3
 *    -----
 *      2
 * Since in the case where the ends are on each edge the ends are
 * assumed to be connected like \\
 *
 * @param [in] E0 Value that defines this surface
 *
 */
void Square::ends_of_edge(double E0)
{
  //Edge 1
  if ( (points[0][0]-E0)*(points[1][0]-E0) < 0 ) {
    cuts[Ncuts][0] = (points[0][0]-E0)/(points[0][0]-points[1][0])*dx[x1];
    cuts[Ncuts][1] = 0;
    Ncuts++;
  } else if ( points[0][0] == E0 && points[1][0] < E0 ) {
    cuts[Ncuts][0] = 1e-9*dx[x1];
    cuts[Ncuts][1] = 0;
    Ncuts++;
  } else if ( points[1][0] == E0 && points[0][0] < E0 ) {
    cuts[Ncuts][0] = (1.0-1e-9)*dx[x1];
    cuts[Ncuts][1] = 0;
    Ncuts++;
  }
  //Edge 2
  if ( (points[0][0]-E0)*(points[0][1]-E0) < 0 ) {
    cuts[Ncuts][0] = 0;
    cuts[Ncuts][1] = (points[0][0]-E0)/(points[0][0]-points[0][1])*dx[x2];
    Ncuts++;
  } else if ( points[0][0] == E0 && points[0][1] < E0 ) {
    cuts[Ncuts][0] = 0;
    cuts[Ncuts][1] = 1e-9*dx[x2];
    Ncuts++;
  } else if ( points[0][1] == E0 && points[0][0] < E0 ) {
    cuts[Ncuts][0] = 0;
    cuts[Ncuts][1] = (1.0-1e-9)*dx[x2];
    Ncuts++;
  }
  //Edge 3
  if ( (points[1][0]-E0)*(points[1][1]-E0) < 0 ) {
    cuts[Ncuts][0] = dx[x1];
    cuts[Ncuts][1] = (points[1][0]-E0)/(points[1][0]-points[1][1])*dx[x2];
    Ncuts++;
  } else if ( points[1][0] == E0 && points[1][1] < E0 ) {
    cuts[Ncuts][0] = dx[x1];
    cuts[Ncuts][1] = 1e-9*dx[x2];
    Ncuts++;
  } else if ( points[1][1] == E0 && points[1][0] < E0 ) {
    cuts[Ncuts][0] = dx[x1];
    cuts[Ncuts][1] = (1.0-1e-9)*dx[x2];
    Ncuts++;
  }
  //Edge 4
  if ( (points[0][1]-E0)*(points[1][1]-E0) < 0 ) {
    cuts[Ncuts][0] = (points[0][1]-E0)/(points[0][1]-points[1][1])*dx[x1];
    cuts[Ncuts][1] = dx[x2];
    Ncuts++;
  } else if ( points[0][1] == E0 && points[1][1] < E0 ) {
    cuts[Ncuts][0] = 1e-9*dx[x1];
    cuts[Ncuts][1] = dx[x2];
    Ncuts++;
  } else if ( points[1][1] == E0 && points[0][1] < E0 ) {
    cuts[Ncuts][0] = (1.0-1e-9)*dx[x1];
    cuts[Ncuts][1] = dx[x2];
    Ncuts++;
  }
  if ( Ncuts != 0 && Ncuts != 2 && Ncuts != 4 ) {
    cout << "Error in EndsOfEdge: Ncuts " << Ncuts << endl;
    exit(1);
  }
}

/**
 *
 * Finds the points showing the outside direction. 
 *
 * @param [in] E0 Value that defines this surface
 *
 */
void Square::find_outside(double E0)
{
  if ( Ncuts == 4 ) { //Here the surface is ambiguous
    ambiguous = 1;
    //Let's calculate the value in the middle of the square
    double Emid = 0;
    for (int i=0; i < 2; i++) {
      for (int j=0; j < 2; j++) {
        Emid += 0.25*points[i][j];
      }
    }
    //The default is that cuts are connected as \\ here
    //If both Emid and (0,0) are above or below the criterion
    //the cuts should be like // and we have to switch order in cuts
    if ( (points[0][0] < E0 && Emid < E0) || (points[0][0] > E0 && Emid > E0) ) {
      for (int i=0; i < 2; i++) {
        double temp = cuts[1][i];
        cuts[1][i] = cuts[2][i];
        cuts[2][i] = temp;
      }
    }
    //The center is below, so the middle point is always outside
    //the surface
    if ( (Emid-E0) < 0 ) {
      for (int i=0; i < 2; i++) {
        for (int j=0; j < 2; j++) {
          if ( j == 0 )
            out[i][j] = 0.5*dx[x1];
          else
            out[i][j] = 0.5*dx[x2];
        }
      }
    } else { // The center is above
       // Cuts are \\ here so bl and tr corners are outside
      if ( (points[0][0]-E0) < 0 ) {
        for (int i=0; i < 2; i++) {
          out[0][i] = 0;
          if ( i == 0 )
            out[1][i] = dx[x1];
          else
            out[1][i] = dx[x2];
        }
      // Cuts are // here so br and tl corners are outside
      } else {
        out[0][0] = dx[x1];
        out[0][1] = 0;
        out[1][0] = 0;
        out[1][1] = dx[x2];
      }
    }
  } else { //This is the normal case (not ambiguous)
    for (int i=0; i < 2; i++) {
      for (int j=0; j < 2; j++) {
        out[i][j] = 0;
      }
    }
    int Nout = 0;
    for (int i=0; i < 2; i++) {
      for (int j=0; j < 2; j++) {
        if ( points[i][j] < E0 ) {
          out[0][0] += i*dx[x1];
          out[0][1] += j*dx[x2];
          Nout++;
        }
      }
    }
    if ( Nout > 0 ) {
      for (int i=0; i < 2; i++) {
        for (int j=0; j < 2; j++) {
          out[i][j] = out[i][j]/double(Nout);
        }
      }
    }
  }
}

/**
 *
 * Tells if the case is ambiguous in this square.
 *
 * @return Zero if not ambiguous, one if ambiguous
 *
 */
int Square::is_ambiguous()
{
  return ambiguous;
}

/**
 *
 * Returns the number of lines in this square.
 * 
 * @return Number of lines in this square
 *
 */
int Square::get_Nlines()
{
  return Nlines;
}

/**
 *
 * Returns a table containing the lines found from this square.
 *
 * @return List of lines found from square
 *
 */
Line* Square::get_lines()
{
  return lines;
}


/**
 *
 * This class handles 3d-cubes. Splits them into squares and collects
 * polygons from the lines given from squares.
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
 * Constructor allocated memory.
 *
 */
Cube::Cube()
{
  cube = new double**[STEPS];
  for (int i=0; i < STEPS; i++) {
    cube[i] = new double*[STEPS];
    for (int j=0; j < STEPS; j++) {
      cube[i][j] = new double[STEPS];
    }
  }
  lines = new Line*[NSQUARES*2]; //Each square may have max. 2 lines
  polygons = new Polygon[MAX_POLY];
  squares = new Square[NSQUARES];
}

/**
 *
 * Destructor frees memory allocated in constructor.
 *
 */
Cube::~Cube()
{
  for (int i=0; i < STEPS; i++) {
    for (int j=0; j < STEPS; j++) {
      delete[] cube[i][j];
    }
    delete[] cube[i];
  }
  delete[] cube;
  delete[] lines;
  delete[] polygons;
  delete[] squares;
}

/**
 *
 * Initializes the cube. Can be used several times to replace the old cube
 * wit a new one.
 *
 * @param [in] c     Values at the corners of cube
 * @param [in] c_i   Index that is constant at this square
 * @param [in] c_val Value for the constant index
 * @param [in] dex   Lenghts of the sides
 *
 */
void Cube::init(double ***&c, int c_i, double c_v, double *&dex)
{
  const_i = c_i;
  const_value = c_v;
  dx = dex;
  //Here we fix the non-zero indices
  x1 = -1;
  x2 = -1;
  x3 = -1;
  for (int i=0; i < DIM; i++) {
    if ( i != c_i ) {
      if ( x1 < 0 ) {
        x1 = i;
      } else if ( x2 < 0 ) {
        x2 = i;
      } else {
        x3 = i;
      }
    }
  }
  for (int i=0; i < STEPS; i++) {
    for (int j=0; j < STEPS; j++) {
      for (int k=0; k < STEPS; k++) {
        cube[i][j][k] = c[i][j][k];
      }
    }
  }
  //We need to set these to zero so when this function is called
  //again we know that nothing has been done for this cube
  Nlines = 0;
  Npolygons = 0;
  ambiguous = 0;
}

/**
 *
 * Splits given cube to squares.
 *
 */
void Cube::split_to_squares()
{
  double **sq = new double*[STEPS];
  int *c_i = new int[STEPS];
  double *c_v = new double[STEPS];
  for (int i=0; i < STEPS; i++) {
    sq[i] = new double[STEPS]; 
  }
  int Nsquares = 0;
  for (int i=0; i < DIM; i++) {
    //i is the index which is kept constant, thus we ignore the index which
    //is constant in this cube
    if ( i != const_i ) {
      c_i[0] = const_i;
      c_i[1] = i;
      for (int j=0; j < STEPS; j++) {
        c_v[0] = const_value;
        c_v[1] = j*dx[i];
        for (int ci1=0; ci1 < STEPS; ci1++) {
          for (int ci2=0; ci2 < STEPS; ci2++) {
            if ( i == x1 ) {
              sq[ci1][ci2] = cube[j][ci1][ci2];
            } else if ( i == x2 ) {
              sq[ci1][ci2] = cube[ci1][j][ci2];
            } else {
              sq[ci1][ci2] = cube[ci1][ci2][j];
            }
          }
        }
        squares[Nsquares].init(sq,c_i,c_v,dx);
        Nsquares++;
      }
    }
  }
  for (int i=0; i < STEPS; i++) {
    delete[] sq[i];
  }
  delete[] sq; 
  delete[] c_i;
  delete[] c_v;
}

/**
 *
 * Here we construct polygons from the lines. If the surface cannot be
 * ambiguous, all lines are just added to single polygon, but if the
 * surface might be ambiguous we order the lines and determine how many
 * polygons we have.
 *
 */
void Cube::construct_polygons(double value0)
{
  //Let's start by splitting the cube to squares and finding
  //the lines from them
  split_to_squares();
  for (int i=0; i < NSQUARES; i++) {
    squares[i].construct_lines(value0);
  } 
  //Then we make a table which contains pointers to the lines
  Nlines = 0;
  for (int i=0; i < NSQUARES; i++) {
    int Nline = squares[i].get_Nlines();
    Line *l = squares[i].get_lines();
    for (int j=0; j < Nline; j++) {
      lines[Nlines] = &l[j];
      double *p = lines[Nlines]->get_start();
      p = lines[Nlines]->get_end();
      Nlines++; 
    }
    l = NULL;
  }
  //If no lines were found we may exit. This can happen only in 4D case
  if ( Nlines == 0 ) {
    return;
  }
  //Then we check if the surface is ambiguous and continue accordingly
  check_ambiguous(Nlines);
  if ( ambiguous > 0 ) {
    //Surface is ambiguous, so let's connect the lines to polygons and see how
    //many polygons we have
    int *not_used = new int[Nlines];
    for (int i=0; i < Nlines; i++) {
      not_used[i] = 1;
    }
    int used = 0; //this keeps track how many lines we have used
    //We may found several polygons with this
    do {
      if ( Nlines - used < 3 ) {
        cout << "Error: cannot construct polygon from " << Nlines -used << " lines" << endl;
        exit(1);
      }
      //First we initialize polygon
      polygons[Npolygons].init(const_i);
      //Then we go through all lines and try to add them to polygon
      for (int i=0; i < Nlines; i++) {
        if ( not_used[i] ) {
          //add_line returns true if line is succesfully added
          if ( polygons[Npolygons].add_line(lines[i],0) ) {
            not_used[i]=0;
            used++;
            //if line is succesfully added we start the loop from the
            //beginning
            i=0;
          }
        }
      }
      //When we have reached this point one complete polygon is formed
      Npolygons++;
    } while ( used < Nlines );
    delete[] not_used;
  } else {
    //Surface is not ambiguous, so we have only one polygons and all lines
    //can be added to it without ordering them
    polygons[Npolygons].init(const_i);
    for (int i=0; i < Nlines; i++) {
      polygons[Npolygons].add_line(lines[i],1);
    }
    Npolygons++;
  }
}

/**
 *
 * Checks if the surface is ambiguous in this cube.
 *
 */
void Cube::check_ambiguous(int Nlines)
{
  //First we check if any squares may have ambiguous elements
  for (int i=0; i < NSQUARES; i++) {
    if ( squares[i].is_ambiguous() ) {
      ambiguous++;
    }
  }
  //If the surface is not ambigous already, it is still possible to
  //have a ambiguous case if we have exatcly 6 lines, i.e. the surface
  //elements are at the opposite corners
  if ( ambiguous == 0 && Nlines == 6 ) {
    ambiguous++;
  }
}

/**
 *
 * Tells if the surface elements is ambiguous in this cube. Zero if not
 * ambiguous, non-zero if ambiguous.
 *
 * @return Zero if not ambiguous, non-zero if ambiguous
 *
 */
int Cube::is_ambiguous()
{
  return ambiguous;
}

/**
 *
 * Returns the number of the lines in this cube.
 *
 * @return Number of the lines
 *
 */
int Cube::get_Nlines()
{
  return Nlines;
}

/**
 *
 * Returns the number of the polygons in this cube.
 * 
 * @return Number of the polygons in this cube
 *
 */
int Cube::get_Npolygons()
{
  return Npolygons;
}

/**
 *
 * Returns a table containing the polygons found from this cube.
 *
 * @return List of the polygons found from this cube
 *
 */
Polygon* Cube::get_polygons()
{
  return polygons;
}

/**
 *
 * This class handles 4d-cubes. Splits them into cubes and collects the
 * polygons from the cubes and form polyhedrons from these. 
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
 * Constructor allocates memory.
 *
 */
Hypercube::Hypercube()
{
  hcube = new double***[STEPS];
  for (int i=0; i < STEPS; i++) {
    hcube[i] = new double**[STEPS];
    for (int j=0; j < STEPS; j++) {
      hcube[i][j] = new double*[STEPS];
      for (int k=0; k < STEPS; k++) {
        hcube[i][j][k] = new double[STEPS];
      }
    }
  }
  polyhedrons = new Polyhedron[MAX_POLY];
  polygons = new Polygon*[NCUBES*10];
  cubes = new Cube[NCUBES];
}

/**
 *
 * Destructor frees memory allocated in the constructor.
 *
 */
Hypercube::~Hypercube()
{
  for (int i=0; i < STEPS; i++) {
    for (int j=0; j < STEPS; j++) {
      for (int k=0; k < STEPS; k++) {
        delete[] hcube[i][j][k];
      }
      delete[] hcube[i][j];
    }
    delete[] hcube[i];
  }
  delete[] hcube;
  delete[] polyhedrons;
  delete[] polygons;
  delete[] cubes;
}

/**
 *
 * Initialized the hypercube. Can be used several times to replace the old
 * hypercube with a new one
 *
 * @param [in] c     Values at the corners of cube
 * @param [in] dex   Lenghts of the sides
 *
 */
void Hypercube::init(double ****c, double *dex)
{
  dx = dex;
  //Here we fix the non-zero indices
  x1 = 0;
  x2 = 1;
  x3 = 2;
  x4 = 3;
  for (int i=0; i < STEPS; i++) {
    for (int j=0; j < STEPS; j++) {
      for (int k=0; k < STEPS; k++) {
        for (int l=0; l < STEPS; l++) {
          hcube[i][j][k][l] = c[i][j][k][l];
        }
      }
    }
  }
  //We need to set these to zero so we can use this function to initialize
  //many hypercubes
  Npolyhedrons = 0;
  ambiguous = 0;
}

/**
 *
 * Splits given hypercube to cubes.
 *
 */
void Hypercube::split_to_cubes()
{
  double ***cu = new double**[STEPS];
  for (int i=0; i < STEPS; i++) {
    cu[i] = new double*[STEPS]; 
    for (int j=0; j < STEPS; j++) {
      cu[i][j] = new double[STEPS]; 
    }
  }
  int Ncubes = 0;
  for (int i=0; i < DIM; i++) {
    //i is the index which is kept constant, thus we ignore the index which
    //is constant in this cube
    int c_i = i;
    for (int j=0; j < STEPS; j++) {
      double c_v = j*dx[i];
      for (int ci1=0; ci1 < STEPS; ci1++) {
        for (int ci2=0; ci2 < STEPS; ci2++) {
          for (int ci3=0; ci3 < STEPS; ci3++) {
            if ( i == x1 ) {
              cu[ci1][ci2][ci3] = hcube[j][ci1][ci2][ci3];
            } else if ( i == x2 ) {
              cu[ci1][ci2][ci3] = hcube[ci1][j][ci2][ci3];
            } else if ( i == x3 ) {
              cu[ci1][ci2][ci3] = hcube[ci1][ci2][j][ci3];
            } else {
              cu[ci1][ci2][ci3] = hcube[ci1][ci2][ci3][j];
            }
          }
        }
      }
      cubes[Ncubes].init(cu,c_i,c_v,dx);
      Ncubes++;
    }
  }
  for (int i=0; i < STEPS; i++) {
    for (int j=0; j < STEPS; j++) {
      delete[] cu[i][j];
    }
    delete[] cu[i];
  }
  delete[] cu;
}

/**
 *
 * Here we construct polyhedrons from the polygons. First we check if the surface
 * is ambiguous and if it is not amibiguous we add all polygons to a single
 * polyhedron. If surface is ambiguous we need to connect the polygons one by one
 * and see in the end how many polyhedrons we have. 
 *
 */
void Hypercube::construct_polyhedrons(double value0)
{
  split_to_cubes();
  for (int i=0; i < NCUBES; i++) {
    cubes[i].construct_polygons(value0);
  }
  int Npolygons = 0;
  //then polygons are loaded to table
  for (int i=0; i < NCUBES; i++) {
    int Npoly = cubes[i].get_Npolygons();
    Polygon *p = cubes[i].get_polygons();
    for (int j=0; j < Npoly; j++) {
      polygons[Npolygons] = &p[j];
      Npolygons++; 
    }
  }
  check_ambiguous(value0);
  if ( ambiguous > 0 ) {
    //Here surface might be ambiguous and we need to connect the polygons and
    //see how many polyhedrons we have
    int *not_used = new int[Npolygons];
    for (int i=0; i < Npolygons; i++) {
      not_used[i] = 1;
    }
    int used = 0; //this keeps track how many lines we have used
    //We may found several polyhedrons with this
    do {
      //First we initialize polyhedron
      polyhedrons[Npolyhedrons].init();
      //Then we go through all polygons and try to add them to polyhedron
      for (int i=0; i < Npolygons; i++) {
        if ( not_used[i] ) {
          //add_polygon returns true if the polygon is succesfully added
          if ( polyhedrons[Npolyhedrons].add_polygon(polygons[i],0) ) {
            not_used[i]=0;
            used++;
            //if polygon is succesfully added we start the loop from the
            //beginning
            i=0;
          }
        }
      }
      //When we have reached this point one complete polyhedron is formed
      Npolyhedrons++;
    } while ( used < Npolygons );
    delete[] not_used;
    /*if ( ambiguous == 0 && Npolyhedrons != 1 ) {
      cout << "error" << endl;
    }*/
  } else {
    //Here surface cannot be ambiguous and all polygons can be added to
    //the polygehdron without ordering them
    polyhedrons[Npolyhedrons].init();
    for (int i=0; i < Npolygons; i++) {
      polyhedrons[Npolyhedrons].add_polygon(polygons[i],1);
    }
    //When we have reached this point one complete polyhedron is formed
    Npolyhedrons++;
  }
}

/**
 *
 * Checks if the surface is ambiguous in this hypercube.
 *
 */
void Hypercube::check_ambiguous(double value0)
{
  for (int i=0; i < NCUBES; i++) {
    if ( cubes[i].is_ambiguous() ) {
      ambiguous++;
    }
  }
  if ( ambiguous == 0 ) {
    int Nlines = 0;
    for (int i=0; i < NCUBES; i++) {
      Nlines += cubes[i].get_Nlines();
    }
    int points = 0;
    for (int i1=0; i1 < STEPS; i1++) {
      for (int i2=0; i2 < STEPS; i2++) {
        for (int i3=0; i3 < STEPS; i3++) {
          for (int i4=0; i4 < STEPS; i4++) {
            if ( hcube[i1][i2][i3][i4] < value0 ) {
              points++;
            }
          }
        }
      }
    }
    if ( points > 8 ) {
      points = 16-points;
    }
    /* these are not needed
    if ( Nlines == 42 ) {
      //ambiguous++;
    } else if ( Nlines == 36 && ( points == 5 || points == 4 ) ) {
      //ambiguous++;
    } else if ( Nlines == 30 && points == 3 ) {
      //ambiguous++;
    }*/
    if ( Nlines == 24 && points == 2 ) {
      ambiguous++;
    }
  }
}

/**
 *
 * Returns the number of polyhedrons in this cube.
 * 
 * @return Number of polyhedrons in this hypercube
 *
 */
int Hypercube::get_Npolyhedrons()
{
  return Npolyhedrons;
}

/**
 *
 * Returns a table containing the polyhedrons found from this square.
 *
 * @return List of polyhedrons found from hypercube
 *
 */
Polyhedron* Hypercube::get_polyhedrons()
{
  return polyhedrons;
}

/**
 *
 * A class for finding a constant value surface from a 2-4 dimensional
 * cube. The values at corners and length of the side are needed as
 * an input.
 *
 * Algorithm by Pasi Huovinen. This code is based on the original FORTRAN
 * code by Pasi Huovinen.
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

/**
 *
 * Constructor allocates some memory.
 *
 */
Cornelius::Cornelius()
{
  Nelements = 0;
  initialized = 0;
  print_initialized = 0;
  normals = new double*[MAX_ELEMENTS];
  centroids = new double*[MAX_ELEMENTS];
  for (int i=0; i < MAX_ELEMENTS; i++) {
    normals[i] = new double[DIM];
    centroids[i] = new double[DIM];
  }
}

/**
 *
 * Destructor frees memory and closes printing file is necessary.
 *
 */
Cornelius::~Cornelius()
{
  for (int i=0; i < MAX_ELEMENTS; i++) {
    delete[] normals[i];
    delete[] centroids[i];
  }
  delete[] normals;
  delete[] centroids;
  if ( initialized ) {
    delete[] dx;
  }
  //If file for printing was opened we close it here.
  if ( print_initialized ) {
    output_print.close();
  }
}

/**
 *
 * Initializes Cornelius. Can be used several times without any problems. This
 * might be useful if cube size varies during evolution.
 *
 * @param [in] d   Dimension of the problem.
 * @param [in] v0  Value which defines the surface.
 * @param [in] dex Length of the sides of the cube. Must contains as many elements as
 *                 the dimension of the problem (dx1,dx2,...).
 *
 */
void Cornelius::init(int d, double v0, double *dex)
{
  cube_dim = d;
  value0 = v0;
  if ( initialized != 1 ) {
    dx = new double[DIM];
  }
  for (int i=0; i < DIM; i++) {
    if ( i < DIM-cube_dim ) {
      dx[i] = 1;
    } else {
      dx[i] = dex[i-(DIM-cube_dim)];
    }
  }
  initialized = 1;
}

/**
 *
 * This initializes the printing of the surface elements into file. Final elements
 * are not printed with this, but instead this is used to print the triangles,
 * which are used to construct the final elements.
 *
 * @param [in] filename Filename of the file where the information is printed.
 *
 */
void Cornelius::init_print(string filename)
{
  output_print.open(filename.c_str());
  print_initialized = 1;
}

/**
 *
 * Finds the surface elements in 2-dimensional case.
 *
 * @param [in] cube Values at the corners of the cube as a 2d table so that value
 *                  [0][0] is at (0,0,0) and [1][1] is at (dx1,dx2).
 *
 */
void Cornelius::find_surface_2d(double **cube)
{
  if ( !initialized || cube_dim != 2 ) {
    cout << "Cornelius not initialized for 2D case" << endl;
    exit(1);
  }
  int *c_i = new int[cube_dim];
  double *c_v = new double[cube_dim];
  c_i[0] = 0;
  c_i[1] = 1;
  c_v[0] = 0;
  c_v[1] = 0;
  cu2d.init(cube,c_i,c_v,dx);
  cu2d.construct_lines(value0);
  Nelements = cu2d.get_Nlines();
  Line *l = cu2d.get_lines();
  for (int i=0; i < Nelements; i++) {
    for (int j=0; j < DIM; j++) {
      normals[i][j] = l[i].get_normal()[j];
      centroids[i][j] = l[i].get_centroid()[j];
    }
  }
  delete[] c_i;
  delete[] c_v;
}

/**
 *
 * Finds the surface elements in 3-dimensional case.
 *
 * @param [in] cube Values at the corners of the cube as a 3d table so that value
 *                  [0][0][0] is at (0,0,0) and [1][1][1] is at (dx1,dx2,dx3).
 *
 */
void Cornelius::find_surface_3d(double ***cube)
{
  double *pos = NULL;
  surface_3d(cube,pos,0);
}

/**
 *
 * Finds the surface elements in 3-dimensional case and prints the actual triangles
 * which are found by the algorithm.
 *
 * @param [in] cube Values at the corners of the cube as a 3d table so that value
 *                  [0][0][0] is at (0,0,0) and [1][1][1] is at (dx1,dx2,dx3).
 * @param [in] pos  Absolute position at the point [0][0][0] in form (0,x1,x2,x3).
 *
 */
void Cornelius::find_surface_3d_print(double ***cube, double *pos)
{
  surface_3d(cube,pos,1);
}

/**
 *
 * Finds the surface elements in 3-dimensional case.
 *
 * @param [in] cube     Values at the corners of the cube as a 3d table so that value
 *                      [0][0][0] is at (0,0,0) and [1][1][1] is at (dx1,dx2,dx3).
 * @param [in] pos      Absolute position at the point [0][0][0] in form (0,x1,x2,x3).
 * @param [in] do_print 1 if triangles are printed, otherwise 0
 *
 */
void Cornelius::surface_3d(double ***cube, double *pos, int do_print)
{
  if ( !initialized || cube_dim != 3 ) {
    cout << "Cornelius not initialized for 3D case" << endl;
    exit(1);
  }
  //First we check if the cube actually contains surface elements.
  //If all or none of the elements are below the criterion, no surface
  //elements exist.
  int above = 0;
  for (int i=0; i < STEPS; i++) {
    for (int j=0; j < STEPS; j++) { 
      for (int k=0; k < STEPS; k++) {
        if ( cube[i][j][k] >= value0 )
          above++;
      }
    }
  }
  if ( above == 0 || above == 8 ) {
    //No elements in this cube
    Nelements = 0;
    return;
  }
  //This cube has surface elements, so let's start by constructing
  //the cube.
  int c_i = 0;
  double c_v = 0;
  cu3d.init(cube,c_i,c_v,dx);
  //Then we find the elements
  cu3d.construct_polygons(value0);
  //Now let's get the information about the elements
  Nelements = cu3d.get_Npolygons();
  Polygon *p = cu3d.get_polygons();
  for (int i=0; i < Nelements; i++) {
    //Here we always work with 4-dimensions
    for (int j=0; j < DIM; j++) {
      normals[i][j] = p[i].get_normal()[j];
      centroids[i][j] = p[i].get_centroid()[j];
    }
    //If we want to print the actual triangles, which are found, we must do it here
    //before polygons are removed from the memory.
    if ( print_initialized && do_print ) {
      p[i].print(output_print,pos);
    }
  }
}

/**
 *
 * Finds the surface elements in 4-dimensional case.
 *
 * @param [in] cube Values at the corners of the cube as a 4d table so that value
 *                  [0][0][0][0] is at (0,0,0,0) and [1][1][1][1] is at
 *                  (dx1,dx2,dx3,dx4).
 *
 */
void Cornelius::find_surface_4d(double ****cube)
{
  if ( !initialized || cube_dim != 4 ) {
    cout << "Cornelius not initialized for 4D case" << endl;
    exit(1);
  }
  //First we check if the cube actually contains surface elements.
  //If all or none of the elements are below the criterion, no surface
  //elements exist.
  int above = 0;
  for (int i=0; i < STEPS; i++) {
    for (int j=0; j < STEPS; j++) { 
      for (int k=0; k < STEPS; k++) {
        for (int l=0; l < STEPS; l++) {
          if ( cube[i][j][k][l] >= value0 )
            above++;
        }
      }
    }
  }
  if ( above == 0 || above == 16 ) {
    //No elements in this cube
    Nelements = 0;
    return;
  }
  //This cube has surface elements, so let's start by constructing
  //the hypercube.
  cu4d.init(cube,dx);
  //Then we find the elements
  cu4d.construct_polyhedrons(value0);
  //Now let's get the information about the elements
  Nelements = cu4d.get_Npolyhedrons();
  Polyhedron *p = cu4d.get_polyhedrons();
  for (int i=0; i < Nelements; i++) {
    for (int j=0; j < DIM; j++) {
      centroids[i][j] = p[i].get_centroid()[j];
      normals[i][j] = p[i].get_normal()[j];
    }
  }
}

/**
 *
 * Returns the number of the surface elements in the given cube.
 *
 * @return Number of surface elements in the given cube.
 *
 */
int Cornelius::get_Nelements()
{
  return Nelements;
}

/**
 *
 * Returns the centroid vectors as a 2d table with the following number of indices
 * [number of elements][4]. If the dimension of the problem is smaller than four,
 * first (4-dimension) elements are zero. Note that this function allocates memory and
 * the user must free it!
 *
 * @return Table with dimensions [number of elements][4] containing the normal
 *         vectors of the surface elements.
 *
 */
double** Cornelius::get_centroids_4d()
{
  double **vect = new double*[Nelements];
  for (int i=0; i < Nelements; i++) {
    vect[i] = new double[DIM];
    for (int j=0; j < DIM; j++) {
      vect[i][j] = centroids[i][j];
    }
  }
  return vect;
}

/**
 *
 * Returns the normal vectors as a 2d table with the following number of indices
 * [number of elements][4]. If the dimension of the problem is smaller than four,
 * first (4-dimension) elements are zero. This gives \sigma_\mu without
 * factors(sqrt(-g)) from the metric. Note that this function allocates memory and
 * the user must free it!
 *
 * @return Table with dimensions [number of elements][4] containing the normal
 *         vectors of the surface elements.
 *
 */
double** Cornelius::get_normals_4d()
{
  double **vect = new double*[Nelements];
  for (int i=0; i < Nelements; i++) {
    vect[i] = new double[DIM];
    for (int j=0; j < DIM; j++) {
      vect[i][j] = normals[i][j];
    }
  }
  return vect;
}

/**
 *
 * Returns the centroid vectors as a 2d table with the following number of indices
 * [number of elements][dimension of the problem]. Note that this function allocates
 * memory and the user must free it!
 *
 * @return Table with dimensions [number of elements][dimensions] containing the centroid
 *         vectors of the surface elements.
 *
 */
double** Cornelius::get_centroids()
{
  double **vect = new double*[Nelements];
  for (int i=0; i < Nelements; i++) {
    vect[i] = new double[cube_dim];
    for (int j=0; j < cube_dim; j++) {
      vect[i][j] = centroids[i][j+(DIM-cube_dim)];
    }
  }
  return vect;
}

/**
 *
 * Returns the normal vectors as a 2d table with the following number of indices
 * [number of elements][dimension of the problem]. This gives \sigma_\mu without
 * factors(sqrt(-g)) from the metric. Note that this function allocates memory and
 * the user must free it!
 *
 * @return Table with dimensions [number of elements][dimensions] containing the normal
 *         vectors of the surface elements.
 *
 */
double** Cornelius::get_normals()
{
  double **vect = new double*[Nelements];
  for (int i=0; i < Nelements; i++) {
    vect[i] = new double[cube_dim];
    for (int j=0; j < cube_dim; j++) {
      vect[i][j] = normals[i][j+(DIM-cube_dim)];
    }
  }
  return vect;
}

/**
 *
 * Returns an element of the centroid vector.
 *
 * @param [in] i Number of surface element whose centroid one wants to get. Valid
 *               values are [0,number of elements in this cube].
 * @param [in] j Index of the element of centroid. Valid values are
 *               [0,dimension of the problem].
 * @return       Element j of the centroid of the surface element i.
 *
 */
double Cornelius::get_centroid_elem(int i, int j)
{
  if ( i >= Nelements || j >= cube_dim ) {
    cout << "Cornelius error: asking for an element which does not exist." << endl;
    exit(1);
  }
  return centroids[i][j+(DIM-cube_dim)];
}

/**
 *
 * Returns an element of the normal vector. This gives \sigma_\mu without factors
 * (sqrt(-g)) from the metric.
 *
 * @param [in] i Number of surface element whose normal one wants to get. Valid
 *               values are [0,number of elements in this cube].
 * @param [in] j Index of the element of centroid. Valid values are
 *               [0,dimension of the problem].
 * @return       Element j of the normal of the surface element i.
 *
 */
double Cornelius::get_normal_elem(int i, int j)
{
  if ( i >= Nelements || j >= cube_dim ) {
    cout << "Cornelius error: asking for an element which does not exist." << endl;
    exit(1);
  }
  return normals[i][j+(DIM-cube_dim)];
}
