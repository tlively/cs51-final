/********************************************************************
 * a file of vector math functions
 *
 * more of Annaleah's shenanigans
 ********************************************************************/
#include "vector_math.h"

/*vector subtraction */
po_vector vect_minus(po_vector vect1, po_vector vect2)
{
  po_vector res;
  res.x = vect1.x - vect2.x;
  res.y = vect1.y - vect2.y;
  return res;
}
/* forms a line vector from two coordinate points 
 * origin is subtracted from destination */
po_vector vect_from_points(po_vector origin, po_vector destination) {
  po_vector line;
  line.x = destination.x - origin.x;
  line.y = destination.y - origin.y;
  return line;
}

/* gets the slope from two points */
float vect_slope(po_vector p1, po_vector p2){
  po_vector slope = vect_from_points(p1,p2);
  return slope.y / slope.x;
}

/* gets the right hand normalized normal vector */
po_vector vect_normal(po_vector vect) {
  float temp = vect.x;
  vect.x = vect.y;
  vect.y = -temp;
  return vect;
}

/* turn this into a vector of length 1 */
po_vector vect_unit(po_vector vect){
  float mag = sqrt(vect_mag_squared(vect));
  vect.x = vect.x/mag;
  vect.y = vect.y/mag;
  return vect;
}

/* calculate the dot product of two vectors */
float vect_dot_prod(po_vector a, po_vector b){
  return a.x * b.x + a.y * b.y;
}

/* get the magnitude of a vector squared*/
float vect_mag_squared(po_vector vect) {
  return vect_dot_prod(vect, vect);
}

/* the projection of vector a onto vector p */
po_vector vect_project(po_vector a, po_vector p){
  float coeff = vect_dot_prod(a,p) / vect_mag_squared(a);
  po_vector proj;
  proj.x = coeff * a.x;
  proj.y = coeff * a.y;
  return proj;
}

/* get corresponding axis from two points; vertex1 should be the orgin */
po_vector vect_axis(po_vector vertex1, po_vector vertex2) {
  return vect_normal(vect_from_points(vertex1,vertex2));
}

/* creates a matrix of the form
 * {a,b}
 * {c,d} */
matrix vect_create_matrix (float a, float b, float c, float d){
  matrix a_matrix;
  a_matrix.rows[0].x = a;
  a_matrix.rows[0].y = b;
  a_matrix.rows[1].x = c;
  a_matrix.rows[1].y = d;
  return a_matrix;
}

/* multiplies a vector with a 2x2 matrix
 * matrix[n] for n is an int, will grab the n row of the matrix. 
 * matrix[n].x corresponds to the first column in that row, 
 * matrix[n].y to the second column in that row */
po_vector vect_matrix_mult(po_vector vert, matrix a_matrix){
  float temp = vert.x;
  vert.x = a_matrix.rows[0].x * temp + a_matrix.rows[0].y * vert.y;
  vert.y = a_matrix.rows[1].x * temp + a_matrix.rows[1].y * vert.y;
  return vert;
}

/* compute the 2d cros prod of two vectors; returns a number */
float vect_cross_prod(po_vector vect1, po_vector vect2){
  return vect1.x * vect2.y - vect2.x * vect1.y;
}
