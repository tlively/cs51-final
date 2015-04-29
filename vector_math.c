/********************************************************************
 * a file of vector math functions
 *
 * more of Annaleah's shenanigans
 ********************************************************************/
#include "vector_math.h"

/* forms a line vector from two coordinate points 
 * origin is subtracted from destination */
po_vector vect_from_points(po_vector origin, po_vector destination) {
  po_vector line;
  line.x = destination.x - origin.x;
  line.y = destination.y - origin.y;
  return line;
}

/* gets the right hand normal vector */
po_vector vect_normal(po_vector vect) {
  float temp = vect.x;
  vect.x = vect.y;
  vect.y = -temp;
  return vect;
}

/* calculate the dot product of two vectors */
float vect_dot_prod(po_vector a, po_vector b){
  return a.x * b.x + a.y * b.y;
}

/* the projection of vector a onto vector p */
po_vector vect_project(po_vector a, po_vector p){
  float coeff = vect_dot_prod(a,p) / vect_dot_prod(a,a);
  po_vector proj;
  proj.x = coeff * a.x;
  proj.y = coeff * a.y;
  return proj;
}

/* get corresponding axis from two points; vertex1 should be the orgin */
po_vector vect_axis(po_vector vertex1, po_vector vertex2) {
  return vect_normal(vect_from_points(vertex1,vertex2));
}

po_vector vect_matrix_mult(po_vector vert, po_vector* matrix){
  float temp = vert.x;
  /* matrix[n] for n is an int, will grab the n row of the matrix. matrix.x
   * and matrix.y extracts the columns of the matrix*/
  vert.x = matrix[0].x * temp + matrix[0].y * vert.y;
  vert.y = matrix[1].x * temp + matrix[1].y * vert.y;
  return vert;
}
