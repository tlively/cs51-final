/********************************************************************
 * the header file for some handy vector math
 ********************************************************************/
#ifndef VECTOR_INCLUDED
#define VECTOR_INCLUDED

/* A generic struct for handing two dimensional vectors. Used in many
 * situations with different meanings. */
typedef struct po_vector {
  float x;
  float y;
} po_vector;

/* a 2x2 matrix */
typedef struct matrix {
  po_vector rows[2];
} matrix;

/* forms a line vector from two coordinate points */
po_vector vect_from_points(po_vector origin, po_vector destination);

/* gets the right hand normal vector */
po_vector vect_normal(po_vector vect);

/* calculate the dot product of two vectors */
float vect_dot_prod(po_vector a, po_vector b);

/* get the maginitude of a vector squared*/
float vect_mag_squared(po_vector vect);

/* the projection of vector a onto vector p */
po_vector vect_project(po_vector a, po_vector p);

/* get corresponding axis from two points; vertex1 should be the orgin */
po_vector vect_axis(po_vector vertex1, po_vector vertex2);

/* makes 2x2 a matrix of the form 
 * {a,b}
 * {c,d} */
matrix vect_create_matrix(float a, float b, float c, float d);

/* multiplies a vector with a 2x2 matrix
 * matrix[n] for n is an int, will grab the n row of the matrix. 
 * matrix[n].x corresponds to the first column in that row, 
 * matrix[n].y to the second column in that row */
po_vector vect_matrix_mult(po_vector vert, matrix a_matrix);

/* cross product of two vectors in 2D */
float vect_cross_prod(po_vector vect1, po_vector vect2);

#endif
