/**************************************************************
 * Physics.c - the actual game engine
 *
 * Creates/removed physics objects
 * Updates physics world  
 **************************************************************/
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include "physics.h"
#include "dynamic_array.h"
#include "collisions.h"

int MAX_HASH_LEN = 1021;

// number of pixels per bucket in the spatial hash
#define BUCKET_SIZE 500;
#define INIT_SIZE 10;

/* the actual implementation of a physics object structure */
typedef struct po_imp {
  // change in x, y, and rotations
  float dx;
  float dy;
  float dr;

  // location (origin (x,y) in global coordinates)
  float x;
  float y;

  // positive r is counterclockwise
  float r;

  // the actual shape
  po_geometry shape;
  
  // used for collision detection and resolution; 0 if extrema have been defined
  int extrema_set;

  // centroid of the object in local coordinates
  po_vector centroid;

  // farthest distance from the centroid of the object
  float max_delta;
  
  // allows for linked lists within the hash table
  po_handle next;
} po_imp;

/* dat spatial hash though */
typedef struct world_t {
  // a dynamic array of dynamic arrays (hashed with the y vals)
  // this is basically a column pointing to rows
  dynamic_array* rows;
} world_t;

/* create a new world 
 * returns world on success, NULL on failure */
world_handle new_world () {
  // make new world
  world_handle world;
  world->rows = dynamic_array_create();
  return world;
}

/* add object to the physics world */
po_handle add_object (world_handle world, po_geometry* geom, 
		      float x, float y, float r) {
  // make a new object
  po_handle new_obj;
  new_obj->x = x;
  new_obj->y = y;
  new_obj->r = r;
  new_obj->dx = 0;
  new_obj->dy = 0;
  new_obj->dr = 0;
  new_obj->extrema_set = 1;
  new_obj->shape = *geom;
  
  // variables to store our x and y index
  int kx = x/BUCKET_SIZE;
  int ky = y/BUCKET_SIZE;
  
  // get array at that row number and figure out what's there
  dynamic_array* row_k = dynamic_array_get(world->rows,ky);

  if (row_k == NULL) {
    // update next pointer
    new_obj->next = NULL;
    
    // make a row here and place object at the xth place in that row
    row_k = dynamic_array_create();
    dynamic_array_add(row_k, kx, new_obj);

    // add this row to the yth index of rows
    dynamic_array_add(world->rows, ky, row_k);
  }
  else {
    // get whatever's at the kxth column of the row
    new_obj->next = dynamic_array_get(row_k, kx);

    // add updated object to the row at index x
    dynamic_array_add(row_k, kx, new_obj);
  }
}

int remove_object (world_handle world, po_handle obj){
  if (world == NULL || obj == NULL) {
    return 1;
  } 
  // get indexes
  int kx = obj->x/BUCKET_SIZE;
  int ky = obj->y/BUCKET_SIZE;

  // get row in the world index, remove stuff from that index from array
  dynamic_array* row_k = dynamic_array_get(world->rows,ky);
  po_handle obj_list = dynamic_array_remove(row_k, kx);

  // run through linked list to get the object we want to remove
  if (obj_list == obj) {
    // first node equality is a corner case; update the array
    dynamic_array_add(row_k, kx, obj_list->next);
    return 0;
  }
  else if (obj_list->next != NULL) {
    // the initializers for our for loop
    po_handle prev = obj_list;
    po_handle current = obj_list->next;
    
    // loop through the linked list
    while (current != NULL) {
      if (current == obj) {
	// update pointers, free current, and re-add the list sans our object
	prev->next = current->next;
	free(current);
	dynamic_array_add(row_k, kx, obj_list);
	return 0;
      }
      // update variables
      prev = current;
      current = current->next;
    }
  }
  // if we get here, we found nothing
  return 1;
}
/* Updates object's global position based on velocity
 * Future versions may include more sophistocated algorthims using acceleration
 * error checking should happen prior to passing things in */
void integrate (po_handle obj, float dx, float dy, float dr, float time_step) {
  // apply euler's method (the most logical choice since we have no accel)
  obj->x = obj->x + (dx * time_step);
  obj->y = obj->y + (dy * time_step);
  obj->r = obj->r + (dr * time_step);
}

int set_location (po_handle obj, float x, float y) {
  // check input
  if (obj == NULL) { 
    return 1;
  }

  // update location, return success
  obj->x = x;
  obj->y = y;
  return 0;
}

int set_rotation (po_handle obj, float r) {
  // check input
  if (obj == NULL) {
    return 1;
    }

  // update rotation, return success
  obj->r = r;
  return 0;
}

int set_velocity (po_handle obj, float dx, float dy) {
  // check input
  if (obj == NULL) {
    return 1;
    }

  // update velocity, return success
  obj->dx = dx;
  obj->dy = dy;
  return 0;
}

int set_angular_vel (po_handle obj, float dr) {
  // check input
  if (obj == NULL) {
    return 1;
    }

  // update rotation, return success
  obj->dr = dr;
  return 0;

}

/* currently only resolves collisions for circles
 * returns 0 on succes, 1 on failure */
int resolve_collision (po_handle obj1, po_handle obj2){
  // check inputs
  if (obj1 == NULL || obj2 == NULL) {
    return 1;
  }
  // change in x and y
  float d_x = (obj1->x - obj2->x)/2.0;
  float d_y = (obj1->y - obj2->y)/2.0;
  
  // reverse velocities, set location, check for error
  // TODO hopefully, 1s are true, 0s are false
  if (set_velocity(obj1, obj1->dx * -1, obj1->dy * -1) ||
      set_velocity(obj2, obj2->dx * -1, obj2->dy * -1) ||
      set_location(obj1, obj1->x + d_x, obj1->y + d_y) ||
      set_location(obj2, obj2->x + d_x, obj2->y + d_y)) {
    // something went wrong
    return 1;
  }
  // otherwise, no errors
  return 0;
}

void coll_broadphase (world_handle world) {
// min and max keys for the outer array determined by y vals
  int ky_min = dynamic_array_min(world->rows);
  int ky_max = dynamic_array_max(world->rows);
  for (int i = ky_min; i <= ky_max; i++){
    // get the current row
    dynamic_array* cur_row = dynamic_array_get(world->rows, i);
    if (cur_row != NULL) {
      // find our bounds
      int cur_min = dynamic_array_min(cur_row);
      int cur_max = dynamic_array_max(cur_row);
      dynamic_array* next_row = dynamic_array_get(world->rows, i+1);
      if (next_row != NULL) {
	// if both rows contain things, run dection on them
	check_rows(cur_row, cur_min, cur_max, next_row);
      }
      else {
	// we just need to compare this row to itself
	check_row(cur_row, cur_min, cur_max);
      }
    }
  }
}

/* get the area of a polygon */
/* polygon must not be self intersecting */
float poly_area(po_poly polygon){
  float area = 0;
  int max_index = polygon.nvert - 1;
  po_vector* vertex = polygon.vertices;

  // calculate area using summing
  for (int i = 0; i < max_index; i++){
    area += (vertex[i].x * vertex [i+1].y - vertex[i+1].x * vertex[i].y);
  }
  return 0.5 * area;
}

float distance_squared(po_vector point1, po_vector point2){
  return pow(point1.x - point2.x, 2) + pow(point1.y - point2.y, 2);
}

// TODO : best practices for local variables?
// finds and sets the extrema in local coordinates for a given physics object
// polygon may not contain circles
// return 0 on succes, 1 on failure
int get_centroid(po_handle obj) {
  // if our shape is a circle, calculations are trivial
  if (obj->shape.shape_type = 0) {
    obj->centroid = obj->shape.circ.center;
    obj->max_delta = obj->shape.circ.radius;
  }
  else {
    // our shape is a polygon
    po_poly polygon = obj->shape.poly;

    // store the vertices in a slightly more user friendly thing
    po_vector* vertex = polygon.vertices;
    if (vertex == NULL){
      // failure
      return 1;
    }
    // these will hold our centroid sums
    float sum_x = 0;
    float sum_y = 0;

    // the max index of our vertices array
    int max_index = obj->shape.poly.nvert - 1;

    /* calculate the centroid: 
     * cent_x = (1/6A) * Sum((x[i] + x[i+1])*(x[i] * y[i+1] - x[i+1] * y[i])) from 0 to n-1 
     * basically ditto for the y component of centroid */
    for (int i = 0; i < max_index; i++){ 
      // the area component of our sum
      int ar_sum = (vertex[i].x * vertex [i+1].y - vertex[i+1].x * vertex[i].y);
      sum_x += (vertex[i].x + vertex[i+1].x) * ar_sum;
      sum_y += (vertex[i].y + vertex[i+1].y) * ar_sum;	
    }
    // necessaray for the centroid formula
    float six_area_invert = 1 / (6 * poly_area(polygon));
    obj->centroid.x = six_area_invert * sum_x;
    obj->centroid.y = six_area_invert * sum_y;

    // initialize max distance squared
    float max_delta_squared = distance_squared(obj->centroid, vertex[max_index]);

    // loop through the vertices, checking distance again the centroid
    for (int j = 0; j < max_index; j++){
      float cur_dist_squared = distance_squared(obj->centroid, vertex[j]);
      if (cur_dist_squared > max_delta_squared) {
	// update the max distance
	max_delta_squared = cur_dist_squared;
      }
    }
    obj->max_delta = sqrt(max_delta_squared);
  }
  return 0;
}

void coll_midphase(po_handle bucket1, po_handle bucket2){
}
// takes the objects in the hash buckets passed by broadphase
// draws bounding boxes around these objects
// detects overlap between bounding boxes
// if overlap, call narrowphase



/* for circles only */
void coll_narrowphase(po_handle obj1, po_handle obj2){
  // distance squared
  float d_2 = pow((obj1->x - obj2->x), 2.0) + pow((obj1->x - obj2->x), 2.0);

  // sum of radii squared
  float r_2 = pow(obj1->shape.radius,2.0) + pow(obj2->shape.radius,2.0);

  // there's a collision, resolve it
  if(d_2 <= r_2){
    resolve_collision(obj1, obj2);
  }
}

void check_row(dynamic_array* row_k, int k_min, int k_max){
  for (int i = k_min; i < k_max; i++)
  {
    // get current bucket, check for empty
    po_handle cur_kbucket = dynamic_array_get(row_k, i);
    if (cur_kbucket != NULL){

      // if its not empty, get the next bucket in row, check for empty
      po_handle next_kbucket = dynamic_array_get(row_k, i+1);
      if (next_kbucket != NULL) {

	// if that bucket's not empty, go to midphase on this smaller group
	coll_midphase(cur_kbucket, next_kbucket);
      }
    }
  }
}

/* checks two rows for collision
 * does this by calling check row, then comparing the top left square to lower two */
void check_rows(dynamic_array* row_k, int k_min, int k_max, dynamic_array* row_kplus){
  // do collision detection within the row
  check_row(row_k, k_min, k_max);

  // get maxes and mins of our second array
  int kplus_min = dynamic_array_min(row_kplus);
  int kplus_max = dynamic_array_max(row_kplus);

  // find our minimum and maximum indices for looping through both rows
  int min_index = k_min > kplus_min ? k_min : kplus_min;
  int max_index = k_max < kplus_max ? k_max : kplus_max;

  // loop through the top list, compare top right bucket to adjacent buckets
  for (int i = min_index; i < max_index; i++) {
    // get current bucket, check for empty
    po_handle cur_kbucket = dynamic_array_get(row_k, i);
    if (cur_kbucket != NULL)
    {
      // get bucket directly below, check for empty
      po_handle cur_plusbucket = dynamic_array_get(row_kplus, i);
      if (cur_plusbucket != NULL) {
	// move to the next step in collision resolution with these buckets
	coll_midphase(cur_kbucket, cur_plusbucket);
      }
      // get bucket below and to the right, check for empty
      po_handle next_plusbucket = dynamic_array_get(row_kplus, i+1);
      if (next_plusbucket != NULL){
	coll_midphase(cur_kbucket, next_plusbucket);
      }
    }
  }
}

void detect_collision();
// check collision for every clock tick
// really, check if things are collided or would have been collided

// loop through two rows and find things that could collide
// calls midphase on anything that could
// neither array should be NULL
// seperating axis theorem on objects that might collide
// if collision, call resolve collsion (with two objects)? set collision flag?
