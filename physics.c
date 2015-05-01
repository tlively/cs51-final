/**************************************************************
 * Physics.c - the actual game engine
 *
 * Creates/removed physics objects
 * Updates physics world  
 *
 * structure of file designed to maintain maximum abstraction
 **************************************************************/
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include "physics.h"
#include "dynamic_array.h"

#define DEBUG
#ifdef DEBUG
#define LOG(args...) do {printf(args);}while(0);
#else
#define LOG(args...)
#endif

// make some macros
#define SHAPE_TYPE(obj) (obj->shape.shape_type)
#define POLY(obj) (obj->shape.poly)
#define CIRC(obj) (obj->shape.circ)
#define VERTEX(obj) (obj->shape.poly.vertices)
#define NVERTS(obj) (obj->shape.poly.nvert)

// number of pixels per bucket in the spatial hash
#define BUCKET_SIZE 500
#define MY_PI 3.1415926535

// for some reasont this is being a struggle,

/*********************************************************
 * Structures
 ********************************************************/

/* the actual implementation of a physics object structure */
typedef struct po_imp {

  // linear components (origin in global coords)
  po_vector origin;
  po_vector vel;
  po_vector force;

  // angualar positive r (rotations) is counterclockwise
  float r;
  float dr;

  // force and torque vectors
  float torque;

  // the actual shape
  po_geometry shape;

  // area of an obj
  float area;
  
  // moment of inertia
  float moment;
 
  // centroid of the object in local coordinates of poly
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

/*************************************************************
 * Helpers Header
 ************************************************************/

/* create a circle */
po_circle create_circ(po_vector center, float radius){
  po_circle circ;
  circ.center = center;
  circ.radius = radius;
  return circ;
}

/* create a poly */
po_poly create_poly(po_vector* vertices, int nvert){
  po_poly poly;
  poly.vertices = vertices;
  poly.nvert = nvert;
  return poly;
}

/* create geometry with polygon, hide our dirty laundry */
po_geometry create_geom_poly(po_poly poly, float density){
  po_geometry geom;
  geom.shape_type = 1;
  geom.poly = poly;
  geom.density = density;
  return geom;
}

/* create geometry with circle */

po_geometry create_geom_circ(po_circle circ, float density){
  po_geometry geom;
  geom.shape_type = 0;
  geom.circ = circ;
  geom.density = density;
  return geom;
}

/* get the distance between two points squared */
float distance_squared(po_vector point1, po_vector point2);

/* calculate the centroid of a polygon */
po_vector get_poly_centroid (po_handle poly);

/* sets the centroid, max distanc from centroid in local coordinates, and area
 * return 0 on succes, 1 on failure */
int set_centroid_area(po_handle obj);

/* converts a (polygon) centroid to global coordinates
 * accepts a centroid in local coords and origin in global coords */
po_vector get_centroid_global(po_vector cent, po_vector origin);

/* gets the moment of inertia for a polygon */
int moment_of_inertia(po_handle obj);

/* gets the area for a polygon */
float poly_area(po_poly obj);

/* checks concavity and order of points
 * returns 0 if convex, 1 if oriented improperly or concave */
int check_concavity (po_handle obj);

/* returns the vertices of an polygon in global coordinates */
void get_global_coord (po_handle obj, po_vector** global_vertices);

/*************************************************************
 * User Related - Initiation
 ************************************************************/

/* create a new world 
 * returns world on success, NULL on failure */
world_handle new_world () {
  // make new world
  world_handle world = malloc(sizeof(world_t));
  world->rows = dynamic_array_create();
  return world;
}

/* add object to the physics world */
po_handle add_object (world_handle world, po_geometry* geom, 
		      float x, float y, float r) {
  // make a new object
  po_handle new_obj = malloc(sizeof(po_imp));
  new_obj->origin.x = x;
  new_obj->origin.y = y;
  new_obj->r = r;
  new_obj->vel.x = 0;
  new_obj->vel.y = 0;
  new_obj->dr = 0;
  new_obj->force.x = 0;
  new_obj->force.y = 0;
  new_obj->shape = *geom;

  if(geom->shape_type && (check_concavity(new_obj) || set_centroid_area(new_obj)))
  {
    // strugs - either fails concavity failure to set cetroid
    return NULL;
  }

  //reject object arbitrarily if it is too large
  if(new_obj->max_delta > BUCKET_SIZE/2) 
    {
      return NULL;
    }

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
  return new_obj;
}

int set_location (po_handle obj, float x, float y) {
  // check input
  if (obj == NULL) {
    LOG("NULL pointer exception in physics.c\n NULL pointer passed in function set_location"); 
    return 1;
  }

  // update location, return success
  obj->origin.x = x;
  obj->origin.y = y;
  return 0;
}

int set_rotation (po_handle obj, float r) {
  // check input
  if (obj == NULL) {
    LOG("NULL pointer exception in physics.c: NULL pointer passed in function set_rotation"); 
    return 1;
    }

  // update rotation, return success
  obj->r = r;
  return 0;
}

int set_velocity (po_handle obj, float dx, float dy) {
  // check input
  if (obj == NULL) {
    LOG("NULL pointer exception in physics.c: NULL pointer passed in function set_velocity"); 
    return 1;
    }
  // update velocity, return success
  obj->vel.x = dx;
  obj->vel.y = dy;
  return 0;
}
//gets velocity of object
po_vector get_velocity (po_handle obj){
  po_vector velocity;
  velocity.x = obj->vel.x;
  velocity.y = obj->vel.y;
  return velocity;
}
//gets position in global coordinates
po_vector get_position (po_handle obj){
  po_vector position;
  position.x = obj->vel.x;
  position.y = obj->vel.y;
  return position;
}

int set_angular_vel (po_handle obj, float dr) {
  // check input
  if (obj == NULL) {
    LOG("NULL pointer exception in physics.c\n NULL pointer passed in function set_angular_vel"); 
    return 1;
    }

  // update rotation, return success
  obj->dr = dr;
  return 0;

}

/* remove object from the world */
int remove_object (world_handle world, po_handle obj){
  if (world == NULL) {
    LOG("NULL pointer exception in physics.c: remove_object (world_handle world, po_handle obj): world is NULL"); 
    return 1;
  }
  else if(obj == NULL){
    LOG("NULL pointer exception in physics.c: remove_object (world_handle world, po_handle obj): obj is NULL"); 
    return 1;
  }
  // get indexes
  int kx = obj->origin.x / BUCKET_SIZE;
  int ky = obj->origin.y / BUCKET_SIZE;

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

/************************************************************
 * Integration
 ************************************************************/

/* Updates object's global position based on velocity
 * Future versions may include more sophistocated algorthims using acceleration
 * error checking should happen prior to passing things in */
void integrate (po_handle obj, float time_step) {
  // apply euler's method (the most logical choice since we have no accel)

  // use F=ma to get acceleration: a = m / F
  po_vector a = vect_scaled(vect_recip(obj->force), obj->area * obj->shape.density);
  
  // update velocities: f = ma, dr = moment * torque 
  obj->vel = vect_scaled(obj->force, 
			      time_step / (obj->area * obj->shape.density));
  obj->dr += obj->torque * (1.0 / obj->moment) * time_step;

  // update position
  obj->origin = vect_add(obj->origin, vect_scaled(obj->vel, time_step));
  obj->r += obj->dr * time_step;
}
//TODO :update world
int update (world_handle world, float dt){

  for(int i = dynamic_array_min(world->rows), maxi = dynamic_array_max(world->rows); i <= maxi; i++)
  {
    dynamic_array* rows = dynamic_array_get(world->rows,i);
    if (rows == NULL) {
      continue;
    }
    else 
    {
      for(int j = dynamic_array_min(rows), maxj = dynamic_array_max(rows); j <= maxj; j++)
      {
        // removed an object from a row
        po_handle current_obj = dynamic_array_remove(rows,j);
        integrate(current_obj,dt);
        po_handle next_obj = current_obj -> next;
        while (next_obj != NULL) 
        {
          integrate(next_obj, dt);
          next_obj = next_obj->next;
        }
      }
    }
  }
}
/*************************************************************
 * Collisions - Header
 *************************************************************/

/************ collision detection **************************/
void coll_broadphase(world_handle world);
void coll_narrowphase(po_handle obj1, po_handle obj2);
void coll_midphase(po_handle bucket1, po_handle bucket2);

/************ collision resolution *****************************/
int resolve_coll_circs (po_handle circ1, po_handle circ2);
int resolve_coll_polys (po_handle circ1, po_handle circ2, int run_once);
int resolve_coll_mixed (po_handle poly, po_handle circ);

/************************************************************
 * Collsion Detection - Broadphase, Spatial Hashing
 ************************************************************/
 
/* helper: checks for objects in adjacent buckts along row */
void check_row(dynamic_array* row_k, int k_min, int k_max);

/* helper: checks adjacent buckets in two rows for collision */
void check_rows(dynamic_array* row_k, int k_min, int k_max, 
    dynamic_array* row_kplus);

/* goes through all the things in our spatial hash */
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

/**************** Broadphase Helpers ************************/

/* helper function for broadphase collision detect 
 *  checks for objects in adjacent buckts along row */
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

/* helper function for collision broadphase
 * checks adjacent buckets in two rows for collision */
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

/************************************************************
 * Collision Detection - Midphase, Bounding Boxes
 ************************************************************/

/* helper: returns 1 if bounding boxes of two objects collide */
int check_bounding (po_handle obj1, po_handle obj2);

/* use bounding boxes to narrow down collisions further */
void coll_midphase(po_handle bucket1, po_handle bucket2) {
  po_handle cur_obj = bucket1;
  po_handle secondary_list = bucket2;
  while(cur_obj != NULL){
    // check collision with other things in the bucket
    for (po_handle next_obj = cur_obj->next; next_obj != NULL; 
	 next_obj = next_obj->next){
      if (check_bounding(cur_obj, next_obj)){
	coll_narrowphase(cur_obj, next_obj);
      }
    }  
    // check for collisions with stuff in the other bucket
    for (po_handle next_b2 = bucket2; next_b2 != NULL; 
	 next_b2 = next_b2->next) {
      if (check_bounding(cur_obj, next_b2)){
	coll_narrowphase(cur_obj, next_b2);
      }
    }
    cur_obj = cur_obj->next;
  }
}

/****************** Midphase Helpers ***********************/

/* midphase helper function returns 1 if bounding boxes collide */
int check_bounding (po_handle obj1, po_handle obj2) {

  // get our max widths/heights
  float summed_deltas = 2 * (obj1->max_delta + obj1->max_delta);
  
  // convert centroid to global coordinates if it's a poly
  po_vector cent1 = SHAPE_TYPE(obj1) ? 
    get_centroid_global(obj1->centroid, obj1->origin) 
    : obj1->centroid;
  po_vector cent2 = SHAPE_TYPE(obj2) ? 
    get_centroid_global(obj2->centroid, obj2->origin) 
    : obj2->centroid;
  
  // use bounding boxes to do collisiion detection
  return (abs(cent1.x - cent2.x) * 2 < summed_deltas 
	  && abs(cent1.y - cent2.y) * 2);
}

/************************************************************
 * Collsion Detection - Narrowphase, Parallel Axis Theorem
 ************************************************************/

/* helper: find the points associated with min and max dot product with axis
 * first value in array is min, second is max
 * updated pointers that are passed in to point to min and max vals */
int vect_dot_extrema(po_handle obj, po_vector axis, float* min, float* max);

/* helper: checks for collision, returns 1 on collision, 0 on none*/
int sep_axis(po_handle obj1, po_handle obj2);

/* helper: for detecting collsions on circles and polygons 
 * returns 1 on collision, 0 on none */
int coll_poly_circ(po_handle poly, po_handle circ);

/* detects tiny collisions depending on shape
 * uses separating axis theorem for polygons */ 
void coll_narrowphase(po_handle obj1, po_handle obj2) {
  // check the object types and act accordingly
  if (!SHAPE_TYPE(obj1) && !SHAPE_TYPE(obj2)) {  
    // sum of radii squared
    float r_2 = pow(CIRC(obj1).radius,2.0) + pow(CIRC(obj2).radius,2.0);

    // there's a collision, resolve it
    if(distance_squared(obj1->centroid,obj2->centroid) <= r_2){
      resolve_coll_circs(obj1, obj2);
    } 
  }
  else if (SHAPE_TYPE(obj1) && SHAPE_TYPE(obj2)) {
    // two polygons are colliding
    if (sep_axis(obj1, obj2) || sep_axis(obj2, obj1)) {
      resolve_coll_polys(obj1,obj2,0);
    }
  }
  else {
    // we have a poly - circle collision
    if(SHAPE_TYPE(obj1)){
      // first object is a poly, second is a circle
      if (coll_poly_circ(obj1, obj2)){
	resolve_coll_mixed(obj1,obj2);
      }
    }
    else {
      // first object is a circ, second poly
      if (coll_poly_circ(obj2, obj1)){
	resolve_coll_mixed(obj2,obj1);
      }
    }
  }
}

/************** Narrowphase Helpers *************************/

/* find the points associated with min and max dot product with axis
 * first value in array is min, second is max
 * updated pointers that are passed in to point to min and max vals */
int vect_dot_extrema(po_handle obj, po_vector axis,float* min, float* max) {
  if(obj == NULL){
    LOG("NULL pointer exception in physics.c");
    return 1;
  }
  // initialize extrema (and do a lot of pointer magic)
  po_vector* global_vertex;
  get_global_coord(obj, &global_vertex);
  *min = vect_dot_prod(global_vertex[0], axis);
  *max = *min;

  // got through and find min and max coords, updating as we go
  for (int i = 1; i < NVERTS(obj); i++){
    float dot_product = vect_dot_prod(axis, global_vertex[i]);
    if (dot_product < *min){
      *min = dot_product;
    }
    else if(dot_product > *max){
      *max = dot_product;
    }
  }
  return 0;
}

/* checks for collision, returns 1 on collision, 0 on none 
 * uses all axis associated with obj1 for the parallel axis theorem */
int sep_axis(po_handle obj1, po_handle obj2) {
  // go through all the axis on our stuffs
  for (int i = 0, j = 1; i < NVERTS(obj1); i++, j = (j+1) % NVERTS(obj1)) {
  
    // get the normal to one of the sides on obj1
    po_vector axis = vect_axis(VERTEX(obj1)[i],(VERTEX(obj1)[j]));

    // get the min and max projections
    float min1, max1, min2, max2;
    vect_dot_extrema(obj1, axis, &min1, &max1);  
    vect_dot_extrema(obj2, axis, &min2, &max2);

    // if there's a space between...
    if (max1 < min2 || min1 > max2) {
      // the objects have definitely not collided
      return 0;
    }			       
  }  
  // there was no separation, the objects have collided
  return 1; 
}

/* for detecting collsions on circles and polygons 
 * returns 1 on collision, 0 on none */
int coll_poly_circ(po_handle poly, po_handle circ){
  if(poly == NULL || circ == NULL)
  {
    LOG("NULL pointer exception in physics.c");
    return 1;
  }
  // get all the sides of the poly
  po_vector* global_verts;
  get_global_coord(poly, &global_verts);

  // go through all sides of the polygon!
  for (int i = 0, j = 1; i < NVERTS(poly); i++, j = (j+1) % NVERTS(poly)) {
    // get the vector of the current side
    po_vector side = vect_from_points(global_verts[i], global_verts[j]);

    // get the vector from the vertex to the circle's center (hypotenus)
    po_vector circ_to_vert = vect_from_points(global_verts[i], circ->centroid);

    // get the projection of this on to the side (base)
    po_vector proj = vect_project(circ_to_vert, side);

    // now we're doing pythagorean's theorem to get dist from point to line
    if (vect_mag_squared(circ_to_vert) - vect_mag_squared(proj) 
	> pow(CIRC(circ).radius,2)) {
      // length of the final side is greater than the radius, no collision
      return 0;
    }
  }
  // there's a collision
  return 1;
}

/************************************************************
 * Collision Resolution
 ************************************************************/

/* helper function to get unit normal vectors */
void get_normals (po_vector* verts, int size, po_vector** normals) {
  *normals[size];
  for (int i = 0, j = 1; i < size; i++, j = (j+1) % size){
    *normals[i] = vect_unit(vect_axis(verts[i], verts[j]));
  }
} 

/* finds point of collision in global coords
   takes point and vector in global coords */
po_vector get_coll_pt(po_vector point, po_vector side_origin, po_vector side_end) {
  // get the projection of the vector connecting the colliding vertex onto the line
  po_vector proj = vect_project(vect_from_points(side_origin, point), 
				vect_from_points(side_origin, side_end));

  // get the point this hits on the side in global coords
  po_vector intersect_point;
  intersect_point.x = proj.x + side_origin.x;
  intersect_point.y = proj.y + side_origin.y;
  
  return intersect_point;
}

/* takes point and vector in global coords
 * this vector will point in the direction of force on the point */
po_vector get_force_vector(po_vector point, po_vector intersect_point) {
  // get the force vector! (it's the vector from the vertice to the closest part of the side)
  return vect_from_points(point, intersect_point);
}

/* gets the torque  on an object
 * accepts an object and the point of collision */
float get_torque(po_vector point, po_handle poly) {
  // the cross prod of the vector from the center of the poly to the point 
  // with the angular velocity 
  po_vector r = vect_from_points(get_centroid_global(poly->centroid, poly->origin), point);
  return sqrt(vect_mag_squared(r)) * vect_cross_prod(r, poly->vel);
}

/* go through the sides of poly1 comparing with the verts of poly2 
 * to get the vertex that is poking through 
 * returns 1 on failure, 0 on success
 * if we don't find anything, we need to switch inputs and try again
 * if run_once = true, the function will terminate after one round, 
 * else it will switch input and go again (to prove we can recurse even if not normally in c) */
int resolve_coll_polys (po_handle po_pts, po_handle po_sides, int run_once) {

  // the polygon we're doing corner stuff with 
  po_vector* vert_pts;
  get_global_coord(po_pts, &vert_pts);
  
  // the polygon we're doing line stuff with
  po_vector* vert_sides;
  po_vector* normals;
  get_global_coord(po_sides, &vert_sides);
  get_normals(vert_sides, NVERTS(po_sides), &normals);

  // the outer loop is for the vertices of the first poly
  for (int i = 0; i < NVERTS(po_pts); i++){
    // some initializers 
  int i_s = 0, min_dot_prod = 0, max_j = NVERTS(po_sides) - 1;

  // go through the vertices of po_pts (i_s is the index of min dot prod)
  for (int j = 0; j < NVERTS(po_sides); j++) {      
    // get the vector from a corner to 
    po_vector corner_to_point = vect_from_points(vert_pts[i], vert_sides[j]);
    float cur_dot_prod = vect_dot_prod(normals[j], corner_to_point);
    if (0 > cur_dot_prod){
      // no intersection, skip the rest of the dot prods
      break;
    }
    else if (-cur_dot_prod > min_dot_prod){
      // we've found a new min value!
      i_s = j;
      min_dot_prod = -cur_dot_prod;
    }
    
    // we've made it to the end without breaking...
    if (j == max_j) {     
      // get point on sides_poly where collision is happening
      po_vector side_point = get_coll_pt(vert_pts[i],vert_sides[i_s], 
					 vert_sides [(i_s+1) % NVERTS(po_sides)]);
      
      // update force information
      po_pts->force = get_force_vector(vert_pts[i], side_point);;
      po_sides->force = vect_scaled(po_pts->force,-1);

      // update update torque
      po_pts->torque = get_torque(vert_pts[i], po_pts);
      po_sides->torque = get_torque(side_point, po_sides);
      
      // update position information
      po_pts->origin = vect_add(po_pts->origin, vect_scaled(po_pts->force,0.5)); 
      po_sides->origin = vect_add(po_sides->origin, vect_scaled(po_sides->force,0.5));
      
      // success!
      return 0;
    }
  }
    // else, carry on
  }
  if (run_once == 0) {
     // we didn't get anything the first time so we need switch and go again
    resolve_coll_polys(po_sides, po_pts, 1);
  }
  return 1;
}

/* make this a thing: takes a poly and a circ and resolves the collision */
int resolve_coll_mixed (po_handle poly, po_handle circ){
  po_vector origin = get_centroid_global(circ->centroid, circ->origin);
  
  // the polygon we're doing line stuff with
  po_vector* vertex;
  po_vector* normals;

  get_global_coord(poly, &vertex);
  get_normals(vertex, NVERTS(poly), &normals);
  float rad_squared = pow(CIRC(circ).radius, 2);

  for (int i = 0, j = 1; i < NVERTS(poly); i++, j = (j+1) % NVERTS(poly)){

    po_vector coll_pt = get_coll_pt(origin, vertex[i], vertex[j]);
    if (distance_squared(origin, coll_pt) < rad_squared){
      // we have collision!
            // update force information
      circ->force = get_force_vector(circ->force, coll_pt);
      poly->force = vect_scaled(poly->force,-1);

      // update update torque
      poly->torque = get_torque(coll_pt, poly);
      
      // update position information
      circ->origin = vect_add(circ->origin, vect_scaled(circ->force,0.5)); 
      poly->origin = vect_add(poly->origin, vect_scaled(poly->force,0.5));
    }
  }

  // if we get a one, we're successful!
  //return !(update_resolve_polys(circ, poly, circ->origin, vert_sides, normals));
}

/* resolves collision between two circles 
 * returns 0 on success, 1 on failure */ 
int resolve_coll_circs (po_handle circ1, po_handle circ2){
  // check inputs
  if (circ1 == NULL || circ2 == NULL) {
    LOG("NULL pointer exception in physics.c");
    return 1;
  }
  // get centroids in global coords
  po_vector cent1 = get_centroid_global(circ1->centroid, circ1->origin);
  po_vector cent2 = get_centroid_global(circ1->centroid, circ1->origin);
  
  // get the vector conecting them
  circ1->force = vect_add_scalar(vect_from_points(cent1, cent1), 
				 - CIRC(circ1).radius - CIRC(circ2).radius);
  circ2->force = vect_scaled(circ1->force, -1);

  // move them apart
  circ1->origin = vect_add(circ1->origin, vect_scaled(circ1->force, 0.5));
  circ2->origin = vect_add(circ2->origin, vect_scaled(circ2->force, 0.5));
}

/***************************************************************
 * Helper Functions
 **************************************************************/

/* get the distance between two points squared */
float distance_squared(po_vector point1, po_vector point2){
  return pow(point1.x - point2.x, 2) + pow(point1.y - point2.y, 2);
}

/* calculate the centroid and area of a polygon 
 * cent_x = (1/6A) * Sum((x[i] + x[i+1])*(x[i] * y[i+1] - x[i+1] * y[i])) 
 * from [0,n-1] 
 * basically ditto for the y component of centroid */
po_vector get_poly_centroid (po_handle poly){

  // our centroid sum and area sums, initiallized at 0
  po_vector sum;
  sum.x = 0;
  sum.y = 0;
  poly->area = 0;

  // do our sum 
  for (int i = 0, max_index = NVERTS(poly) - 1; i < max_index; i++){ 
    // the area component of our sum
    float ar_sum = (VERTEX(poly)[i].x * VERTEX(poly)[i+1].y 
		  - VERTEX(poly)[i+1].x * VERTEX(poly)[i].y);
    sum.x += (VERTEX(poly)[i].x + VERTEX(poly)[i+1].x) * ar_sum;
    sum.y += (VERTEX(poly)[i].y + VERTEX(poly)[i+1].y) * ar_sum;
    poly->area += ar_sum;
  }
  // necessaray for the centroid formula
  float sixth_inverted_area = 1 / (6 * poly->area);
  sum.x = sixth_inverted_area * sum.x;
  sum.y = sixth_inverted_area * sum.y;
}

// sets the centroid and max distanc from centroid in local coordinates
// polygon may not contain circles
// return 0 on succes, 1 on failure
int set_centroid_area(po_handle obj) {
  // if our shape is a circle, calculations are trivial
  if (SHAPE_TYPE(obj) = 0) {
    obj->centroid = CIRC(obj).center;
    obj->max_delta = CIRC(obj).radius;
    obj->area = MY_PI * pow(CIRC(obj).radius, 2);
  }
  else {
    // our shape is a polygon
    if (VERTEX(obj) == NULL) {
      // something went horribly wrong
      LOG("NULL pointer exception in physics.c");
      return 1;
    }
    // get the centroid! (in local coords)
    obj->centroid = get_poly_centroid(obj);

    // initialize max distance squared, loop through vertices, get max distance
    float max_delta_squared = distance_squared(obj->centroid, VERTEX(obj)[0]);

    for (int i = 1, max_index = obj->shape.poly.nvert; i < max_index; i++){
      float cur_dist_squared = distance_squared(obj->centroid, VERTEX(obj)[i]);
      if (cur_dist_squared > max_delta_squared) {
	// update the max distance
	max_delta_squared = cur_dist_squared;
      }
    }
    // update our object
    obj->max_delta = sqrt(max_delta_squared);
  }
  return 0;
}

/* converts a centroid to global coordinates
 * accepts a centroid in local coords and origin in global coords */
po_vector get_centroid_global(po_vector cent, po_vector origin) {
  cent.x = cent.x + origin.x;
  cent.y = cent.y + origin.y;
  return cent;
}

/* makes sure the shape is convex points are order in a counter clockwise direction
 * returns 0 if convex, 1 if oriented improperly or concave */
int check_concavity (po_handle obj){
  //iterates through vertices generating two vectors
  for(int i = 0, j = 1, k = 2; i < NVERTS(obj); 
      i++, j = (j+1) % NVERTS(obj), k = (k+1) % NVERTS(obj)) {
    // makes sure the vector cross products are positive (idicating correctness)
    po_vector side1 = vect_from_points(VERTEX(obj)[i], VERTEX(obj)[j]);
    po_vector side2 = vect_from_points(VERTEX(obj)[j], VERTEX(obj)[k]);

    if (vect_cross_prod(side1, side2) <= 0) {
      // either point are out of order, or this is concave
      return 1;
    }
  }
  return 0;
}


/* returns the vertices of an obj in global coordinates */
void get_global_coord (po_handle obj, po_vector** global_vertices){

  // creates the rotation matrix with specified r
  matrix rotation_matrix = 
    vect_create_matrix(cos(obj->r), -sin(obj->r), sin(obj->r), cos(obj->r));
  
  // give our array a size!
  *global_vertices[NVERTS(obj)];

  po_vector global_centroid = get_centroid_global(obj->centroid, obj->origin);
  for(int i =0; i < NVERTS(obj); i++)
  {
    // gets the coordinates of the vertice with the centroid as the origin
    po_vector point = vect_from_points(obj->centroid, obj->shape.poly.vertices[i]);

    // "rotates the current point with the transformation matrix"
    point = vect_matrix_mult(point, rotation_matrix);

    //converts the local vertices to global coordinates
    point.x = global_centroid.x + point.x;
    point.y = global_centroid.y + point.y;
    *global_vertices[i] = point;
  }
}
