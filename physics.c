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

int MAX_HASH_LEN = 1021;

// number of pixels per bucket in the spatial hash
int BUCKET_SIZE = 500;
int INIT_SIZE = 10;

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
  po_geometry object;
  
  // allows for linked lists within the hash table
  po_handle next;
} po_imp;

/* dat spatial hash though */
typedef struct world_t {
  // the quadrants of our world - arrays of linked lists
  // these are spatially hashed and correspond to quadrants on an xy plane
  po_handle* quad1;
  po_handle* quad2;
  po_handle* quad3;
  po_handle* quad4;
} world_t;

// accepts the origin of a physics object and returns a pointer to a bucket
// in the world
po_handle* spatial_hash (int x, int y, world_handle world) {
  if (x > 0 && y >= 0) // object's in the first quadrant
    {} 
  else if (x <= 0 && y >= 0) // it's in the second quadrant
    {}
  else if (x <= 0 && y < 0) // it's in the third quadrant
    {}
  else // it's in the fourth quadrant
    {}
}

/* create a new world */
/* basically, just an empty piece of memory */
world_handle new_world ()
{
  // make new world
  world_handle world;
  
  // figure out amount of memory to malloc
  int size = INIT_SIZE * sizeof(po_handle);

  // initialize the size of our quadrants
  world->quad1 = malloc(size);
  world->quad2 = malloc(size);
  world->quad3 = malloc(size);
  world->quad4 = malloc(size);

  // return this handle
  return world;
}

/* add object to the physics world */
po_handle add_object (world_handle world, po_geometry* geom, 
		      float x, float y, float r)
{
  // make a new object
  po_handle new_obj;
  new_obj->x = x;
  new_obj->y = y;
  new_obj->r = r;
  new_obj->dx = 0;
  new_obj->dy = 0;
  new_obj->dr = 0;
  new_obj->object = *geom; 
  new_obj->next = NULL;

}

int remove_object (world_handle world, po_handle obj){
  if (world == NULL || obj == NULL) {
    return 1;
    }
}
/* Updates object's global position based on velocity
 * Future versions may include more sophistocated algorthims using acceleration */
void integrate (float dx, float dy, float dr, float time_step, po_handle obj) {
  // check input
  if (obj == NULL) { 
    return; // if you want to return an error value then change the function declaration
  }

  // apply euler's method
  obj->x = obj->x + (dx * time_step);
  obj->y = obj->y + (dy * time_step);
  obj->r = obj->r + (dr * time_step);
  return; // ditto
}

int set_location (float x, float y, po_handle obj) {
  // check input
  if (obj == NULL) { 
    return 1;
  }

  // update location, return success
  obj->x = x;
  obj->y = y;
  return 0;
}

int set_rotation (float r, po_handle obj) {
  // check input
  if (obj == NULL) {
    return 1;
    }

  // update rotation, return success
  obj->r = r;
  return 0;
}

int set_velocity (float dx, float dy, po_handle obj) {
  // check input
  if (obj == NULL) {
    return 1;
    }

  // update velocity, return success
  obj->dx = dx;
  obj->dy = dy;
  return 0;
}

int set_angular_vel (float dr, po_handle obj) {
  // check input
  if (obj == NULL) {
    return 1;
    }

  // update rotation, return success
  obj->dr = dr;
  return 0;

}
/* currently only resolves collisions for circles */
int resolve_collision (po_handle obj1, po_handle obj2){
  if (obj1 == NULL || obj2 == NULL) {
    return 1;
    }
  obj1-> dx * -1;
  obj1-> dy * -1;
  obj2-> dx * -1;
  obj2-> dy * -1;

  float delta_x = (obj1->x - obj2->x)/2.0; 
  float delta_y = (obj1->y - obj2->y)/2.0;
  obj1->x += delta_x;
  obj1->y += delta_y;
  obj2->x -= delta_x;
  obj2->y -= delta_y; 
  return 0;
}

void detect_collision();
// check collision for every clock tick
// really, check if things are collided or would have been collided


void coll_broadphase (world_handle world);
// create new, smaller world -> a bucket if you will
// (what exactly this is depends on world implementation)
// only possible things that will collide
// if no object in bucket, or one object in bucket, skip
// call midphase on each of the small worlds we've created

void coll_midphase();
// takes the objects in the hash buckets passed by broadphase
// draws bounding boxes around these objects
// detects overlap between bounding boxes
// if overlap, call narrowphase 

void coll_narrowphase(po_handle obj1, po_handle obj2){
  float d_2 = pow((obj1->x - obj2->x), 2.0) + pow((obj1->x - obj2->x), 2.0);
  float r_2 = pow(obj1->object.radius,2.0) + pow(obj2->object.radius,2.0);
  if(d_2 <= r_2){
    resolve_collision(obj1, obj2);
  }
}
// seperating axis theorem on objects that might collide
// if collision, call resolve collsion (with two objects)? set collision flag?
