/**************************************************************
 * Physics.c - the actual game engine
 *
 * Creates/removed physics objects
 * Updates physics world  
 **************************************************************/
#include <stdlib.h>
#include <stddef.h>
#include "physics.h"

int MAX_HASH_LEN = 1021;

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

  po_geometry object;
  
  // allows for linked lists within the hash table
  po_handle next;
} po_imp;

/* for our silly linked list version of world 
 * this will not be a thing later, hopefully
/* actual implementation of physics world structure */
typedef struct world_t {
  po_handle object;
  
  // buckets are 500 pixels in size
  int bucket_size;
  struct world_t* next;
  //po_handle quadrant1[];
  //po_handle quadrant2[];
  //po_handle quadrant3[];
  //po_handle quadrant4[];
  
  // a place holding hash table while we sort out the quadrants
  po_handle hash_table[1000];
} world_t;

// takes coordinates in global coords, returns pointer to proper bucket
int spatial_hash (int x, int y, world_handle world) {
  // do the thing that returns us a pointer in our quadrant
  // need to check signs of x and y?

  return 1;

}

void add_to_world (world_handle world, po_handle object){
  int x = object->x;
  int y = object->y;
  
  if (x > 0 && y >= 0){
    //hashing will lead to quadrant 1
  }
  else if (x <= 0 && y >= 0){
    // hashing leads to quadrant 2
  } 
  else if (x <= 0 && y < 0){
    // hashing leads to quad 3
  }
  else{
    // hashing leads to quad 4
  }

  // determine bucket the object belongs in
  int key = spatial_hash(object->x, object->y, world);

  // make a temporart object to hold value at that location
  po_handle temp = world->hash_table[key];
  if (temp == NULL){
    // nothing there, insert at that location
    world->hash_table[key] = object;
  }
  else {
    // place object at beginning of linked list
    object->next = temp;
    world->hash_table[key] = object;
  }
}

/* create a new world */
/* basically, just an empty piece of memory */
world_handle new_world ()
{
  //world_t* new_first = malloc(sizeof(world_t));
  //new_first->object = NULL;
  //new_first->next = NULL;
  return NULL;
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
int resolve_collision (po_handle obj1, po_handle obj2){
  obj1-> dx * -1;
  obj1-> dy * -1;
  obj2-> dx * -1;
  obj2-> dy * -1;
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

void coll_narrowphase();
// parallel axis theorem on objects that might collide
// if collision, call resolve collsion (with two objects)? set collision flag?
