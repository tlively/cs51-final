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
  po_geometry object;

  // linked lists stuffs -> the next object in list
  po_handle next;
} po_imp;

/* dat spatial hash though */
typedef struct world_t {
  dynamic_array* contents;
} world_t;

/* accepts the origin of a physics object in global coords
 * returns a pointer to a bucket
 * use the cantor pairing function to map tuples to single ints */
po_handle spatial_hash (float x, float y, world_handle world) {
  float k1 = x/BUCKET_SIZE;
  float k2 = y/BUCKET_SIZE;
  
  // cantor key uniquely maps two values to a single value
  // in our case, values within BUCKET_SIZE chunks will be mapped to the same bucket
  int key = .5*(k1+k2)*(k1+k2+1)+k2;
  
  return dynamic_array_get(world->contents,key);
}

/* create a new world 
 * returns world on success, NULL on failure */
world_handle new_world ()
{
  // make new world
  world_handle world;
  world->contents = dynamic_array_create();

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
  
  // the pointer struggles are really real here...
  // add to world
  po_handle at_index = spatial_hash(x, y, world);
  if (at_index == NULL){
    // insert at this index in array
    at_index = new_obj;
    new_obj->next = NULL;
  }
  else {
    // I'm not sure about how to get this nonsence to work like a normal table/linked list thing
    // insert new_obj at head of linked list
    po_handle temp = at_index;
    new_obj->next = at_index;
    at_index = new_obj;
    // hashtable[key] = new_obj;
    
  }
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
