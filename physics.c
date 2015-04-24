/**************************************************************
 * Physics.c - the actual game engine
 *
 * Creates/removed physics objects
 * Updates physics world  
 **************************************************************/
#include <stddef.h>
#include "physics.h"

/* the actual implementation of a physics object structure */
struct po_imp {
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
}po_imp;

/* for our silly linked list version of world 
 * this will not be a thing later, hopefuly */
/* actual implementation of physics world structure */
struct world_t {
  /*struct po_handle object;
    struct world_t* next;*/
  int placeholder;
} world_t;

typedef world_t* world_handle;
typedef po_imp* po_handle;

/* create a new world */
world_handle new_world(){

  return NULL;
}

/* add object to the physics world */
/*
po_handle add_object (world_handle world, po_geometry* geom, 
		      float x, float y, float r)
{
  // make a new object
  struct po_handle new_obj;
  new_obj.x = x;
  new_obj.y = y;
  new_obj.r = r;
  new_obj.dx = 0;
  new_obj.dy = 0;
  new_obj.dr = 0;
  new_obj.object = geom; 
  
  struct world_handle new_first;
  new_first.object = new_obj;
  new_first.next = world;
  return new_world;
}
*/

/* Updates object's global position based on velocity
 * Future versions may include more sophistocated algorthims using acceleration */
void integrate (float dx, float dy, float dr, float time_step, po_handle obj) {
  // apply euler's method
  
  obj.x = obj.x + (dx * time_step);
  obj.y = obj.y + (dy * time_step);
  obj.r = obj.r + (dr * time_step);
}

int set_location (float x, float y, po_handle obj) {
  /* // check input
  if (obj == NULL) { 
    return 1;
  } */

  // update location, return success
  *obj.x = x;
  *obj.y = y;
  return 0;
}

int set_rotation (float r, po_handle obj) {
  /*  // check input
  if (obj == NULL) {
    return 1;
    } */

  // update rotation, return success
  obj.r = r;
  return 0;
}

int set_velocity (float dx, float dy, po_handle obj) {
  /*  // check input
  if (obj == NULL) {
    return 1;
    } */

  // update velocity, return success
  obj.dx = dx;
  obj.dy = dy;
  return 0;
}

int set_angular_vel (float dr, po_handle obj) {
  /*  // check input
  if (obj == NULL) {
    return 1;
    } */

  // update rotation, return success
  obj.dr = dr;
  return 0;

}

