/****************************************************************************
 * Physics.h - the public interface to the physics engine
 * 
 * Takes physics object geometry data from the user and simulates rigid body
 * kinematics with a discrete timestep
 ******************************************************************************/
#ifndef PHYSICS_INCLUDED
#define PHYSICS_INCLUDED

#include "vector_math.h"

struct world_t;
struct po_imp;

/* An opaque handle to a physics world. This allows multiple separate physics
 * worlds to be simulated simultaneously. */
typedef struct world_t* world_handle;

/* An opaque handle to a physics object. The client keeps track of physics
 * objects, which are handled internally in the physics engine, using these
 * handles. */
typedef struct po_imp* po_handle;

/* define circle
 * origin dependent on context */
typedef struct po_circle {
    po_vector center;
    float radius;
    float density;
} po_circle;

/* define polygons
 * origin dependent on context */
// TDOD: ask about origin of po_polys
typedef struct po_poly {
  po_vector* vertices;
  int nvert;

} po_poly;

/* define total geomentry of an object
 * origin defined in global coordinates */
typedef struct po_geometry {
  // 0 if circle, 1 if polygon
  int shape_type;
  
  // one of these will be ignored by the program depending on shape type
  po_poly poly;
  po_circle circ;

} po_geometry; 

/* create a new world */
world_handle new_world ();

/* Places a physics object in the given world
 * origin set at global coordinates (x,y)
 * rotation is defined as rotation in radians 
 * returns handle for new phys obj or NULL on failure */
po_handle add_object (world_handle world, po_geometry* geom, 
		      float x, float y, float r);

/* flags an phys object for removal
 * object stops participating in collsions 
 * returns 0 on success, 1 on failure */
int remove_object (world_handle world, po_handle obj);

/* simulates dt seconds of physics in the given physics world
 * detects and resolves collisions */
int update (world_handle world, float dt);

/* attaches callback force to physics object
 * called when force applied > than min
 * returns 0 on success, 1 on failure */
int set_force_callback (po_handle obj, float min, 
			void (*callback)(po_handle, po_vector, po_vector));

/* allows user to set the global location of the object's origin 
 * returns 0 on success, 1 on failure*/
int set_location (po_handle obj, float x, float y);

/* allows user to set the initial rotation of the object 
 * returns 0 on success, 1 on failure */
int set_rotation (po_handle obj, float r);

/* allows user to set the initial velocity of the object 
 * returns 0 on success, 1 on failure */
int set_velocity (po_handle obj, float dx, float dy);

/* allows user to set the initial angular velocity of the object 
 * returns 0 on success, 1 on failure */
int set_angular_vel (po_handle obj, float dr);

/* two objects passed will be resolved when collison is detected (for now reversal of velocity) */ 
int resolve_collision (po_handle obj1, po_handle obj2);

#endif
