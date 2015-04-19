/******************************************************************************
 * Physics.h - the public interface to the physics engine
 * 
 * Takes physics object geometry data from the user and simulates rigid body
 * kinematics with a discrete timestep
 ******************************************************************************/

struct world_t;
struct po_t;

/* An opaque handle to a physics world. This allows multiple separate physics
 * worlds to be simulated simultaneously. */
typedef struct world_t world_handle;

/* An opaque handle to a physics object. The client keeps track of physics
 * objects, which are handled internally in the physics engine, using these
 * handles. */
typedef struct po_t po_handle;


/* A generic struct for handing two dimensional vectors. Used in many
 * situations with different meanings. */
typedef struct po_vector {
  float x;
  float y;
} po_vector;

/* define circle
 * origin dependent on context */
typedef struct po_circle {
    po_vector center;
    float radius;
    float density;
} po_circle;

/* define polygons
 * origin dependent on context */
typedef struct po_poly {
    po_vector* vertices;
    po_circle* circs;
    int npolys;
    int ncircs;
} po_poly;

/* define total geomentry of an object
 * origin defined in global coordinates */
typedef struct po_geometry {
    po_poly* polys;
    po_circle* circs;
    int npolys;
    int ncircs;
} po_geometry; 





