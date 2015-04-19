#include "physics.h"

struct world_t {
  int test_data;
};

struct po_implementation {
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
};


/* Updates object's global position based on velocity
 * Future versions may include more sophistocated algorthims using acceleration */
void integrate (float dx, float dy, float dr, float time_step, po_implementation my_object) {
  my_object.x = my_object.x + (dx * time_step);
  my_object.y = my_object.y + (dy * time_step);
  my_object.r = my_object.r + (dr * time_step);
}

int main() {

};
