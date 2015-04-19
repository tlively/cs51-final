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


int main() {

};
