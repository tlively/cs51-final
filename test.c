#include "SDL.h"
#include "physics.h"
#include "graphics.h"
#include "dynamic_array.h"
#include "vector_math.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define PI 3.14159265358979323846

void test_colors() {
  printf("testing colors\n");
  int r = 25;
  int g = 225;
  int b = 200;
  int a = 255;
  color c = get_color(r,g,b,a);
  assert(c == 0x19E1C8FF);
  assert(get_red(c) == r);
  assert(get_green(c) == g);
  assert(get_blue(c) == b);
  assert(get_alpha(c) == a);
  assert(set_red(c,0xAA) == 0xAAE1C8FF);
  assert(set_green(c,0xAA) == 0x19AAC8FF);
  assert(set_blue(c,0xAA) == 0x19E1AAFF);
  assert(set_alpha(c,0xAA) == 0x19E1C8AA);
}
void test_vectors() {
  po_vector v;
  v.x = 5.0; v.y = 5.0;
  po_vector o;
  o.x = 0.0; o.y = 0.0;
  po_vector a;
  a.x = 3.0; a.y = 4.0;
  assert(vect_mag_squared(vect_from_points(v,v))
	 == 0.0);
  assert(vect_normal(v).x == 5.0
	 && vect_normal(v).y == -5.0);
  assert(vect_mag_squared(a) == 25.0);
  assert(vect_dot_prod(a,a) == 25.0);
  assert(vect_mag_squared(vect_project(v,o))
	 == vect_mag_squared(o));
  assert(vect_cross_prod(v,a) == 5.0);
}

void test_dynamic_array() {
  printf("testing dynamic arrays\n");
  dynamic_array* da = dynamic_array_create();
  assert(da != NULL);
  assert(dynamic_array_length(da) == 0);
  int a = 1;
  int b = 2;
  int c = 3;
  int d = 4;
  int e = 5;
  assert(dynamic_array_min(da) == 0);
  assert(dynamic_array_max(da) == 0);
  dynamic_array_add(da, 0, &a);
  dynamic_array_add(da, 12, &b);
  dynamic_array_add(da, -12, &c);
  dynamic_array_add(da, 50, &d);
  dynamic_array_add(da, -666, &e);
  assert(dynamic_array_min(da) == -666);
  assert(dynamic_array_max(da) == 50);
  assert(dynamic_array_length(da) == 5);
  assert(dynamic_array_remove(da, 0) == &a);
  assert(dynamic_array_remove(da, 0) == NULL);
  assert(dynamic_array_remove(da, 1) == NULL);
  assert(dynamic_array_remove(da, 10000) == NULL);
  assert(dynamic_array_remove(da, -999) == NULL);
  assert(dynamic_array_length(da) == 4);
  assert(dynamic_array_get(da, 12) == &b);
  assert(dynamic_array_get(da, 12) == &b);
  assert(dynamic_array_get(da, -12) == &c);
  assert(dynamic_array_get(da, 50) == &d);
  assert(dynamic_array_remove(da, -666) == &e);
  assert(dynamic_array_min(da) == -12);
  assert(dynamic_array_max(da) == 50);
  dynamic_array_remove(da, 12);
  dynamic_array_remove(da, -12);
  dynamic_array_remove(da, 50);
  assert(dynamic_array_length(da) == 0);
  assert(dynamic_array_min(da) == 0);
  assert(dynamic_array_max(da) == 0);
  dynamic_array_free(da);
  
}
void test_physics(){
  printf("testing physics\n");
  world_handle hello = new_world();
  po_vector center;
  center.x = 0.0;
  center.y = 0.0;
  po_circle circle1 = create_circ(center, 5.0, 1.0);
  po_geometry geo_circle = create_geom_circ(circle1);
  add_object (hello, &geo_circle, 0.0, 0.0, 0.0);
}
int main() {

  test_colors();
  test_dynamic_array();
  test_vectors();
  test_physics();

  int pixels[] = {0x0000ffff, 0x0000ffff, 0x0000ffff, 0x0000ffff, 0x0000ffff, 0x0000ffff,
		  0x0000ffff, 0x0000ffff, 0x0000ffff, 0x0000ffff, 0x0000ffff, 0x0000ffff,
		  0x0000ffff, 0x0000ffff, 0x0000ffff, 0x0000ffff, 0x0000ffff, 0x0000ffff,
		  0x0000ffff, 0x0000ffff, 0x0000ffff, 0x0000ffff, 0x0000ffff, 0x0000ffff,
		  0x0000ffff, 0x0000ffff, 0x0000ffff, 0x0000ffff, 0x0000ffff, 0x0000ffff,
		  0x0000ffff, 0x0000ffff, 0x0000ffff, 0x0000ffff, 0x0000ffff, 0x0000ffff};
  
  color pixels2[10000];
  for (int i = 0; i < 10000; i++) pixels2[i] = 0xE650E1FF;

  int window_width = 800;
  int window_height = 600;
  
  init_graphics();

  renderer_handle rend = create_window(window_width, window_height, "TEST", 0);
  texture_handle tex = load_texture_data(rend, pixels, 6, 6);
  texture_handle tex2 = load_texture_data(rend, pixels2, 100, 100);
  set_clear_color(rend, get_color(0,0,0,0xFF));

  while (SDL_GetTicks() < 10000) {
    clear(rend);
    double theta = (double) SDL_GetTicks() / 300;
    double r = 200 * cos(3*theta/5);

    int x1 = r*cos(theta) - 3 + window_width/2;
    int y1 = r*sin(theta) - 3 + window_height/2;
    int x2 = r*cos(theta + 2*PI/3) - 3 + window_width/2;
    int y2 = r*sin(theta + 2*PI/3) - 3 + window_height/2;
    int x3 = r*cos(theta - 2*PI/3) - 3 + window_width/2;
    int y3 = r*sin(theta - 2*PI/3) - 3 + window_height/2;

    double scale = .75 - .25 * cos(6*theta/5);
    int s = 100 * scale;

    draw(rend, tex2, (window_width - s) / 2, (window_height - s) / 2, -theta, s/2, s/2, 0, 0, scale);
    draw(rend, tex, x1, y1, 0, 0, 0, 0, 0, 1);
    draw(rend, tex, x2, y2, 0, 0, 0, 0, 0, 1);
    draw(rend, tex, x3, y3, 0, 0, 0, 0, 0, 1);
    show(rend);
    SDL_Delay(16);
  }

  destroy_window(rend);
  cleanup();
  
  return 0;
}

