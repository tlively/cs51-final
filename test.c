#include "SDL.h"
#include "physics.h"
#include "graphics.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define PI 3.1415926536

void test_colors() {
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


int main() {

  printf("testing colors\n");
  test_colors();

  int pixels[] = {0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff,
		  0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff,
		  0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff,
		  0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff,
		  0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff,
		  0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff};
  
  int window_width = 800;
  int window_height = 600;

  SDL_SetMainReady();

  renderer_handle rend = init(window_width, window_height, 0);
  texture_handle tex = load_texture_data(rend, pixels, 6, 6);
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

    draw(rend, tex, x1, y1, 0, 0, 0, 0, 0);
    draw(rend, tex, x2, y2, 0, 0, 0, 0, 0);
    draw(rend, tex, x3, y3, 0, 0, 0, 0, 0);
    show(rend);
    SDL_Delay(16);
  }

  //  SDL_DestroyTexture(texture);
  //  SDL_DestroyWindow(window);

  return 0;
}

