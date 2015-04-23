#define SDL_MAIN_HANDLED

#include "graphics.h"
#include "SDL.h"


/* Graphics engine internal data types */

typedef struct texture {
  int placeholder;
} texture;

typedef struct renderer {
  int placeholder;
} renderer;

/* External data types */
typedef renderer* renderer_handle;
typedef texture* texture_handle;


/* Color functions */

int get_red(color c) { return (int)((c & 0xFF000000) >> 24); }
int get_green(color c) { return (int)((c & 0x00FF0000) >> 16); }
int get_blue(color c) { return (int)((c & 0x0000FF00) >> 8); }
int get_alpha(color c) { return (int)(c & 0x000000FF); }

color get_color(int r, int g, int b, int a) {
  return (color) (a | (b << 8) | (g << 16) | (r << 24));
}

color set_red(color c, int r) { return (c & 0x00FFFFFF) | (r << 24); }
color set_green(color c, int g) { return (c & 0xFF00FFFF) | (g << 16); }
color set_blue(color c, int b) { return (c & 0xFFFF00FF) | (b << 8); }
color set_alpha(color c, int a) { return (c & 0xFFFFFF00) | a; }

/* Core functionality */

renderer_handle init(int width, int height, int fullscreen) {
  return NULL;
}


texture_handle load_texture_data(int* pixels, int width, int height) {
  return NULL;
}

void pack_textures() {
  return;
}

color get_clear_color(renderer_handle renderer) {
  return 0;
}

void set_clear_color(renderer_handle renderer, color c) {
  return;
}


/* Graphical functions */

void clear(renderer_handle renderer) {
  return;
}


void draw(texture_handle tex, int x, int y, double r, int u, int v, int flip_h, int flip_v) {
  return;
}

void show() {
  return;
}
