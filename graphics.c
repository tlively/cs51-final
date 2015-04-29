#define SDL_MAIN_HANDLED

#include "graphics.h"
#include "SDL.h"


/* Graphics engine internal data types */

typedef struct texture {
  SDL_Texture* tex;
  SDL_Rect src_rect;
  int flip_h;
  int flip_v;
} texture;

typedef struct renderer {
  SDL_Window* win;
  SDL_Renderer* rend;
  void (*close_callback)(void);
} renderer;

/* External data types */
typedef renderer* renderer_handle;
typedef texture* texture_handle;


/* global state */
static int sdl_initialized = 0;


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

void init_graphics() {
  if (!sdl_initialized) {
    SDL_SetMainReady();
    SDL_Init(SDL_INIT_VIDEO);
  }
}

renderer_handle create_window(int width, int height, const char* title, int fullscreen) {
  renderer_handle rend = malloc(sizeof(renderer));
  rend->win = SDL_CreateWindow(title, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, width, height,
			       SDL_WINDOW_INPUT_GRABBED | ((fullscreen) ? SDL_WINDOW_FULLSCREEN_DESKTOP : 0));
  if (rend->win == NULL) {
    free(rend);
    return NULL;
  }
  rend->rend = SDL_CreateRenderer(rend->win, -1, 0);
  if (rend->rend == NULL) {
    SDL_DestroyWindow(rend->win);
    free(rend);
    return NULL;
  }

  rend->close_callback = NULL;

  return rend;
}

void destroy_window(renderer_handle win) {
  SDL_DestroyRenderer(win->rend);
  SDL_DestroyWindow(win->win);
  free(win);
}

int get_window_width(renderer_handle win) {
  int w, h;
  SDL_GetWindowSize(win->win, &w, &h);
  return w;
}

int get_window_height(renderer_handle win) {
  int w, h;
  SDL_GetWindowSize(win->win, &w, &h);
  return h;
}

/* The callback passed to SDL that will then call the user's callback */

void set_close_callback(renderer_handle rend, void (*callback)(void)) {
  rend->close_callback = callback;
}

void cleanup() {
  SDL_Quit();
}

texture_handle load_texture_data(renderer_handle renderer, int* pixels, int width, int height) {
  SDL_Surface* surface = SDL_CreateRGBSurfaceFrom(pixels, width, height,
						  32, width*sizeof(int),
						  0xFF000000, 0x00FF0000,
						  0x0000FF00, 0x000000FF);
  texture_handle tex = malloc(sizeof(texture));
  tex->tex = SDL_CreateTextureFromSurface(renderer->rend, surface);
  SDL_FreeSurface(surface);
  if (tex->tex == NULL) {
    free(tex);
    return NULL;
  }

  tex->src_rect = (SDL_Rect) {.x = 0, .y = 0, .w = width, .h = height};
  tex->flip_h = 0;
  tex->flip_v = 0;
  return tex;
}

void pack_textures() {
  return;
}

color get_clear_color(renderer_handle renderer) {
  Uint8 r,g,b,a;
  SDL_GetRenderDrawColor(renderer->rend, &r, &g, &b, &a);
  return get_color(r,g,b,a);
}

void set_clear_color(renderer_handle renderer, color c) {
  SDL_SetRenderDrawColor(renderer->rend, get_red(c), get_green(c),
			 get_blue(c), get_alpha(c));
}


/* Graphical functions */

void clear(renderer_handle renderer) {
  SDL_RenderClear(renderer->rend);
}


void draw(renderer_handle rend, texture_handle tex, int x, int y, double r,
	  int u, int v, int flip_h, int flip_v, double scale) {
  SDL_Rect dest = {.x = x, .y = y, .w = tex->src_rect.w * scale, .h = tex->src_rect.h * scale};
  SDL_Point center = {.x = u, .y = v};
  int flip = (flip_h && !tex->flip_h || !flip_h && tex->flip_h) ? SDL_FLIP_HORIZONTAL : 0;
  flip |= (flip_v && !tex->flip_v || !flip_v && tex->flip_v) ? SDL_FLIP_VERTICAL : 0;
  if (!flip) flip = SDL_FLIP_NONE;
  SDL_RenderCopyEx(rend->rend, tex->tex, &tex->src_rect, &dest, r * 180 / 3.1419526536, &center, flip);
}

void show(renderer_handle renderer) {
  SDL_RenderPresent(renderer->rend);
}
