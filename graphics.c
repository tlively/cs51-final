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

/* input data */
static int keys[127];/* = {-1,-1,-1,-1,-1,-1,-1,-1,-1,0,
                        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                        -1,-1,0,-1,-1,-1,-*/
static bool mouse_pressed;
static bool shift_pressed;
static int mouse_x;
static int mouse_y;

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
  
  SDL_SetWindowData(rend->win,"renderer",rend);

    /* setup for ascii key array */
    for(int a = 0; a < sizeof(a); a++) {
      if((a>96 && a<122) || (a>47 && a<58) || a == 9 || a == 32 || a == 127)
        keys[a] = 0;
      else
        keys[a] = -1;
    } 
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
  SDL_Event event;
  while(SDL_PollEvent(&event)){
    switch(event.type){
    case SDL_KEYDOWN:
      //toggle bool true for each key
      switch(event.key.keysym.scancode){
      case SDL_SCANCODE_0:
	keys[48] = 1;
        break;
      case SDL_SCANCODE_1:
	keys[49] = 1;
        break;
      case SDL_SCANCODE_2:
	keys[50] = 1;
        break;
      case SDL_SCANCODE_3:
	keys[51] = 1;
        break;
      case SDL_SCANCODE_4:
	keys[52] =1;
        break;
      case SDL_SCANCODE_5:
	keys[53] = 1;
        break;
      case SDL_SCANCODE_6:
	keys[54] = 1;
        break;
      case SDL_SCANCODE_7:
	keys[55] = 1;
        break;
      case SDL_SCANCODE_8:
	keys[56] = 1;
        break;
      case SDL_SCANCODE_9:
	keys[57] = 1;
        break;
      case SDL_SCANCODE_A:
	keys[97] = 1;
        break;
      case SDL_SCANCODE_B:
        keys[98] = 1;
        break;
      case SDL_SCANCODE_C:
	keys[99] = 1;
        break;
      case SDL_SCANCODE_D:
	keys[100] = 1;
        break;
      case SDL_SCANCODE_E:
	keys[101] = 1;
        break;
      case SDL_SCANCODE_F:
	keys[102] = 1;
        break;
      case SDL_SCANCODE_G:
	keys[103] = 1;
        break;
      case SDL_SCANCODE_H:
	keys[104] = 1;
        break;
      case SDL_SCANCODE_I:
	keys[105] = 1;
        break;
      case SDL_SCANCODE_J:
	keys[106] = 1;
        break;
      case SDL_SCANCODE_K:
	keys[107] = 1;
        break;
      case SDL_SCANCODE_L:
	keys[108] = 1;
        break;
      case SDL_SCANCODE_M:
	keys[109] = 1;
        break;
      case SDL_SCANCODE_N:
	keys[110] = 1;
        break;
      case SDL_SCANCODE_O:
	keys[111] = 1;
        break;
      case SDL_SCANCODE_P:
	keys[112] = 1;
        break;
      case SDL_SCANCODE_Q:
	keys[113] = 1;
        break;
      case SDL_SCANCODE_R:
	keys[114] = 1;
        break;
      case SDL_SCANCODE_S:
	keys[115] = 1;
        break;
      case SDL_SCANCODE_T:
	keys[116] = 1;
        break;
      case SDL_SCANCODE_U:
	keys[117] = 1;
        break;
      case SDL_SCANCODE_V:
	keys[118] = 1;
        break;
      case SDL_SCANCODE_W:
	keys[119] = 1;
        break;
      case SDL_SCANCODE_X:
	keys[120] = 1;
        break;
      case SDL_SCANCODE_Y:
	keys[121] = 1;
        break;
      case SDL_SCANCODE_Z:
	keys[122] = 1;
        break;
      case SDL_SCANCODE_TAB:
	keys[9] = 1;
	break;
      case SDL_SCANCODE_SPACE:
	keys[32] = 1;
	break;
      case SDL_SCANCODE_DELETE:
	keys[127] = 1;
	break;
      case SDL_SCANCODE_LSHIFT:
      case SDL_SCANCODE_RSHIFT:
	shift_pressed = true;
	break;
      }
      break;
    case SDL_KEYUP:
      //toggle bool true for each key
      switch(event.key.keysym.scancode){
      case SDL_SCANCODE_0:
	keys[48] = 0;
        break;
      case SDL_SCANCODE_1:
	keys[49] = 0;
        break;
      case SDL_SCANCODE_2:
	keys[50] = 0;
        break;
      case SDL_SCANCODE_3:
	keys[51] = 0;
        break;
      case SDL_SCANCODE_4:
	keys[52] = 0;
        break;
      case SDL_SCANCODE_5:
	keys[53] = 0;
        break;
      case SDL_SCANCODE_6:
	keys[54] = 0;
        break;
      case SDL_SCANCODE_7:
	keys[55] = 0;
        break;
      case SDL_SCANCODE_8:
	keys[56] = 0;
        break;
      case SDL_SCANCODE_9:
	keys[57] = 0;
        break;
      case SDL_SCANCODE_A:
	keys[97] = 0;
        break;
      case SDL_SCANCODE_B:
        keys[98] = 0;
        break;
      case SDL_SCANCODE_C:
	keys[99] = 0;
        break;
      case SDL_SCANCODE_D:
	keys[100] = 0;
        break;
      case SDL_SCANCODE_E:
	keys[101] = 0;
        break;
      case SDL_SCANCODE_F:
	keys[102] = 0;
        break;
      case SDL_SCANCODE_G:
	keys[103] = 0;
        break;
      case SDL_SCANCODE_H:
	keys[104] = 0;
        break;
      case SDL_SCANCODE_I:
	keys[105] = 0;
        break;
      case SDL_SCANCODE_J:
	keys[106] = 0;
        break;
      case SDL_SCANCODE_K:
	keys[107] = 0;
        break;
      case SDL_SCANCODE_L:
	keys[108] = 0;
        break;
      case SDL_SCANCODE_M:
	keys[109] = 0;
        break;
      case SDL_SCANCODE_N:
	keys[110] = 0;
        break;
      case SDL_SCANCODE_O:
	keys[111] = 0;
        break;
      case SDL_SCANCODE_P:
	keys[112] = 0;
        break;
      case SDL_SCANCODE_Q:
	keys[113] = 0;
        break;
      case SDL_SCANCODE_R:
	keys[114] = 0;
        break;
      case SDL_SCANCODE_S:
	keys[115] = 0;
        break;
      case SDL_SCANCODE_T:
	keys[116] = 0;
        break;
      case SDL_SCANCODE_U:
	keys[117] = 0;
        break;
      case SDL_SCANCODE_V:
	keys[118] = 0;
        break;
      case SDL_SCANCODE_W:
	keys[119] = 0;
        break;
      case SDL_SCANCODE_X:
	keys[120] = 0;
        break;
      case SDL_SCANCODE_Y:
	keys[121] = 0;
        break;
      case SDL_SCANCODE_Z:
	keys[122] = 0;
        break;
      case SDL_SCANCODE_TAB:
	keys[9] = 0;
	break;
      case SDL_SCANCODE_SPACE:
	keys[32] = 0;
	break;
      case SDL_SCANCODE_DELETE:
	keys[127] = 0;
	break;
      case SDL_SCANCODE_LSHIFT:
      case SDL_SCANCODE_RSHIFT:
	shift_pressed = false;
	break;
      }
      break;
    case SDL_WINDOWEVENT:
      if(event.window.event == SDL_WINDOWEVENT_CLOSE){
        SDL_Window* window = SDL_GetWindowFromID(event.window.windowID);
	  renderer_handle rend = SDL_GetWindowData(window,"renderer");
      rend->close_callback();
    }
      break;
    case SDL_MOUSEMOTION:
      //store mouse position
      mouse_x = event.motion.x;
      mouse_y = event.motion.y;
      break;
    case SDL_MOUSEBUTTONDOWN:
      //toggle mouse button clicked
      mouse_pressed = true;
      break;
    case SDL_MOUSEBUTTONUP:
      //toggle mouse button not clicked
      mouse_pressed = false;
      break;
    }
  }
}

bool shift_is_pressed() {
  return shift_pressed;
}

bool mouse_is_pressed() {
  return mouse_pressed;
}

//returns -1 for not handled, 0 for up, 1 for down
int key_status_for_ascii(int index) {
  if(index < 0 || index > 127)
    return -1;
  return keys[index];
}

//mouse position in place
void get_mouse_position(int* x, int* y) {
  *x = mouse_x;
  *y = mouse_y;
}
