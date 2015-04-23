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

  SDL_SetMainReady();
  if (SDL_Init(SDL_INIT_VIDEO|SDL_INIT_TIMER) != 0) {
    fprintf(stderr, "Failed to initialize SDL: %s\n", SDL_GetError());
    return 1;
  }
  atexit(SDL_Quit);

  SDL_Window* window = NULL;
  SDL_Renderer* renderer = NULL;
  SDL_Texture* texture = NULL;
  SDL_Surface* surface = NULL;


  int pixels[] = {0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff,
		  0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff,
		  0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff,
		  0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff,
		  0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff,
		  0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff, 0xff0000ff};
  
  int window_width = 800;
  int window_height = 600;
  if (!(window = SDL_CreateWindow("test application",
				  SDL_WINDOWPOS_CENTERED,
				  SDL_WINDOWPOS_CENTERED,
				  800, 600, SDL_WINDOW_RESIZABLE))) {
    fprintf(stderr, "Failed to create SDL window: %s\n", SDL_GetError());
    return 1;
  }

  if (!(renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED))) {
    fprintf(stderr, "Failed to create SDL renderer: %s\n", SDL_GetError());
    return 1;
  }

  if (!(surface = SDL_CreateRGBSurfaceFrom(&pixels, 6, 6, 32, 24, 0, 0, 0, 0))) {
    fprintf(stderr, "Failed to create surface: %s\n", SDL_GetError());
    return 1;
  }

  if (!(texture = SDL_CreateTextureFromSurface(renderer, surface))) {
    fprintf(stderr, "Failed to create texture: %s\n", SDL_GetError());
    return 1;
  }

  SDL_FreeSurface(surface);
  surface = NULL;

  if (SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0xff)) {
    fprintf(stderr, "Failed to set background color: %s\n", SDL_GetError());
    return 1;
  }


  while (SDL_GetTicks() < 10000) {
    SDL_RenderClear(renderer);
    double theta = (double) SDL_GetTicks() / 300;
    double r = 200 * cos(3*theta/5);

    SDL_Rect dst1 = {.x = r*cos(theta) - 3 + window_width/2,
		     .y = r*sin(theta) - 3 + window_height/2,
		     .w = 6,
		     .h = 6 };
    SDL_Rect dst2 = {.x = r*cos(theta + 2*PI/3) - 3 + window_width/2,
		     .y = r*sin(theta + 2*PI/3) - 3 + window_height/2,
		     .w = 6,
		     .h = 6 };
    SDL_Rect dst3 = {.x = r*cos(theta - 2*PI/3) - 3 + window_width/2,
		     .y = r*sin(theta - 2*PI/3) - 3 + window_height/2,
		     .w = 6,
		     .h = 6 };
    
    SDL_RenderCopy(renderer, texture, NULL, &dst1);
    SDL_RenderCopy(renderer, texture, NULL, &dst2);
    SDL_RenderCopy(renderer, texture, NULL, &dst3);
    SDL_RenderPresent(renderer);
    SDL_Delay(16);
  }

  SDL_DestroyTexture(texture);
  SDL_DestroyWindow(window);

  return 0;
}

