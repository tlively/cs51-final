/****************************************************************************
 * graphics.h - the public interface to the graphics engine
 *
 * creates texture objects from the contents of provided PNG files, manages
 * windowing, and renders graphics
 ****************************************************************************/

#ifndef GRAPHICS_INCLUDED
#define GRAPHICS_INCLUDED

struct texture;
struct renderer;

/* An opaque handle to a texture. This allows the user to manipulate textures
 * using this API without worrying about their internal representation. */
typedef struct texture* texture_handle;

/* An opaque handle to a renderer. Usually this is a window, but it doesn't
 * need to be, since the supported rendering operations are independent of
 * rendering target. */
typedef struct renderer* renderer_handle;

/* A useful alias for color, allowing improved code clarity. */
typedef unsigned int color;

/*
// Sets up a window with the given dimensions
renderer_handle init(int width, int height, int fullscreen);
*/

/* Performs any necessary initialization of the graphics system. Must be called before
 * any other function in this module. */
void init_graphics(void);

/* Creates a new window object and returns a handle to its renderer */
renderer_handle create_window(int width, int height, const char* title, int fullscreen);

/* Destroys the window and its associated resources, including its renderer */
void destroy_window(renderer_handle win);

/* Cleans up any loose ends. Call this right before exiting the program */
void cleanup(void);

/* Creates a texture object with the given data. Pixel format is rrggbbaa. */
texture_handle load_texture_data(renderer_handle renderer, int* pixels, int width, int height);

/* Asks the graphics engine to pack textures internally to improve performance.
 * This function may be slow. The interface to the graphics engine is unaffected. */
void pack_textures(void);

/* Returns the current screen clearing color. */
color get_clear_color(renderer_handle renderer);

/* Sets the current screen clearing color. */
void set_clear_color(renderer_handle renderer, color c);

/* Clear the screen. Should be called at the beginning of each frame. */
void clear(renderer_handle renderer);

/* Draw a texture at a given location (x,y) in the renderer target with a given
 * rotation. The texture will be rotated about the given center point (u,v)
 * relative to the texture's local origin. The texture will be flipped vertically
 * if flip_v is true and horizontally if flip_h is true. */
void draw(renderer_handle rend, texture_handle tex, int x, int y, double r,
	  int u, int v, int flip_h, int flip_v, double scale);

/* Swap screen buffers to show the contents drawn since that last call to show().
 * Should be called at the end of each frame. */
void show(renderer_handle renderer);

/* Utility functions for dealing with colors */
int get_red(color c);
int get_green(color c);
int get_blue(color c);
int get_alpha(color c);
color get_color(int r, int g, int b, int a);
color set_red(color c, int r);
color set_green(color c, int g);
color set_blue(color c, int b);
color set_alpha(color c, int a);

#endif
