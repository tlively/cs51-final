/***************************************************************************
 * dynamic_array.h - solving all of your dynamic indexing needs since 2015 *
 ***************************************************************************/
#ifndef DYNAMIC_ARRAY_INCL
#define DYNAMIC_ARRAY_INCL

struct dynamic_array;
typedef struct dynamic_array dynamic_array;

// set up a newly created dynamic array with all of the appropriate internals.
dynamic_array* dynamic_array_create(void);

// add the data pointer to the array at index "index". Indices can be negative.
void dynamic_array_add(dynamic_array* da, int index, void* data);

// returns the item at the given index or NULL if it doesn't exist.
void* dynamic_array_get(dynamic_array* da, int index);

// same as dynamic_array_get but also removes the item as a side effect.
void* dynamic_array_remove(dynamic_array* da, int index);

// reduce total memory usage by removing padding on the edges of da.
void dynamic_array_shrink(dynamic_array* da);

// safely and cleanly frees all memory associated with da.
void dynamic_array_free(dynamic_array* da);

#endif
