#include <stdlib.h>
#include <limits.h>
#include <string.h>

#define INITIAL_CAPACITY 16

typedef struct dynamic_array {
  int offset;
  int size;
  int min;
  int max;
  int capacity;
  void** arr;
} dynamic_array;

dynamic_array* dynamic_array_create() {
  dynamic_array* da = malloc(sizeof(dynamic_array));
  if (da == NULL) return NULL;
  da->offset = INITIAL_CAPACITY / 2;
  da->size = 0;
  da->min = 0;
  da->max = 0;
  da->capacity = INITIAL_CAPACITY;
  da->arr = malloc(INITIAL_CAPACITY * sizeof(void*));
  return da;
}

void dynamic_array_add(dynamic_array* da, int index, void* data) {
  if (da == NULL || data == NULL) return;
  // index out of bounds to the right
  if (index >= da->offset + da->capacity) {
    int new_capacity = da->capacity * 2;
    if (index >= da->offset + new_capacity) new_capacity = index - da->offset + 1;
    void** new_arr = malloc(new_capacity * sizeof(void*));
    memset(new_arr, (int)NULL, new_capacity * sizeof(void*));
    memcpy(new_arr, da->arr, da->capacity * sizeof(void*));
    free(da->arr);
    da->arr = new_arr;
    da->capacity = new_capacity;
  }
  // index out of bounds to the left
  else if (index < da->offset) {
    int new_capacity = da->capacity * 2;
    if (index < da->offset - da->capacity)
      new_capacity = da->capacity + da->offset - index;
    void** new_arr = malloc(new_capacity * sizeof(void*));
    memset(new_arr, (int)NULL, new_capacity * sizeof(void*));
    memcpy(&new_arr[new_capacity - da->capacity], da->arr,
	   da->capacity * sizeof(void*));
    free(da->arr);
    da->arr = new_arr;
    da->offset = da->offset + da->capacity - new_capacity;
    da->capacity = new_capacity;
  }
  // add item and update max and min
  da->arr[index - da->offset] = data;
  da->size++;
  if (index < da->min || da->size == 1) da->min = index;
  if (index > da->max || da->size == 1) da->max = index;
  return;
}

void* dynamic_array_get(dynamic_array* da, int index) {
  if (da == NULL || index < da->offset || index >= da->offset + da->capacity
      || da->size == 0)
    return NULL;
  return da->arr[index - da->offset];
}

void* dynamic_array_remove(dynamic_array* da, int index) {
  if (da == NULL || index < da->offset || index >= da->offset + da->capacity
      || da->size == 0)
    return NULL;
  void* ret = da->arr[index - da->offset];
  if (da->arr[index - da->offset]) da->size--;
  else return NULL;
  da->arr[index - da->offset] = NULL;
  if (da->size == 0) {
    da->min = 0;
    da->max = 0;
  }
  else {
    if (index == da->min) {
      for (int i = da->min + 1; i <= da->max; i++) {
	if (da->arr[i - da->offset] != NULL) {
	  da->min = i;
	  break;
	}
      }
    }
    if (index == da->max) {
      for (int i = da->max - 1; i >= da->min; i--) {
	if (da->arr[i - da->offset] != NULL) {
	  da->max = i;
	  break;
	}
      }
    }
  }
  return ret;
}

size_t dynamic_array_length(dynamic_array* da) {
  return (size_t) da->size;
}

int dynamic_array_min(dynamic_array* da){
  return da->min;
}

int dynamic_array_max(dynamic_array* da){
  return da->max;
}

void dynamic_array_shrink(dynamic_array* da) {
  // haha jokes
}

void dynamic_array_free(dynamic_array* da) {
  if(da != NULL) free(da->arr);
  free(da);
}
