/***********************************************************
 * Collisions.c
 * header for functions associated with colliding things into things
 ***********************************************************/

/* broadphase of collision detection: spatial hashing */
void coll_broadphase (world_handle world);

/* currently nonfunctional */
/* in the future, will be bounding boxes */
void coll_midphase(po_handle bucket1, po_handle bucket2);

/* currently for circles only */
/* future iterations will be parallel axis theorem */
void coll_narrowphase(po_handle obj1, po_handle obj2);

/* helper function for broadphase collision detection 
 * takes in a single row array and checks adjacent cells for contents */
void check_row(dynamic_array* row_k, int k_min, int k_max);

/* helper function for broadphase collision detection
 * takes in current row and next row and checks each cell in the top row 
 * for contents and contents in the bottom two buckets */
void check_rows(dynamic_array* row_k, int k_min, int k_max, 
		dynamic_array* row_kplus);
		
/* helper function for narrow phase that checks bounding box overlap
 * calls narrowphase */		
void check_bounding (po_handle obj1, po_handle obj2);

/* distance formula squared */
float distance_squared(po_vector point1, po_vector point2);
 
