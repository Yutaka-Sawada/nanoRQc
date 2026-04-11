#include "oblas.h"
#include "octmat.h"

void om_resize(octmat *v, unsigned rows, unsigned cols) {
	// allocate memory and move to the alignment
	uint8_t *aligned;
	unsigned cols_al = ALIGNED_COLS(cols);
	size_t aligned_sz = (size_t)rows * cols_al + OCTMAT_ALIGN;
	if (v->data - v->off_al == NULL){	// allocate and zero fill
		aligned = calloc(aligned_sz, 1);
	} else {	// re-allocate size and keep original data
		aligned = realloc(v->data - v->off_al, aligned_sz);
	}
	if (aligned == NULL)
		return;	// If it fails to allocate memory, it does nothing.

	int off_al = (size_t)aligned % OCTMAT_ALIGN;
	if (off_al != 0){
		off_al = OCTMAT_ALIGN - off_al;
		aligned += off_al;
	}
	v->data = aligned;
	v->rows = rows;
	v->cols = cols;
	v->cols_al = cols_al;
	v->off_al = off_al;	// save offset from the allocated address
}

void om_destroy(octmat *v) {
	v->rows = 0;
	v->cols = 0;
	v->cols_al = 0;
	if (v->data){
		free(v->data - v->off_al);	// restore original address
		v->data = NULL;
	}
	v->off_al = 0;
}
