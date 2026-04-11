#ifndef OCTMAT_H
#define OCTMAT_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef OBLAS_AVX
#define OCTMAT_ALIGN 32
#else
#ifdef OBLAS_SSE
#define OCTMAT_ALIGN 16
#endif
#endif

#ifndef OCTMAT_ALIGN
#define OCTMAT_ALIGN 8
#endif

typedef struct {
  uint8_t *data;
  unsigned rows;
  unsigned cols;
  unsigned cols_al;
  int off_al;
} octmat;

#define OM_INITIAL { .rows = 0, .cols = 0, .cols_al = 0, .data = 0, .off_al = 0 }

// When T (cols or cols_al) is larger than 65536 bytes, x * v may become over than 32-bit.
// Then, cast x with size_t for 64-bit integer multiplication.
#define om_R(v, x) ((v).data + ((size_t)(x) * (v).cols_al))
#define om_P(v) ((v).data)
#define om_A(v, x, y) (om_R(v, x)[(y)])

void om_resize(octmat *v, unsigned rows, unsigned cols);
void om_destroy(octmat *v);

#endif
