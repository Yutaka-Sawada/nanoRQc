#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nanorqc.h"
#include "nanorq_core.h"
#include "nanorq_ops.h"
#include "partition.h"
#include "util.h"

#define LOAD_OK (1 << 25)	// This value is larger than max ESI.
#define INVERT_OK (1 << 25)	// This value is larger than max ESI.

typedef struct {
  uint32_t esi;
  uint8_t *row;
} repair_sym;

typedef struct {
  repair_sym *a;
  size_t n;
  size_t m;
} repair_vec;

typedef struct {
  uint32_t *words;
  size_t size;
} compat_bitmask;

struct block_encoder {
  uint16_t K; /* the number of source symbols */
  int loaded;
  bool inverted;

  // Core low level solver:
  nanorq_core core;

  // Dynamic scratch memory for core solver:
  uint8_t *prep_mem;
  size_t prep_len;
  uint8_t *work_mem;
  size_t work_len;

  size_t matrix_rows; /* allocated capacity retained across resets */
  uint8_t *D;         /* matrix data */
  uint32_t stride;    /* recommended row stride */
  uint8_t *mix_buf;

  schedule S; // Operations schedule

  repair_vec repair_bin;
  compat_bitmask repair_mask;
};

struct nanorq {
  uint32_t T; /* the symbol size in octets */
  params P;
  uint32_t max_esi;

  struct block_encoder encoder[1]; // This supports single block only.

  // Cached precalculated components
  nanorq_core *precalc_core;
  uint8_t *precalc_prep_mem;
  uint8_t *precalc_work_mem;
  schedule *precalc_S;
};

static bool repair_vec_push(repair_vec *v, repair_sym val) {
  if (v->n >= v->m) {
    if (v->m > SIZE_MAX / 2)
      return false;
    size_t new_m = v->m == 0 ? 8 : v->m * 2;
    if (new_m > SIZE_MAX / sizeof(repair_sym))
      return false;
    repair_sym *new_a = (repair_sym *)realloc(v->a, new_m * sizeof(repair_sym));
    if (!new_a)
      return false;
    v->a = new_a;
    v->m = new_m;
  }
  v->a[v->n++] = val;
  return true;
}

static void repair_vec_free(repair_vec *v) {
  if (v->a) {
    for (size_t i = 0; i < v->n; i++) {
      obl_free(v->a[i].row);
    }
    free(v->a);
  }
  v->a = NULL;
  v->n = v->m = 0;
}

static compat_bitmask compat_bitmask_new(size_t size) {
  compat_bitmask b = {0};
  if (size > SIZE_MAX - 31)
    return b;
  size_t num_words = (size + 31) / 32;
  if (num_words > SIZE_MAX / sizeof(uint32_t))
    return b;
  b.words = (uint32_t *)calloc(num_words, sizeof(uint32_t));
  if (b.words)
    b.size = size;
  return b;
}

static void compat_bitmask_free(compat_bitmask *b) {
  free(b->words);
  b->words = NULL;
  b->size = 0;
}

static void compat_bitmask_set(compat_bitmask *b, size_t bit) {
  if (bit < b->size) {
    b->words[bit / 32] |= (1U << (bit % 32));
  }
}

static bool compat_bitmask_check(compat_bitmask *b, size_t bit) {
  if (bit < b->size) {
    return (b->words[bit / 32] & (1U << (bit % 32))) != 0;
  }
  return false;
}

static void compat_bitmask_reset(compat_bitmask *b) {
  size_t num_words = (b->size + 31) / 32;
  if (num_words > 0 && b->words)
    memset(b->words, 0, num_words * sizeof(uint32_t));
}

static size_t compat_bitmask_gaps(compat_bitmask *b, size_t limit) {
  size_t gaps = 0;
  size_t check_limit = limit < b->size ? limit : b->size;
  size_t limit_words = check_limit / 32;
  size_t limit_rem = check_limit % 32;

  for (size_t i = 0; i < limit_words; i++) {
    gaps += 32 - __builtin_popcount(b->words[i]);
  }

  if (limit_rem > 0) {
    uint32_t mask = (1U << limit_rem) - 1;
    uint32_t last_word = b->words[limit_words] & mask;
    gaps += limit_rem - __builtin_popcount(last_word);
  }

  if (limit > b->size) {
    gaps += limit - b->size;
  }

  return gaps;
}

static void transfer_esi(uint32_t esi, uint8_t *ptr, size_t len, uint8_t *src, int out) {
  size_t offset = esi * len;
  if (out){	// write on src
    memcpy(src + offset, ptr, len);
  } else {	// read from src
    memcpy(ptr, src + offset, len);
  }
}

static bool set_block_encoder(nanorq *rq, uint16_t K) {
  if (!rq)
    return false;

  struct block_encoder *enc = rq->encoder;
  enc->K = K;

  enc->repair_mask = compat_bitmask_new(K);
  if (!enc->repair_mask.words) {
    return false;
  }
  /* pack matrix to SIMD alignment */
  uint32_t align = get_align_size();
  enc->stride = ((rq->T + align - 1) / align) * align;
  if (enc->stride == 0) {
    compat_bitmask_free(&enc->repair_mask);
    return false;
  }

  return true;
}

/*
 K: number of source symbols
 T: size of each symbol in bytes
*/
nanorq *nanorq_encoder_new(uint16_t K, uint32_t T) {
  if (K == 0 || K > K_max || T < 1)
    return NULL;

  nanorq_core_init();

  nanorq *rq = (nanorq *)calloc(1, sizeof(nanorq));
  if (!rq)
    return NULL;
  rq->T = T;
  rq->P = params_init(K);

  rq->max_esi = (1 << 24) - 1;	// Is this required for encoder ?

  if (!set_block_encoder(rq, K)){
    free(rq);
    return NULL;
  }

  return rq;
}

void nanorq_free(nanorq *rq) {
  if (!rq)
    return;
  nanorq_encoder_cleanup(rq);

  if (rq->precalc_core) {
    if (rq->precalc_prep_mem)
      free(rq->precalc_prep_mem);
    if (rq->precalc_work_mem)
      free(rq->precalc_work_mem);
    if (rq->precalc_S) {
      if (rq->precalc_S->ops.a)
        free(rq->precalc_S->ops.a);
      free(rq->precalc_S);
    }
    free(rq->precalc_core);
  }

  free(rq);
}

/*
 K: number of source symbols, which MUST be same as K in encoder
 T: size of each symbol in bytes, which MUST be same as T in encoder
*/
nanorq *nanorq_decoder_new(uint16_t K, uint32_t T) {
  if (K == 0 || K > K_max || T < 1)
    return NULL;

  nanorq_core_init();

  nanorq *rq = (nanorq *)calloc(1, sizeof(nanorq));
  if (!rq)
    return NULL;
  rq->T = T;
  rq->P = params_init(K);

  rq->max_esi = (1 << 24) - 1;

  if (!set_block_encoder(rq, K)){
    free(rq);
    return NULL;
  }

  return rq;
}

int nanorq_precalculate(nanorq *rq) {
  if (!rq)
    return NANORQ_SYM_PARA;
  if (rq->precalc_core)
    return NANORQ_SYM_OK; // Already precalculated

  uint16_t K = rq->encoder->K;
  nanorq_core *core = (nanorq_core *)calloc(1, sizeof(nanorq_core));
  uint8_t *prep_mem = NULL;
  uint8_t *work_mem = NULL;
  schedule *S = NULL;
  void *sched_mem = NULL;
  if (!core || !nanorq_core_encoder_new(K, 0, core))
    goto fail;

  size_t prep_len = nanorq_core_calculate_prepare_memory(core);
  prep_mem = (uint8_t *)malloc(prep_len);
  if (!prep_mem || !nanorq_core_prepare(core, prep_mem, prep_len))
    goto fail;

  size_t work_len = nanorq_core_calculate_prepared_work_memory(core);
  work_mem = (uint8_t *)malloc(work_len);
  S = (schedule *)calloc(1, sizeof(schedule));
  size_t sched_bytes = ops_estimate_schedule_bytes(K);
  sched_mem = malloc(sched_bytes);
  if (!work_mem || !S || !schedule_init(S, sched_mem, sched_bytes))
    goto fail;

  nanorq_core_set_op_callback(core, S, ops_push);
  if ((nanorq_core_precalculate(core, work_mem, work_len) <= 0) ||
      S->overflowed || S->cpidx != 2)
    goto fail;

  rq->precalc_core = core;
  rq->precalc_prep_mem = prep_mem;
  rq->precalc_work_mem = work_mem;
  rq->precalc_S = S;
  return NANORQ_SYM_OK;

fail:
  free(sched_mem);
  free(S);
  free(work_mem);
  free(prep_mem);
  free(core);
  return NANORQ_SYM_ERR;
}

/* retain scratch storage across block resets */
static bool reserve_scratch(uint8_t **mem, size_t *capacity, size_t bytes) {
  if (bytes == 0)
    return false;
  if (*mem && *capacity >= bytes)
    return true;
  uint8_t *next = (uint8_t *)realloc(*mem, bytes);
  if (!next)
    return false;
  *mem = next;
  *capacity = bytes;
  return true;
}

static bool reset_schedule(schedule *S, uint16_t K) {
  size_t bytes = ops_estimate_schedule_bytes(K);
  if (bytes == 0)
    return false;
  if (!S->ops.a) {
    S->ops.a = (sched_op *)malloc(bytes);
    if (!S->ops.a)
      return false;
  }
  return schedule_init(S, S->ops.a, bytes);
}

static bool ensure_block_matrix(struct block_encoder *enc) {
  if (enc->D)
    return true;
  if (!nanorq_core_encoder_new(enc->K, 0, &enc->core))
    return false;
  enc->D = (uint8_t *)obl_alloc(nanorq_core_get_pc_rows(&enc->core),
                                enc->stride, nanorq_oblas.align_size);
  if (!enc->D)
    return false;
  enc->matrix_rows = nanorq_core_get_pc_rows(&enc->core);
  return true;
}

static bool ensure_mix_buf(struct block_encoder *enc) {
  if (!enc->mix_buf)
    enc->mix_buf = (uint8_t *)malloc(enc->stride);
  return enc->mix_buf != NULL;
}

/*
 data: input bytes of adding source symbols (size must be esi_count * T)
 esi_start: the first Encoding Symbol ID (range is 0 ~ K-1)
 esi_count: number of adding source symbols (range is 1 ~ K)
*/
int nanorq_encoder_add_symbols(nanorq *rq, void *data, int esi_start, int esi_count) {
  if (!rq || !data)
    return NANORQ_SYM_PARA;
  struct block_encoder *enc = rq->encoder;
  if (esi_start < 0 || esi_start >= enc->K || esi_start + esi_count > enc->K)
    return NANORQ_SYM_PARA;	// ESI is out of range.

  if (enc->loaded >= enc->K)
    return NANORQ_SYM_IGN;	// All source symbols were added already.

  if (!ensure_block_matrix(enc))
    return NANORQ_SYM_ERR;

  int success = NANORQ_SYM_OK;
  uint8_t *src = data;
  for (int esi = esi_start; esi < esi_start + esi_count; esi++) {
    if (compat_bitmask_check(&enc->repair_mask, esi)){
      success = NANORQ_SYM_DUP;
    } else {
      uint8_t *row = nanorq_core_get_symbol_ptr(&enc->core, enc->D, enc->stride, esi);
      memcpy(row, src, rq->T);
      compat_bitmask_set(&enc->repair_mask, esi);
      enc->loaded += 1;	// count how many symbols were added.
    }
    src += rq->T;
  }

  return success;
}

/*
 src: input bytes of all source symbols (total size must be K * T)
*/
int nanorq_generate_symbols(nanorq *rq, void *src) {
  if (!rq)
    return NANORQ_SYM_PARA;
  struct block_encoder *enc = rq->encoder;

  if (enc->inverted)
    return NANORQ_SYM_OK;	// Intermediate symbols were generated already.

  if (enc->loaded < enc->K) {
    if (!src)
      return NANORQ_SYM_PARA;
    if (!ensure_block_matrix(enc))
      return NANORQ_SYM_ERR;

    for (int esi = 0; esi < enc->K; esi++) {
      uint8_t *row = nanorq_core_get_symbol_ptr(&enc->core, enc->D, enc->stride, esi);
      transfer_esi(esi, row, rq->T, src, 0);
    }
    enc->loaded = enc->K;	// All source symbols were added.
  }

  if (rq->precalc_core) {
    ops_run(rq->precalc_core, enc->D, enc->stride, rq->precalc_S);
    enc->inverted = true;
    return NANORQ_SYM_OK;
  }

  if (!nanorq_core_encoder_new(enc->K, 0, &enc->core))
    return false;
  size_t prep_len = nanorq_core_calculate_prepare_memory(&enc->core);
  if (!reserve_scratch(&enc->prep_mem, &enc->prep_len, prep_len) ||
      !nanorq_core_prepare(&enc->core, enc->prep_mem, enc->prep_len))
    return NANORQ_SYM_ERR;

  size_t work_len = nanorq_core_calculate_prepared_work_memory(&enc->core);
  if (!reserve_scratch(&enc->work_mem, &enc->work_len, work_len) ||
      !reset_schedule(&enc->S, enc->K))
    return NANORQ_SYM_ERR;

  nanorq_core_set_op_callback(&enc->core, &enc->S, ops_push);
  if ((nanorq_core_precalculate(&enc->core, enc->work_mem, enc->work_len) <= 0) ||
      enc->S.overflowed || enc->S.cpidx != 2)
    return NANORQ_SYM_ERR;

  ops_run(&enc->core, enc->D, enc->stride, &enc->S);
  enc->inverted = true;

  return NANORQ_SYM_OK;
}

/*
 data: output symbol bytes for the ESI (size must be T)
 esi: Encoding Symbol ID, 0 ~ K-1 for source symbols, K ~ 1<<24-1 for repair symbols
      When esi < K, source symbol is restored like nanorq_restore().
      When esi >= K, repair symbol is created.
 src: input bytes of all source symbols (total size must be K * T)
      If all source symbols were added or intermediate symbols were generated, src can be NULL.
*/
int nanorq_encode(nanorq *rq, void *data, uint32_t esi, void *src) {
  if (!rq || !data)
    return NANORQ_SYM_PARA;
  struct block_encoder *enc = rq->encoder;

  if (esi < enc->K && !enc->inverted) {
    if (!src)
      return NANORQ_SYM_PARA;
    transfer_esi(esi, (uint8_t *)data, rq->T, src, 0);
    return NANORQ_SYM_OK;
  }

  if (esi > ((1 << 24) - 1))
    return NANORQ_SYM_PARA;
  if (!enc->inverted) {
    int success = nanorq_generate_symbols(rq, src);
    if (success != 0)
      return success;
  }
  if (!ensure_mix_buf(enc))
    return NANORQ_SYM_ERR;
  ops_mix(&enc->core, enc->D, enc->stride, esi, enc->mix_buf);
  memcpy(data, enc->mix_buf, rq->T);
  return NANORQ_SYM_OK;
}

void nanorq_encoder_cleanup(nanorq *rq) {
  if (!rq)
    return;
  struct block_encoder *enc = rq->encoder;
  if (enc->D) {
    obl_free(enc->D);
  }
  free(enc->mix_buf);
  if (enc->prep_mem) {
    free(enc->prep_mem);
  }
  if (enc->work_mem) {
    free(enc->work_mem);
  }
  if (enc->S.ops.a) {
    free(enc->S.ops.a);
  }
  repair_vec_free(&enc->repair_bin);
  compat_bitmask_free(&enc->repair_mask);
}

void nanorq_encoder_reset(nanorq *rq) {
  if (!rq)
    return;
  struct block_encoder *enc = rq->encoder;
  enc->loaded = 0;
  enc->inverted = false;
  if (enc->core.overhead != 0) {
    nanorq_core base_core = {0};
    if (nanorq_core_encoder_new(enc->K, 0, &base_core))
      enc->core = base_core;
  }
  if (enc->D) {
    nanorq_core_init_matrix(&enc->core, enc->D, enc->stride);
  }
  repair_vec_free(&enc->repair_bin);
  compat_bitmask_reset(&enc->repair_mask);
}

int nanorq_set_max_esi(nanorq *rq, uint32_t max_esi) {
  if (!rq || max_esi >= (1 << 24) || max_esi < rq->P.Kprime)
    return NANORQ_SYM_PARA;
  rq->max_esi = max_esi;
  return NANORQ_SYM_OK;
}

/*
 data: input symbol bytes for the ESI (size must be T)
 esi: Encoding Symbol ID, 0 ~ K-1 for source symbols, K ~ 1<<24-1 for repair symbols
 src: output bytes of all source symbols (total size must be K * T)
      If output isn't needed, set src to NULL. Or else, speed may become slower.
*/
int nanorq_decoder_add_symbol(nanorq *rq, void *data, uint32_t esi, void *src) {
  if (!rq || !data)
    return NANORQ_SYM_PARA;
  if (esi >= (1 << 24) || esi > rq->max_esi)
    return NANORQ_SYM_PARA;
  struct block_encoder *dec = rq->encoder;

  if (compat_bitmask_gaps(&dec->repair_mask, dec->K) == 0) {
    return NANORQ_SYM_IGN;	// Because all source symbols exist, repair isn't needed.
  }

  if (esi < dec->K) {
    if (compat_bitmask_check(&dec->repair_mask, esi))
      return NANORQ_SYM_DUP;
  } else {
    for (size_t i = 0; i < dec->repair_bin.n; i++) {
      if (dec->repair_bin.a[i].esi == esi) {
        return NANORQ_SYM_DUP;
      }
    }
  }

  if (!ensure_block_matrix(dec))
    return NANORQ_SYM_ERR;

  if (esi < dec->K) {
    if (src != NULL)
      transfer_esi(esi, (uint8_t *)data, rq->T, src, 1);
    if (!dec->inverted)
      nanorq_core_place_symbol(&dec->core, dec->D, dec->stride, esi,
                               (const uint8_t *)data, rq->T);
    compat_bitmask_set(&dec->repair_mask, esi);
  } else {
    repair_sym rs = {0};
    rs.esi = esi;
    rs.row = (uint8_t *)obl_alloc(1, dec->stride, nanorq_oblas.align_size);
    if (!rs.row)
      return NANORQ_SYM_ERR;
    memcpy(rs.row, data, rq->T);
    if (!repair_vec_push(&dec->repair_bin, rs)) {
      obl_free(rs.row);
      return NANORQ_SYM_ERR;
    }
  }
  dec->loaded += 1;	// count how many symbols were added.

  return NANORQ_SYM_OK;
}

int nanorq_num_loaded(nanorq *rq) {
  struct block_encoder *dec = rq->encoder;
  return dec->loaded;
}

int nanorq_num_missing(nanorq *rq) {
  struct block_encoder *dec = rq->encoder;
  return (int)compat_bitmask_gaps(&dec->repair_mask, dec->K);
}

int nanorq_num_repair(nanorq *rq) {
  struct block_encoder *dec = rq->encoder;
  return (int)dec->repair_bin.n;
}

static void write_recovered_block(nanorq *rq, struct block_encoder *dec,
                                  void *src) {
  for (uint32_t gap = 0; gap < dec->K; gap++) {
    if (compat_bitmask_check(&dec->repair_mask, gap))
      continue;
    ops_mix(&dec->core, dec->D, dec->stride, gap, dec->mix_buf);
    transfer_esi(gap, dec->mix_buf, rq->T, src, 1);
  }

  for (uint32_t gap = 0; gap < dec->K; gap++) {
    if (!compat_bitmask_check(&dec->repair_mask, gap))
      compat_bitmask_set(&dec->repair_mask, gap);
  }
}

/*
 src: output bytes of all source symbols (total size must be K * T)
      Missing source symbols will be restored in their positions.
*/
int nanorq_repair(nanorq *rq, void *src) {
  if (!rq || !src)
    return NANORQ_SYM_PARA;
  struct block_encoder *dec = rq->encoder;

  size_t conventional_count = dec->repair_bin.n;

  size_t num_gaps = compat_bitmask_gaps(&dec->repair_mask, dec->K);
  if (num_gaps == 0)
    return NANORQ_SYM_IGN;	// Because all source symbols exist, repair isn't needed.
  if (dec->inverted) {	// When intermediate symbols were restored, it recovers all source symbols.
    if (!ensure_mix_buf(dec))
      return NANORQ_SYM_ERR;
    write_recovered_block(rq, dec, src);
    return NANORQ_SYM_OK;
  }
  if (conventional_count < num_gaps || conventional_count - num_gaps > UINT32_MAX)
    return NANORQ_SYM_LACK;	// added symbols aren't enough for recovery.

  uint32_t overhead = (uint32_t)(conventional_count - num_gaps);
  nanorq_core core = {0};
  bool matrix_modified = false;
  int success = NANORQ_SYM_ERR;

  if (!ensure_block_matrix(dec) ||
      !nanorq_core_encoder_new(dec->K, overhead, &core))
    goto out;

  size_t rows = nanorq_core_get_pc_rows(&core);
  size_t base_rows = nanorq_core_get_pc_rows(&dec->core);
  if (base_rows > rows)
    goto out;
  if (rows > dec->matrix_rows) {
    uint8_t *D =
        (uint8_t *)obl_alloc(rows, dec->stride, nanorq_oblas.align_size);
    if (!D)
      goto out;
    memcpy(D, dec->D, base_rows * dec->stride);
    obl_free(dec->D);
    dec->D = D;
    dec->matrix_rows = rows;
  } else {
    /* reset clears base matrix; clear only active overhead */
    memset(dec->D + base_rows * dec->stride, 0,
           (rows - base_rows) * dec->stride);
  }

  size_t prep_len = nanorq_core_calculate_prepare_memory(&core);
  if (!reserve_scratch(&dec->prep_mem, &dec->prep_len, prep_len) ||
      !nanorq_core_generate(&core, dec->prep_mem, dec->prep_len))
    goto out;

  size_t rep_idx = 0;
  matrix_modified = true;
  for (uint32_t gap = 0; gap < dec->K; gap++) {
    if (compat_bitmask_check(&dec->repair_mask, gap))
      continue;
    if (rep_idx == dec->repair_bin.n)
      goto out;
    repair_sym *rs = &dec->repair_bin.a[rep_idx++];
    nanorq_core_replace_symbol(&core, gap, rs->esi);
    nanorq_core_place_symbol(&core, dec->D, dec->stride, gap, rs->row,
                             rq->T);
  }

  for (uint32_t extra = 0; extra < overhead; extra++) {
    if (rep_idx == dec->repair_bin.n)
      goto out;
    repair_sym *rs = &dec->repair_bin.a[rep_idx++];
    uint32_t row = core.P.Kprime + extra;
    nanorq_core_replace_symbol(&core, row, rs->esi);
    nanorq_core_place_symbol(&core, dec->D, dec->stride, row, rs->row, rq->T);
  }

  if (!nanorq_core_patch_matrix(&core))
    goto out;
  size_t work_len = nanorq_core_calculate_prepared_work_memory(&core);
  if (!ensure_mix_buf(dec) ||
      !reserve_scratch(&dec->work_mem, &dec->work_len, work_len) ||
      !reset_schedule(&dec->S, dec->K))
    goto out;
  nanorq_core_set_op_callback(&core, &dec->S, ops_push);
  int invert = nanorq_core_precalculate(&core, dec->work_mem, dec->work_len);
  if ((invert <= 0) || dec->S.overflowed || dec->S.cpidx != 2){
    if (invert == 0)
      success = NANORQ_SYM_LACK;	// Failed to invert matrix. It may need more symbols.
    goto out;
  }

  dec->core = core;
  ops_run(&dec->core, dec->D, dec->stride, &dec->S);
  dec->inverted = true;
  write_recovered_block(rq, dec, src);
  success = NANORQ_SYM_OK;

out:
  if (matrix_modified && !dec->inverted) {
    for (uint32_t gap = 0; gap < dec->K; gap++) {
      if (!compat_bitmask_check(&dec->repair_mask, gap))
        memset(nanorq_core_get_symbol_ptr(&core, dec->D, dec->stride, gap), 0, dec->stride);
    }
    for (uint32_t extra = 0; extra < overhead; extra++)
      memset(nanorq_core_get_symbol_ptr(&core, dec->D, dec->stride, core.P.Kprime + extra),
             0, dec->stride);
  }
  return success;
}

// similar to nanorq_generate_symbols(), but this is used at decoding.
int nanorq_regenerate_symbols(nanorq *rq) {
  if (!rq)
    return NANORQ_SYM_PARA;
  struct block_encoder *dec = rq->encoder;

  size_t conventional_count = dec->repair_bin.n;

  size_t num_gaps = compat_bitmask_gaps(&dec->repair_mask, dec->K);
  if (num_gaps == 0)
    return NANORQ_SYM_IGN;	// Because all source symbols exist, repair isn't needed.
  if (dec->inverted)
    return NANORQ_SYM_OK;	// intermediate symbols were generated already.

  if (conventional_count < num_gaps || conventional_count - num_gaps > UINT32_MAX)
    return NANORQ_SYM_LACK;	// added symbols aren't enough for recovery.

  uint32_t overhead = (uint32_t)(conventional_count - num_gaps);
  nanorq_core core = {0};
  bool matrix_modified = false;
  int success = NANORQ_SYM_ERR;

  if (!ensure_block_matrix(dec) ||
      !nanorq_core_encoder_new(dec->K, overhead, &core))
    goto out;

  size_t rows = nanorq_core_get_pc_rows(&core);
  size_t base_rows = nanorq_core_get_pc_rows(&dec->core);
  if (base_rows > rows)
    goto out;
  if (rows > dec->matrix_rows) {
    uint8_t *D = (uint8_t *)obl_alloc(rows, dec->stride, nanorq_oblas.align_size);
    if (!D)
      goto out;
    memcpy(D, dec->D, base_rows * dec->stride);
    obl_free(dec->D);
    dec->D = D;
    dec->matrix_rows = rows;
  } else {
    /* reset clears base matrix; clear only active overhead */
    memset(dec->D + base_rows * dec->stride, 0,
           (rows - base_rows) * dec->stride);
  }

  size_t prep_len = nanorq_core_calculate_prepare_memory(&core);
  if (!reserve_scratch(&dec->prep_mem, &dec->prep_len, prep_len) ||
      !nanorq_core_generate(&core, dec->prep_mem, dec->prep_len))
    goto out;

  size_t rep_idx = 0;
  matrix_modified = true;
  for (uint32_t gap = 0; gap < dec->K; gap++) {
    if (compat_bitmask_check(&dec->repair_mask, gap))
      continue;
    if (rep_idx == dec->repair_bin.n)
      goto out;
    repair_sym *rs = &dec->repair_bin.a[rep_idx++];
    nanorq_core_replace_symbol(&core, gap, rs->esi);
    nanorq_core_place_symbol(&core, dec->D, dec->stride, gap, rs->row, rq->T);
  }

  for (uint32_t extra = 0; extra < overhead; extra++) {
    if (rep_idx == dec->repair_bin.n)
      goto out;
    repair_sym *rs = &dec->repair_bin.a[rep_idx++];
    uint32_t row = core.P.Kprime + extra;
    nanorq_core_replace_symbol(&core, row, rs->esi);
    nanorq_core_place_symbol(&core, dec->D, dec->stride, row, rs->row, rq->T);
  }

  if (!nanorq_core_patch_matrix(&core))
    goto out;
  size_t work_len = nanorq_core_calculate_prepared_work_memory(&core);
  if (!reserve_scratch(&dec->work_mem, &dec->work_len, work_len) ||
      !reset_schedule(&dec->S, dec->K))
    goto out;
  nanorq_core_set_op_callback(&core, &dec->S, ops_push);
  int invert = nanorq_core_precalculate(&core, dec->work_mem, dec->work_len);
  if ((invert <= 0) || dec->S.overflowed || dec->S.cpidx != 2){
    if (invert == 0)
      success = NANORQ_SYM_LACK;	// Failed to invert matrix. It may need more symbols.
    goto out;
  }

  dec->core = core;
  ops_run(&dec->core, dec->D, dec->stride, &dec->S);
  dec->inverted = true;
  success = NANORQ_SYM_OK;

out:
  if (matrix_modified && !dec->inverted) {
    for (uint32_t gap = 0; gap < dec->K; gap++) {
      if (!compat_bitmask_check(&dec->repair_mask, gap))
        memset(nanorq_core_get_symbol_ptr(&core, dec->D, dec->stride, gap), 0, dec->stride);
    }
    for (uint32_t extra = 0; extra < overhead; extra++)
      memset(nanorq_core_get_symbol_ptr(&core, dec->D, dec->stride, core.P.Kprime + extra),
             0, dec->stride);
  }
  return success;
}

// similar to nanorq_encode(), but this is used at decoding.
int nanorq_restore_symbol(nanorq *rq, void *data, uint32_t esi) {
  if (!rq || !data)
    return NANORQ_SYM_PARA;
  if (esi >= (1 << 24) || esi > rq->max_esi)
    return NANORQ_SYM_PARA;
  struct block_encoder *dec = rq->encoder;

  if (!dec->inverted){	// check intermediate symbols
    int success = nanorq_regenerate_symbols(rq);
    if (success == NANORQ_SYM_ERR || success == NANORQ_SYM_LACK)
      return success;
  }

  if (!ensure_mix_buf(dec))
    return NANORQ_SYM_ERR;
  ops_mix(&dec->core, dec->D, dec->stride, esi, dec->mix_buf);
  memcpy(data, dec->mix_buf, rq->T);

  return NANORQ_SYM_OK;
}
