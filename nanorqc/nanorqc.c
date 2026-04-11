#include <stdio.h>

#include "bitmask.h"
#include "nanorqc.h"
#include "precode.h"
#include "tuple.h"

#define LOAD_OK (1 << 25)	// This value is larger than max ESI.

#define APPLYROW(D, w, a, k)                                                   \
  do {                                                                         \
    uint8_t *tmp = om_R(*D, w);                                                \
    for (unsigned i = 0; i < k; i++)                                           \
      a[i] ^= tmp[i];                                                          \
  } while (0);

#define APPLYROW8(D, w, a, k)                                                  \
  do {                                                                         \
    uint8_t *tmp = om_R(*D, w);                                                \
    unsigned i, k8 = k / 8;                                                    \
    for (i = 0; i < k8; i++)                                                   \
      ((uint64_t *)a)[i] ^= ((uint64_t *)tmp)[i];                              \
    for (i = i * 8; i < k; i++)                                                \
      a[i] ^= tmp[i];                                                          \
  } while (0);

struct nanorq {
	uint32_t T;
	uint32_t max_esi;
	params P;
	schedule *S;
	// this supports single block only
	uint16_t K;
	int loaded;
	octmat D;
	repair_vec repair_bin;
	bitmask esi_mask;
};

static bool set_block_encoder(nanorq *rq) {
	rq->esi_mask = bitmask_new(rq->max_esi + 1);
	if (rq->esi_mask.a == NULL)
		return false;

	int spare = 0;
	if (rq->max_esi > rq->K)	// need spare for decoding
		spare = rq->max_esi - rq->K;
	/*
	When there is only one repair symbol, max_esi is K. Then, spare is 0.
	When there are two repair symbols, max_esi is K+1. Then, spare is 1. (one and another)
	When there are three repair symbols, max_esi is K+2. Then, spare is 2. (one and other)
	*/

	om_resize(&rq->D, rq->P.L + spare, rq->T);	// number of intermediate symbols is L
	if (om_P(rq->D) == NULL){
		if (rq->esi_mask.a)
			bitmask_free(&rq->esi_mask);
		return false;
	}

	return true;
}

static void load_source_symbols(nanorq *rq, uint8_t *src) {
	for (int esi = 0; esi < rq->K; esi++){
		if (!bitmask_check(&rq->esi_mask, esi)){
			memcpy(om_R(rq->D, rq->P.S + rq->P.H + esi), src + (size_t)esi * rq->T, rq->T);
			bitmask_set(&rq->esi_mask, esi);
			rq->loaded += 1;
		}
	}
}

static void decode_row(params *P, octmat *D, uint32_t row, uint8_t *ptr, uint32_t len) {
	tuple t = gen_tuple(row, P);
	// Because "0 ^ D = D", use memcpy() instead of zero fill and XOR.
	//memset(ptr, 0, len);
	//APPLYROW(D, t.b, ptr, len);
	memcpy(ptr, om_R(*D, t.b), len);

	if (((size_t)ptr) & 7){	// When ptr isn't aligned to 8, XOR per 1-byte.
		for (unsigned j = 1; j < t.d; j++) {
			t.b = (t.b + t.a) % P->W;
			APPLYROW(D, t.b, ptr, len);
		}
		while (t.b1 >= P->P)
			t.b1 = (t.b1 + t.a1) % P->P1;

		APPLYROW(D, P->W + t.b1, ptr, len);
		for (unsigned j = 1; j < t.d1; j++) {
			t.b1 = (t.b1 + t.a1) % P->P1;
			while (t.b1 >= P->P)
				t.b1 = (t.b1 + t.a1) % P->P1;
			APPLYROW(D, P->W + t.b1, ptr, len);
		}

	} else {	// When ptr is multiple of 8, XOR per 8-bytes at first.
		for (unsigned j = 1; j < t.d; j++) {
			t.b = (t.b + t.a) % P->W;
			APPLYROW8(D, t.b, ptr, len);
		}
		while (t.b1 >= P->P)
			t.b1 = (t.b1 + t.a1) % P->P1;

		APPLYROW8(D, P->W + t.b1, ptr, len);
		for (unsigned j = 1; j < t.d1; j++) {
			t.b1 = (t.b1 + t.a1) % P->P1;
			while (t.b1 >= P->P)
				t.b1 = (t.b1 + t.a1) % P->P1;
			APPLYROW8(D, P->W + t.b1, ptr, len);
		}
	}
}

/*
 * data: input symbol bytes for the ESI (size must be T)
 * esi: Encoding Symbol ID, 0 ~ K-1
 */
int nanorq_encoder_add_symbol(nanorq *rq, void *data, unsigned int esi) {
	if (esi >= rq->K)
		return NANORQ_SYM_ERR;

	if (rq->loaded == LOAD_OK)
		return NANORQ_SYM_END;	// got enough symbols already

	if (bitmask_check(&rq->esi_mask, esi))
		return NANORQ_SYM_DUP; // already got this esi

	memcpy(om_R(rq->D, rq->P.S + rq->P.H + esi), data, rq->T);
	bitmask_set(&rq->esi_mask, esi);
	rq->loaded += 1;

	return NANORQ_SYM_ADDED;
}

/*
 * src: input bytes of all source symbols (total size must be K * T)
 */
int nanorq_generate_symbols(nanorq *rq, void *src) {
	if (rq->loaded == LOAD_OK)
		return 0;
	if (rq->loaded < rq->K){
		if (src == NULL)
			return NANORQ_ENC_LACK;
		load_source_symbols(rq, src);	// load all source symbols at once
	}

	schedule *S = NULL;
	if (rq->S) {
		S = rq->S;
	} else {
		spmat *A = precode_matrix_gen(&rq->P, 0);
		if (A == NULL)
			return NANORQ_ENC_ERR;
		S = precode_matrix_invert(&rq->P, A);
	}
	if (S == NULL)
		return NANORQ_ENC_INV;
	precode_matrix_intermediate(&rq->P, &rq->D, S);
	if (rq->S == NULL)
		sched_free(S);
	rq->loaded = LOAD_OK;	// sign of complete intermediate symbols
	return 0;
}

/*
 * K: number of source symbols
 * T: size of each symbol in bytes
 */
nanorq *nanorq_encoder_new(int K, unsigned int T) {
	if ((K < 1) || (K > K_max) || (T < 1))
		return NULL;
	nanorq *rq = calloc(1, sizeof(nanorq));
	if (rq == NULL)
		return NULL;

	rq->T = T;
	rq->K = K;

	rq->P = params_init((uint16_t)rq->K);

	// ESI for source symbols are 0 ~ K-1.
	rq->max_esi = K - 1;

	if (!set_block_encoder(rq)){
		free(rq);
		return NULL;
	}

	return rq;
}

void nanorq_free(nanorq *rq) {
	if (rq) {
		if (rq->S)
			sched_free(rq->S);
		nanorq_encoder_cleanup(rq);
		free(rq);
		rq = NULL;
	}
}

/*
 * K: number of source symbols, which MUST be same as K in encoder
 * T: size of each symbol in bytes, which MUST be same as T in encoder
 * max_esi: max ESI of possible repair symbols (should be larger than K)
 */
nanorq *nanorq_decoder_new(int K, unsigned int T, unsigned int max_esi) {
	if ((K < 1) || (K > K_max) || (T < 1))
		return NULL;
	nanorq *rq = calloc(1, sizeof(nanorq));
	if (rq == NULL)
		return NULL;

	rq->T = T;
	rq->K = K;

	rq->P = params_init((uint16_t)rq->K);

	// If max_esi is strange, set default value automatically.
	if ((max_esi >= (1 << 24) || max_esi < rq->K)){
		rq->max_esi = rq->K * 2;	// redundancy 100% and 1 overhead
	} else {
		rq->max_esi = max_esi;
	}

	if (!set_block_encoder(rq)){
		free(rq);
		return NULL;
	}

	return rq;
}

int nanorq_precalculate(nanorq *rq) {
	params *P = &rq->P;
	spmat *A = precode_matrix_gen(P, 0);
	if (A == NULL)
		return NANORQ_ENC_ERR;
	schedule *S = precode_matrix_invert(P, A);
	if (S == NULL)
		return NANORQ_ENC_INV;
	rq->S = S;
	return 0;
}

/*
 * data: output symbol bytes for the ESI (size must be T)
 * esi: Encoding Symbol ID, 0 ~ K-1 for source symbols, K ~ 1<<24-1 for repair symbols
 *      When esi < K, source symbol is restored like nanorq_restore().
 *      When esi >= K, repair symbol is created.
 * src: input bytes of all source symbols (total size must be K * T)
 *      If all source symbols were added or intermediate symbols were calculated, src can be NULL.
 */
int nanorq_encode(nanorq *rq, void *data, unsigned int esi, void *src) {
	if (esi < rq->K) {	// source symbol
		if (rq->loaded == LOAD_OK) {	// calculate a source symbol from intermediate symbols
			decode_row(&rq->P, &rq->D, esi, data, rq->T);
		} else {	// intermediate symbols aren't generated yet
			if (!bitmask_check(&rq->esi_mask, esi))
				load_source_symbols(rq, src);	// If the source symbol wasn't loaded yet, load all now.
			memcpy(data, om_R(rq->D, rq->P.S + rq->P.H + esi), rq->T);
		}
	} else {	// repair symbol
		if (esi > ((1 << 24) - 1))
			return NANORQ_ENC_ERR;
		// esi is for repair symbol
		if (rq->loaded != LOAD_OK){	// intermediate symbols aren't generated yet
			int err = nanorq_generate_symbols(rq, src);
			if (err != 0)
				return err;
		}
		// calculate a repair symbol from intermediate symbols
		uint32_t isi = esi + (rq->P.Kprime - rq->K);
		decode_row(&rq->P, &rq->D, isi, data, rq->T);
	}
	return 0;
}

/*
 * data: output bytes of all repair symbols (total size must be R * T)
 * R: number of creating repair symbols (their ESIs are K ~ K+R-1).
 * src: input bytes of all source symbols (total size must be K * T)
 *      If all source symbols were added or intermediate symbols were calculated, src can be NULL.
 */
int nanorq_encode_all(nanorq *rq, void *data, int R, void *src) {
	if ((R < 1) || (R > (1 << 24) - rq->K))
		return NANORQ_ENC_ERR;

	if (rq->loaded != LOAD_OK){	// intermediate symbols aren't calculated yet
		int err = nanorq_generate_symbols(rq, src);
		if (err != 0)
			return err;
	}

	// calculate every repair symbol from intermediate symbols
	for (int i = 0; i < R; i++){
		uint32_t isi = i + rq->P.Kprime;
		decode_row(&rq->P, &rq->D, isi, (uint8_t *)data + (size_t)i * rq->T, rq->T);
	}
	return 0;
}

void nanorq_encoder_cleanup(nanorq *rq) {
	om_destroy(&rq->D);
	if (kv_size(rq->repair_bin) > 0) {
		for (int rs = 0; rs < kv_size(rq->repair_bin); rs++)
			om_destroy(&(kv_A(rq->repair_bin, rs).row));
		kv_destroy(rq->repair_bin);
	}
	if (kv_size(rq->esi_mask) > 0)
		bitmask_free(&rq->esi_mask);
}

void nanorq_encoder_reset(nanorq *rq) {
	rq->loaded = 0;
	if (om_P(rq->D))
		memset(om_P(rq->D), 0, (size_t)rq->D.rows * rq->D.cols_al);
	if (kv_size(rq->repair_bin) > 0) {
		for (int rs = 0; rs < kv_size(rq->repair_bin); rs++)
			om_destroy(&(kv_A(rq->repair_bin, rs).row));
		kv_destroy(rq->repair_bin);
		kv_init(rq->repair_bin);
	}
	if (kv_size(rq->esi_mask) > 0)
		bitmask_reset(&rq->esi_mask);
}

/*
 * data: input symbol bytes for the ESI (size must be T)
 * esi: Encoding Symbol ID, 0 ~ K-1 for source symbols, K ~ 1<<24-1 for repair symbols
 * src: output bytes of all source symbols (total size must be K * T)
 *      If output isn't needed, set src to NULL. Or else, speed may become slower.
 */
int nanorq_decoder_add_symbol(nanorq *rq, void *data, unsigned int esi, void *src) {
	if (esi > rq->max_esi)
		return NANORQ_SYM_ERR;

	if (rq->loaded == LOAD_OK)
		return NANORQ_SYM_END;	// got enough symbols already

	if (bitmask_gaps(&rq->esi_mask, rq->K) == 0)
		return NANORQ_SYM_IGN; // no repair needed.

	if (bitmask_check(&rq->esi_mask, esi))
		return NANORQ_SYM_DUP; // already got this esi

	if (esi < rq->K) {
		// write original symbol to decode mat and output to src
		memcpy(om_R(rq->D, rq->P.S + rq->P.H + esi), data, rq->T);
		if (src != NULL)
			memcpy((uint8_t *)src + (size_t)esi * rq->T, data, rq->T);
	} else {
		// save repair symbol for precode patching
		repair_sym rs = {esi, OM_INITIAL};
		om_resize(&rs.row, 1, rq->T);
		if (om_P(rs.row) == NULL)
			return NANORQ_SYM_ERR;
		memcpy(om_R(rs.row, 0), data, rq->T);
		kv_push(repair_sym, rq->repair_bin, rs);
	}
	bitmask_set(&rq->esi_mask, esi);
	rq->loaded += 1;

	return NANORQ_SYM_ADDED;
}

int nanorq_num_added(nanorq *rq) {
	return rq->loaded;
}

int nanorq_num_missing(nanorq *rq) {
	return bitmask_gaps(&rq->esi_mask, rq->K);
}

int nanorq_num_repair(nanorq *rq) {
	return kv_size(rq->repair_bin);
}

static void patch_precode_matrix(params *P, spmat *A, uint16_t K, int num_gaps,
                                 bitmask *mask, repair_vec *repair_bin) {
  unsigned padding = P->Kprime - K;
  int rep_idx = 0;
  for (int gap = 0; gap < P->L && num_gaps > 0; gap++) {
    if (bitmask_check(mask, gap))
      continue;
    int row = gap + P->H + P->S;
    uint32_t esi = kv_A(*repair_bin, rep_idx++).esi + padding;

    spmat_clear_row(A, row);
    params_set_idxs(esi, P, &A->idxs[row]);
    num_gaps--;
  }

  for (int row = P->L; row < (int)(A->rows); row++) {
    uint32_t esi = kv_A(*repair_bin, rep_idx++).esi + padding;
    spmat_clear_row(A, row);
    params_set_idxs(esi, P, &A->idxs[row]);
  }
}

static void fill_symbol_matrix_gaps(params *P, octmat *D, uint16_t K,
                                    bitmask *esi_mask,
                                    repair_vec *repair_bin) {
  int rep_idx = 0, skip = P->S + P->H, num_repair = kv_size(*repair_bin);
  for (int gap = 0; gap < K && rep_idx < num_repair; gap++) {
    if (bitmask_check(esi_mask, gap))
      continue;
    int row = skip + gap;
    repair_sym rs = kv_A(*repair_bin, rep_idx++);
    memcpy(om_R(*D, row), om_P(rs.row), D->cols);
  }

  for (int row = P->L; rep_idx < num_repair; row++) {
    repair_sym rs = kv_A(*repair_bin, rep_idx++);
    memcpy(om_R(*D, row), om_P(rs.row), D->cols);
  }
}

static void decoder_reduce(nanorq *rq) {
	// Reduce memory of spare after intermediate symbols
	if ((om_P(rq->D) != NULL) && (rq->D.rows > rq->P.L))
		om_resize(&rq->D, rq->P.L, rq->D.cols);

	// For nanorq_num_repair(), it doesn't reset repair_bin.
	if (kv_size(rq->repair_bin) > 0) {
		for (int rs = 0; rs < kv_size(rq->repair_bin); rs++)
			om_destroy(&(kv_A(rq->repair_bin, rs).row));
	}
}

/*
 * src: output bytes of all source symbols (total size must be K * T)
 *      Missing source symbols will be restored in their positions.
 */
int nanorq_repair_all(nanorq *rq, void *src) {
	params *P = &rq->P;
	octmat *D = &rq->D;
	bitmask *esi_mask = &rq->esi_mask;

	if (rq->loaded < rq->K)
		return NANORQ_ENC_LACK;	// symbols aren't enough
	int num_gaps = bitmask_gaps(esi_mask, rq->K);
	if (num_gaps == 0)
		return 0;	// no repair needed.

	if (rq->loaded != LOAD_OK){
		int overhead = rq->loaded - rq->K;
		if ((int)(D->rows) < P->L + overhead)
			return NANORQ_ENC_LACK;

		repair_vec *repair_bin = &rq->repair_bin;
		fill_symbol_matrix_gaps(P, D, rq->K, esi_mask, repair_bin);
		spmat *A = precode_matrix_gen(P, overhead);
		if (A == NULL)
			return NANORQ_ENC_ERR;
		patch_precode_matrix(P, A, rq->K, num_gaps, esi_mask, repair_bin);

		schedule *S = precode_matrix_invert(P, A);
		if (S == NULL)
			return NANORQ_ENC_INV;
		if (!precode_matrix_intermediate(P, D, S))
			return NANORQ_ENC_ERR;
		sched_free(S);
		decoder_reduce(rq);
		rq->loaded = LOAD_OK;
	}

	// decode missing rows
	for (int esi = 0; esi < rq->K; esi++) {
		if (bitmask_check(esi_mask, esi))
			continue;	// skip existing source symbols
		decode_row(P, D, esi, (uint8_t *)src + (size_t)esi * rq->T, rq->T);
		bitmask_set(esi_mask, esi);	// the esi was repaired
	}
	return 0;	// all missing source symbols were repaired.
}

/*
 * similar to nanorq_encode(), but this is used at decoding.
 * data: output symbol bytes for the ESI (size must be T)
 * esi: Encoding Symbol ID, 0 ~ K-1 for source symbols, K ~ 1<<24-1 for repair symbols
 *      When esi < K, source symbol is recovered.
 *      When esi >= K, repair symbol is created like nanorq_encode().
 */
int nanorq_restore(nanorq *rq, void *data, unsigned int esi) {
	params *P = &rq->P;
	octmat *D = &rq->D;
	bitmask *esi_mask = &rq->esi_mask;

	if ((esi < rq->K) && (rq->loaded != LOAD_OK) && (bitmask_check(&rq->esi_mask, esi) != 0)){
		// copy source symbol from decode mat
		memcpy(data, om_R(rq->D, rq->P.S + rq->P.H + esi), rq->T);
		return 0;
	}
	if (rq->loaded < rq->K)
		return NANORQ_ENC_LACK;	// symbols aren't enough

	if (rq->loaded != LOAD_OK){
		int overhead = rq->loaded - rq->K;
		if ((int)(D->rows) < P->L + overhead)
			return NANORQ_ENC_LACK;

		repair_vec *repair_bin = &rq->repair_bin;
		fill_symbol_matrix_gaps(P, D, rq->K, esi_mask, repair_bin);
		spmat *A = precode_matrix_gen(P, overhead);
		if (A == NULL)
			return NANORQ_ENC_ERR;
		int num_gaps = bitmask_gaps(esi_mask, rq->K);
		patch_precode_matrix(P, A, rq->K, num_gaps, esi_mask, repair_bin);

		schedule *S = precode_matrix_invert(P, A);
		if (S == NULL)
			return NANORQ_ENC_INV;
		if (!precode_matrix_intermediate(P, D, S))
			return NANORQ_ENC_ERR;
		sched_free(S);
		decoder_reduce(rq);
		rq->loaded = LOAD_OK;
	}

	if (esi < rq->K) {	// recover source symbol
		decode_row(&rq->P, &rq->D, esi, data, rq->T);
	} else {	// recreate repair symbol
		if (esi > ((1 << 24) - 1))
			return NANORQ_ENC_ERR;
		uint32_t isi = esi + (rq->P.Kprime - rq->K);
		decode_row(&rq->P, &rq->D, isi, data, rq->T);
	}

	return 0;
}

// similar to nanorq_generate_symbols(), but this is used at decoding.
int nanorq_recover_symbols(nanorq *rq) {
	if (rq->loaded == LOAD_OK)
		return 0;
	if (rq->loaded < rq->K)
		return NANORQ_ENC_LACK;	// symbols aren't enough

	params *P = &rq->P;
	octmat *D = &rq->D;
	bitmask *esi_mask = &rq->esi_mask;
	int num_gaps = bitmask_gaps(esi_mask, rq->K);
	if (num_gaps == 0)
		return 0;	// intermediate symbols aren't needed.

	int overhead = rq->loaded - rq->K;
	if ((int)(D->rows) < P->L + overhead)
		return NANORQ_ENC_LACK;

	repair_vec *repair_bin = &rq->repair_bin;
	fill_symbol_matrix_gaps(P, D, rq->K, esi_mask, repair_bin);
	spmat *A = precode_matrix_gen(P, overhead);
	if (A == NULL)
		return NANORQ_ENC_ERR;
	patch_precode_matrix(P, A, rq->K, num_gaps, esi_mask, repair_bin);

	schedule *S = precode_matrix_invert(P, A);
	if (S == NULL)
		return NANORQ_ENC_INV;
	if (!precode_matrix_intermediate(P, D, S))
		return NANORQ_ENC_ERR;
	sched_free(S);
	decoder_reduce(rq);
	rq->loaded = LOAD_OK;

	return 0;
}

