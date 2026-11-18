#ifndef NANORQC_H
#define NANORQC_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "io.h"

#ifdef __cplusplus
extern "C" {
#endif

#define NANORQ_SYM_LACK 3	// It needs to add more symbols for rpeair.
#define NANORQ_SYM_DUP  2	// This ESI symbol was added already.
#define NANORQ_SYM_IGN  1	// All source symbols were added.
#define NANORQ_SYM_OK   0
#define NANORQ_SYM_ERR  -1	// Fatal error (mostly insufficient memory)
#define NANORQ_SYM_PARA -2	// Given parameter is bad.

typedef struct nanorq nanorq;

/*
 It returns a new encoder configured with given parameters.
 K: number of source symbols
 T: size of each symbol in bytes
K must be 1 ~ 56403. T doesn't need to be multiple of 4.
*/
nanorq *nanorq_encoder_new(uint16_t K, uint32_t T);

/*
 It adds a series of source symbols to the encoder.
 data: input bytes of adding source symbols (size must be esi_count * T)
 esi_start: the first Encoding Symbol ID (range is 0 ~ K-1)
 esi_count: number of adding source symbols (range is 1 ~ K)
It returns 0 at success, see NANORQ_SYM_*.
This is an optional function.
If this isn't called, source symbols must be added later in other functions.
*/
int nanorq_encoder_add_symbols(nanorq *rq, void *data, int esi_start, int esi_count);

/*
 It generates intermediate symbols from source symbols.
 src: input bytes of all source symbols (total size must be K * T)
      When all source symbols were added already, src can be NULL.
It returns 0 at success, see NANORQ_SYM_*.
This is an optional function.
If this isn't called, it will be generated later.
*/
int nanorq_generate_symbols(nanorq *rq, void *src);

/*
 It frees up any resources used by a decoder/encoder.
After nanorq's task was done, call this function.
*/
void nanorq_free(nanorq *rq);

/*
 It precalculates precode matrix inversion.
It returns 0 at success, see NANORQ_SYM_*
This is an optional function.
If this isn't called, it will be calculated later automatically.
*/
int nanorq_precalculate(nanorq *rq);

/*
 It encodes a symbol for a given esi.
 data: output symbol bytes for the ESI (size must be T)
 esi: Encoding Symbol ID, 0 ~ K-1 for source symbols, K ~ 1<<24-1 for repair symbols
      When esi < K, source symbol is restored.
      When esi >= K, repair symbol is created.
 src: input bytes of all source symbols (total size must be K * T)
      If all source symbols were added or intermediate symbols were generated, src can be NULL.
It returns 0 at success, see NANORQ_SYM_*.
Though this can retrieve a source symbol also, just copying from src would be easier.
So, it will create a repair symbol mostly.
*/
int nanorq_encode(nanorq *rq, void *data, uint32_t esi, void *src);

// cleanup encoder resources
void nanorq_encoder_cleanup(nanorq *rq);

/* reset encoder state for reuse; keep allocations until cleanup */
void nanorq_encoder_reset(nanorq *rq);

/*
 It returns a new decoder initialized with given parameters.
 K: number of source symbols, which MUST be same as K in encoder
 T: size of each symbol in bytes, which MUST be same as T in encoder
*/
nanorq *nanorq_decoder_new(uint16_t K, uint32_t T);

/*
 It sets the largest encoding symbol id allowed.
 max_esi: max ESI of possible repair symbols (must be larger than K)
          Set the largest Encoding Symbol ID allowed at decoding.
          If number of possible repair symbols is R, their ESIs are K ~ K+R-1.
          When total number of source symbols and possible repair symbols is K+R = N,
          max_esi is K+R-1 = N-1.
It returns 0 at success, see NANORQ_SYM_*.
When this isn't called, max_esi = 1<<24-1 by default.
This is an optional utility function. No need to use this mostly.
*/
int nanorq_set_max_esi(nanorq *rq, uint32_t max_esi);

/*
 It adds a symbol to the decoder.
 data: input symbol bytes for the ESI (size must be T)
 esi: Encoding Symbol ID, 0 ~ K-1 for source symbols, K ~ 1<<24-1 for repair symbols
 src: output bytes of all source symbols (total size must be K * T)
      If output isn't needed, set src to NULL. Or else, speed may become slower.
      When adding a source symbol and src isn't NULL, the symbol is copied to src.
It returns 0 at success, see NANORQ_SYM_*.
*/
int nanorq_decoder_add_symbol(nanorq *rq, void *data, uint32_t esi, void *src);

/*
 It returns number of added symbols in encoder/decoder.
This is an optional utility function.
The maximum returned value is K in encoder by nanorq_encoder_add_symbols(),
or it's 1<<24 in decoder by nanorq_decoder_add_symbol().
*/
int nanorq_num_loaded(nanorq *rq);

// returns number of symbol gaps in decoder
int nanorq_num_missing(nanorq *rq);

// returns number of repair symbols in decoder
int nanorq_num_repair(nanorq *rq);

/*
 It repairs all missing source symbols in decoder at once.
 src: output bytes of all source symbols (total size must be K * T)
      Missing source symbols will be restored in their positions.
      While this is useful, src buffer size may be large.
      If missing symbols are fewer than K, nanorq_restore() would consume less memory.
It returns 0 at success, see NANORQ_SYM_*.
After it repaired all missing symbols, nanorq_num_missing() will return 0.
Because all source symbols are repaired, calling this twice is useless.
*/
int nanorq_repair(nanorq *rq, void *src);

/*
 It regenerates intermediate symbols from added symbols in decoder.
It returns 0 at success, see NANORQ_SYM_*.
This is an optional function.
It determines possible repair before calling nanorq_repair() or nanorq_restore_symbol().
If you want to know number of added symbols are enough or not before actual repair,
this function will be useful.
*/
int nanorq_regenerate_symbols(nanorq *rq);

/*
 It restores a symbol for a given esi in decoder.
 data: output symbol bytes for the ESI (size must be T)
 esi: Encoding Symbol ID, 0 ~ K-1 for source symbols, K ~ 1<<24-1 for repair symbols
      When esi < K, source symbol is recovered.
      When esi >= K, repair symbol is created.
It returns 0 at success, see NANORQ_SYM_*.
This function is similar to nanorq_encode() for encoder.
If missing symbols are few, this may be useful than nanorq_repair().
*/
int nanorq_restore_symbol(nanorq *rq, void *data, uint32_t esi);

#ifdef __cplusplus
}
#endif

#endif
