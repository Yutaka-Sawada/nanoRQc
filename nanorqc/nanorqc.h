#ifndef NANORQ_H
#define NANORQ_H

#define NANORQ_SYM_END 3	// It got enough symbols already. Don't add anymore.
#define NANORQ_SYM_DUP 2	// This ESI symbol was added already.
#define NANORQ_SYM_IGN 1	// Repair isn't needed.
#define NANORQ_SYM_ADDED 0
#define NANORQ_SYM_ERR -1

#define NANORQ_ENC_ERR -1	// Fatal error (mostly insufficient memory)
#define NANORQ_ENC_INV -2	// It failed to invert matrix (may require more symbols).
#define NANORQ_ENC_LACK -3	// It needs to add more symbols.

typedef struct nanorq nanorq;

/*
 It returns a new encoder configured with given parameters.
 K: number of source symbols
 T: size of each symbol in bytes
K must be 1 ~ 56403. T isn't needed to be multiple of 4.
*/
nanorq *nanorq_encoder_new(int K, unsigned int T);

/*
 It adds a source symbol to the encoder.
 data: input symbol bytes for the ESI (size must be T)
 esi: Encoding Symbol ID, 0 ~ K-1
It returns 0 at success, see NANORQ_SYM_*.
This is an optional function.
If it's not called, they are added in nanorq_generate_symbols() or nanorq_encode_all().
*/
int nanorq_encoder_add_symbol(nanorq *rq, void *data, unsigned int esi);

/*
 It returns number of added symbols in encoder/decoder.
This is an optional utility function.
The maximum returned value is K in encoder by nanorq_encoder_add_symbol(),
or the max is 1<<24 in decoder by nanorq_decoder_add_symbol().
After intermediate symbols were calculated, it returns 1<<25.
*/
int nanorq_num_added(nanorq *rq);

/*
 It generates intermediate symbols from source symbols.
 src: input bytes of all source symbols (total size must be K * T)
      When all source symbols were added already, src can be NULL.
It returns 0 at success, see NANORQ_ENC_*.
This is an optional function.
If it's not called, it's calculated in nanorq_encode_all() or nanorq_encode().
*/
int nanorq_generate_symbols(nanorq *rq, void *src);

/*
 It frees up any resources used by a decoder/encoder.
After nanorq's task was done, call this function.
*/
void nanorq_free(nanorq *rq);

/*
 It precalculates precode matrix inversion.
It returns 0 at success, see NANORQ_ENC_*
This is an optional function.
If it's not called, it does in nanorq_generate_symbols() or nanorq_encode_all() or nanorq_encode().
*/
int nanorq_precalculate(nanorq *rq);

/*
 It encodes all repair symbols at once.
 data: output bytes of all repair symbols (total size must be R * T)
 R: number of creating repair symbols (their ESIs are K ~ K+R-1).
 src: input bytes of all source symbols (total size must be K * T)
      If all source symbols were added or intermediate symbols were calculated, src can be NULL.
It returns 0 at success, see NANORQ_ENC_*.
While this is useful, data buffer size may be large.
When creating symbols are many, nanorq_encode() would consume less memory.
 */
int nanorq_encode_all(nanorq *rq, void *data, int R, void *src);

/*
 It encodes one symbol for a given esi.
 data: output symbol bytes for the ESI (size must be T)
 esi: Encoding Symbol ID, 0 ~ K-1 for source symbols, K ~ 1<<24-1 for repair symbols
      When esi < K, source symbol is restored.
      When esi >= K, repair symbol is created.
 src: input bytes of all source symbols (total size must be K * T)
      If all source symbols were added or intermediate symbols were calculated, src can be NULL.
It returns 0 at success, see NANORQ_ENC_*.
Though this can retrieve a source symbol also, just copying from src would be easier.
So, it will create a repair symbol mostly.
*/
int nanorq_encode(nanorq *rq, void *data, unsigned int esi, void *src);

/*
 It will cleanup encoder resources.
This is an optional utility function.
nanorq_free() calls this internally.
*/
void nanorq_encoder_cleanup(nanorq *rq);

/*
 It will reset internal state.
This is an optional utility function.
No need to use this mostly.
*/
void nanorq_encoder_reset(nanorq *rq);

/*
 It returns a new decoder initialized with given parameters.
 K: number of source symbols, which MUST be same as K in encoder
 T: size of each symbol in bytes, which MUST be same as T in encoder
 max_esi: max ESI of possible repair symbols (should be larger than K)
          Set the largest Encoding Symbol ID allowed at decoding.
          If number of possible repair symbols is R, their ESIs are K ~ K+R-1.
          When total number of source symbols and possible repair symbols is K+R = N,
          max_esi is K+R-1 = N-1.
*/
nanorq *nanorq_decoder_new(int K, unsigned int T, unsigned int max_esi);

/*
 It adds a symbol to the decoder.
 data: input symbol bytes for the ESI (size must be T)
 esi: Encoding Symbol ID, 0 ~ K-1 for source symbols, K ~ 1<<24-1 for repair symbols
 src: output bytes of all source symbols (total size must be K * T)
      If output isn't needed, set src to NULL. Or else, speed may become slower.
      When adding a source symbol and src isn't NULL, the symbol is copied to src.
It returns 0 at success, see NANORQ_SYM_*.
*/
int nanorq_decoder_add_symbol(nanorq *rq, void *data, unsigned int esi, void *src);

/*
 It returns number of missing source symbols in decoder/encoder.
This is an optional utility function.
Contrary to nanorq_num_added(), nanorq_num_missing() returns number of absent symbols.
After nanorq_repair_all() done successfully, nanorq_num_missing() returns 0.
*/
int nanorq_num_missing(nanorq *rq);

/*
 It returns number of added repair symbols in decoder.
This is an optional utility function.
While nanorq_num_added() returns total number of both source symbols and repair symbols,
nanorq_num_repair() points number of repair symbols only.
*/
int nanorq_num_repair(nanorq *rq);

/*
 It repairs all missing source symbols in decoder at once.
 src: output bytes of all source symbols (total size must be K * T)
      Missing source symbols will be restored in their positions.
      While this is useful, src buffer size may be large.
      If missing symbols are fewer than K, nanorq_restore() would consume less memory.
It returns 0 at success, see NANORQ_ENC_*.
After it repaired all missing symbols, nanorq_num_missing() will return 0.
Because all source symbols are repaired, calling this twice will be useless.
*/
int nanorq_repair_all(nanorq *rq, void *src);

/*
 it restores one symbol for a given esi in decoder.
 data: output symbol bytes for the ESI (size must be T)
 esi: Encoding Symbol ID, 0 ~ K-1 for source symbols, K ~ 1<<24-1 for repair symbols
      When esi < K, source symbol is recovered.
      When esi >= K, repair symbol is created like nanorq_encode().
It returns 0 at success, see NANORQ_ENC_*.
*/
int nanorq_restore(nanorq *rq, void *data, unsigned int esi);

/*
 It recovers intermediate symbols from added symbols in decoder.
It returns 0 at success, see NANORQ_ENC_*.
This is an optional function.
It determines possible repair before calling nanorq_repair_all() or nanorq_restore().
If it's not called, it's calculated in nanorq_repair_all() or nanorq_restore().
If you want to know number of added repair symbols are enough or not before actual repair,
this function would be useful.
*/
int nanorq_recover_symbols(nanorq *rq);

#endif
