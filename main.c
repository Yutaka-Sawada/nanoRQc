
#include <inttypes.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "nanorqc/nanorqc.h"


/*
PCG Random Number Generation
http://www.pcg-random.org

I use global RNGs only.
*/
static uint64_t pcg32_state, pcg32_inc;

uint32_t pcg32_random()
{
    uint64_t oldstate = pcg32_state;
    pcg32_state = oldstate * 6364136223846793005ULL + pcg32_inc;
    uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((uint32_t)(-(int32_t)rot) & 31));
}

void pcg32_srandom(uint64_t seed, uint64_t seq)
{
	pcg32_state = 0U;
	pcg32_inc = (seq << 1u) | 1u;
	pcg32_random();
	pcg32_state += seed;
	pcg32_random();
}

uint32_t pcg32_boundedrand(uint32_t bound)
{
    uint32_t threshold = (uint32_t)(-(int32_t)bound) % bound;

    for (;;) {
        uint32_t r = pcg32_random();
        if (r >= threshold)
            return r % bound;
    }
}


/*
source_count = s#
parity_count = p#
lost_count = l#
each_size = d#
*/
int main(int argc, char* argv[]){
	char *p;
	char *source_buf, *parity_buf, *work_buf;
	int i, j, min, max, rv;
	int source_count = 1000, parity_count = 200, lost_count = 0;
	int input_count, esi;
	int *order_buf;
	unsigned int *int_p;
	unsigned int each_size = 64000;
	size_t source_size, parity_size;
	size_t malloc_size, input_size, output_size;
	double value_kb, value_mb, input_speed, output_speed;
	clock_t start, finish;
	double duration;
	nanorq *rq;

	// read input
	for (i = 1; i < argc; i++){
		p = argv[i];
		//printf("argv[%d] = %s\n", i, p);

		if ( (p[0] == 's') || (p[0] == 'S') ){
			source_count = atoi(p + 1);
		} else if ( (p[0] == 'p') || (p[0] == 'P') ){
			parity_count = atoi(p + 1);
		} else if ( (p[0] == 'l') || (p[0] == 'L') ){
			lost_count = atoi(p + 1);
		} else if ( (p[0] == 'e') || (p[0] == 'E') ){
			each_size = atoi(p + 1);

		} else {	// invalid option
			printf("options: s#, p#, l#, e#\ns# or S# = source count (1 ~ 56403)\np# or P# = parity count\nl# or L# = loss count\ne# or E# = each data size\n");
			return 0;
		}
	}

	// check input
	if (source_count < 1)
		source_count = 1;
	if (source_count > 56403)
		source_count = 56403;
	if (parity_count <= 0)
		parity_count = 1;
	if ( (lost_count <= 0) || (lost_count > parity_count) )
		lost_count = parity_count;
	if (lost_count > source_count)
		lost_count = source_count;

	printf("source count = %d, parity count = %d, loss count = %d\n", source_count, parity_count, lost_count);
	if (each_size < 1)
		each_size = 1;
	value_kb = (double)each_size / 1000;
	value_mb = (double)each_size / 1000000.f;
	printf("each data size = %d bytes = %g KB = %g MB\n", each_size, value_kb, value_mb);
	start = clock();
	printf("\n allocate memory and setup data ...\n");

	// allocate buffer
	source_size = (size_t)each_size * source_count;
	source_buf = malloc(source_size);
	if (source_buf == NULL){
		printf("Failed: malloc, %zd\n", source_size);
		return 1;
	}
	value_kb = (double)source_size / 1000;
	value_mb = (double)source_size / 1000000.f;
	printf("source data size = %zd bytes = %g KB = %g MB\n", source_size, value_kb, value_mb);

	// fill source data with random
	max = each_size / 4;
	for (i = 0; i < source_count; i++){
		int_p = (int *)(source_buf + (size_t)each_size * i);
		pcg32_srandom(0, i);
		for (j = 0; j < max; j++)
			int_p[j] = pcg32_random();
		if (each_size & 3){	// each_size may not be multiple of 4.
			rv = pcg32_random();
			if (each_size % 4 == 1){
				source_buf[(size_t)each_size * i + (each_size & ~3)] = (char)rv;
			} else if (each_size % 4 == 2){
				source_buf[(size_t)each_size * i + (each_size & ~3)    ] = (char)rv;
				source_buf[(size_t)each_size * i + (each_size & ~3) + 1] = (char)(rv >> 8);
			} else if (each_size % 4 == 3){
				source_buf[(size_t)each_size * i + (each_size & ~3)    ] = (char)rv;
				source_buf[(size_t)each_size * i + (each_size & ~3) + 1] = (char)(rv >> 8);
				source_buf[(size_t)each_size * i + (each_size & ~3) + 2] = (char)(rv >> 16);
			}
		}
	}

	// fill parity data with zero
	parity_size = (size_t)each_size * parity_count;
	parity_buf = calloc(parity_size, 1);
	if (parity_buf == NULL){
		printf("Failed: malloc, %zd\n", parity_size);
		free(source_buf);
		return 1;
	}
	value_kb = (double)parity_size / 1000;
	value_mb = (double)parity_size / 1000000.f;
	printf("parity data size = %zd bytes = %g KB = %g MB\n", parity_size, value_kb, value_mb);

	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf(" ... %.3f seconds\n", duration);

	printf("\n encode ...\n");
	start = clock();

	// initialize nanorq encoder
	rq = nanorq_encoder_new(source_count, each_size);
	if (rq == NULL){
		printf("Failed: nanorq_encoder_new\n");
		free(source_buf);
		free(parity_buf);
		return 1;
	}

/*
	// (Optionally) add every source data one by one
	p = source_buf;
	for (i = 0; i < source_count; i++){
		esi = i;
		rv = nanorq_encoder_add_symbol(rq, p, esi);
		if (rv != 0){
			printf("Failed: nanorq_encoder_add_symbol, %d, ESI = %d\n", rv, esi);
			free(source_buf);
			free(parity_buf);
			nanorq_free(rq);
			return 1;
		}
		p += each_size;
	}
*/
/*
	// (Optionally) precalculate precode matrix inversion
	// If this is omitted, it will be done at nanorq_generate_symbols().
	rv = nanorq_precalculate(rq);
	if (rv != 0){
		printf("Failed: nanorq_precalculate, %d\n", rv);
		free(source_buf);
		free(parity_buf);
		nanorq_free(rq);
		return 1;
	}
*/
/*
	// (Optionally) generate intermediate data from source data at once
	// If this is omitted, it will be done in nanorq_encode_all() or nanorq_encode().
	rv = nanorq_generate_symbols(rq, source_buf);
	if (rv != 0){
		printf("Failed: nanorq_generate_symbols, %d\n", rv);
		free(source_buf);
		free(parity_buf);
		nanorq_free(rq);
		return 1;
	}
*/

	// create all parity data at once
	rv = nanorq_encode_all(rq, parity_buf, parity_count, source_buf);
	if (rv != 0){
		printf("Failed: nanorq_encode_all, %d\n", rv);
		free(source_buf);
		free(parity_buf);
		nanorq_free(rq);
		return 1;
	}

/*
	// create every parity data one by one
	p = parity_buf;
	for (i = 0; i < parity_count; i++){
		esi = source_count + i;
		rv = nanorq_encode(rq, p, esi, source_buf);
		if (rv != 0){
			printf("Failed: nanorq_encode, %d, ESI = %d\n", rv, esi);
			free(source_buf);
			free(parity_buf);
			nanorq_free(rq);
			return 1;
		}
		p += each_size;
	}
*/

	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf(" ... %.3f seconds\n", duration);

	nanorq_free(rq);	// No need encoder

	output_size = (size_t)each_size * parity_count;
	value_mb = (double)source_size / 1000000.f;
	input_speed = value_mb / duration;
	output_speed = (double)output_size / 1000000.f / duration;
	printf("nanoRQc Encoder(%g MB in %d pieces, %d parities): Input= %g MB/s, Output= %g MB/s\n",
			value_mb, source_count, parity_count, input_speed, output_speed);

/*
{
	FILE *fp;
	if (fopen_s(&fp, "sourceN.bin", "wb") == 0){
		fwrite(source_buf, 1, source_size, fp);
		fclose(fp);
	}
	if (fopen_s(&fp, "parityN.bin", "wb") == 0){
		fwrite(parity_buf, 1, parity_size, fp);
		fclose(fp);
	}
}
*/

	start = clock();
	printf("\n set %d random error and allocate memory ...\n", lost_count);

	// shuffle order of source and parity data
	malloc_size = (source_count + parity_count) * sizeof(int);	// available source and parity
	order_buf = malloc(malloc_size);
	if (order_buf == NULL){
		printf("Failed: malloc, %zd\n", malloc_size);
		free(source_buf);
		free(parity_buf);
		return 1;
	}
	max = source_count + parity_count;
	for (i = 0; i < max; i++)
		order_buf[i] = i;

	// swap index to select available data
	pcg32_srandom(0, 1);
	// Fisher–Yates shuffle [0 ~ source_count - 1]
	i = source_count - 1;
	while (i > 0){
		j = pcg32_boundedrand(i + 1);
		rv = order_buf[j];
		order_buf[j] = order_buf[i];
		order_buf[i] = rv;
		i--;
	}
	// shuffle [source_count ~ source_count + lost_count - 1]
	min = source_count;
	i = source_count + lost_count - 1;
	while (i > min){
		j = min + pcg32_boundedrand(i - min + 1);
		rv = order_buf[j];
		order_buf[j] = order_buf[i];
		order_buf[i] = rv;
		i--;
	}
	if (parity_count > lost_count){
		// shuffle [source_count + lost_count ~ source_count + parity_count - 1]
		min = source_count + lost_count;
		i = source_count + parity_count - 1;
		while (i > min){
			j = min + pcg32_boundedrand(i - min + 1);
			rv = order_buf[j];
			order_buf[j] = order_buf[i];
			order_buf[i] = rv;
			i--;
		}
	}

	// allocate working buffer
	// nanorq decoder writes recovered source data on this.
	work_buf = calloc(source_size, 1);
	if (parity_buf == NULL){
		printf("Failed: malloc, %zd\n", source_size);
		free(source_buf);
		free(parity_buf);
		free(order_buf);
		return 1;
	}
	value_kb = (double)source_size / 1000;
	value_mb = (double)source_size / 1000000.f;
	printf("work data size = %zd bytes = %g KB = %g MB\n", source_size, value_kb, value_mb);

	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf(" ... %.3f seconds\n", duration);

	printf("\n decode ...\n");
	start = clock();

	// initialize nanorq decoder
	rq = nanorq_decoder_new(source_count, each_size, source_count + parity_count - 1);
	if (rq == NULL){
		printf("Failed: nanorq_decoder_new\n");
		free(source_buf);
		free(parity_buf);
		free(order_buf);
		free(work_buf);
		return 1;
	}

	// input available source data at first
	max = source_count - lost_count;
	for (i = 0; i < max; i++){
		esi = order_buf[i];	// index of available source data
		//rv = nanorq_decoder_add_symbol(rq, source_buf + (size_t)each_size * esi, esi, work_buf);
		// Set NULL, when no need to copy source data.
		rv = nanorq_decoder_add_symbol(rq, source_buf + (size_t)each_size * esi, esi, NULL);
		if (rv == NANORQ_SYM_ERR){
			printf("Failed: nanorq_decoder_add_symbol, %d, ESI = %d\n", rv, esi);
			free(source_buf);
			free(parity_buf);
			free(order_buf);
			free(work_buf);
			nanorq_free(rq);
			return 1;
		}
	}
	printf("input %d available source data\n", max);

	input_count = 0;
	// input available parity data until possible repair
	for (i = 0; i < parity_count; i++){
		input_count++;
		esi = order_buf[source_count + i];	// index of available parity data
		rv = nanorq_decoder_add_symbol(rq, parity_buf + (size_t)each_size * (esi - source_count), esi, NULL);
		if (rv == NANORQ_SYM_ERR){
			printf("Failed: nanorq_decoder_add_symbol, %d, ESI = %d\n", rv, esi);
			free(source_buf);
			free(parity_buf);
			free(order_buf);
			free(work_buf);
			nanorq_free(rq);
			return 1;
		}
		if (input_count >= lost_count){	// test possible repair
			rv = nanorq_recover_symbols(rq);
			if (rv == 0){	// got enough data to recover
				break;
			} else if (rv == NANORQ_ENC_ERR){
				printf("Failed: nanorq_recover_symbols, %d\n", rv);
			}
		}
	}
	printf("input %d available parity data\n", input_count);

	// fail at: s400 p123, s500 p145, s600 p255
	if (rv != 0){
		printf("Failed: need more parity to recover, %d\n", input_count);
		free(source_buf);
		free(parity_buf);
		free(order_buf);
		free(work_buf);
		nanorq_free(rq);
		return 1;
	}

	// recover all lost source data at once
	rv = nanorq_repair_all(rq, work_buf);
	if (rv != 0){
		printf("Failed: nanorq_repair_all, %d\n", rv);
		free(source_buf);
		free(parity_buf);
		free(order_buf);
		free(work_buf);
		nanorq_free(rq);
		return 1;
	}
/*
	// restore every source data one by one
	p = work_buf;
	for (i = 0; i < source_count; i++){
		rv = nanorq_restore(rq, p, i);
		if (rv != 0){
			printf("Failed: nanorq_restore, %d, ESI = %d\n", rv, i);
			free(source_buf);
			free(parity_buf);
			free(order_buf);
			free(work_buf);
			nanorq_free(rq);
			return 1;
		}
		p += each_size;
	}
*/
/*
	// recover every lost source data one by one
	for (i = 0; i < lost_count; i++){
		esi = order_buf[i + source_count - lost_count];	// index of lost source data
		rv = nanorq_restore(rq, work_buf + (size_t)each_size * esi, esi);
		if (rv != 0){
			printf("Failed: nanorq_restore, %d, ESI = %d\n", rv, esi);
			free(source_buf);
			free(parity_buf);
			free(order_buf);
			free(work_buf);
			nanorq_free(rq);
			return 1;
		}
	}
*/

	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf(" ... %.3f seconds\n", duration);

	free(source_buf);
	free(parity_buf);
	nanorq_free(rq);	// No need decoder

	output_size = (size_t)each_size * lost_count;	// bytes
	value_kb = (double)output_size / 1000;
	value_mb = (double)output_size / 1000000.f;
	printf("lost data size = %zd bytes = %g KB = %g MB\n", output_size, value_kb, value_mb);
	output_speed = value_mb / duration;

	input_count += source_count - lost_count;
	input_size = (size_t)each_size * input_count;
	value_kb = (double)input_size / 1000;
	value_mb = (double)input_size / 1000000.f;
	printf("input data size = %zd bytes = %g KB = %g MB\n", input_size, value_kb, value_mb);
	input_speed = value_mb / duration;
	printf("nanoRQc Decoder(%g MB in %d pieces, %d losses): Input= %g MB/s, Output= %g MB/s\n",
			value_mb, input_count, lost_count, input_speed, output_speed);

/*
{
	FILE *fp;
	if (fopen_s(&fp, "workN.bin", "wb") == 0){
		fwrite(work_buf, 1, source_size, fp);
		fclose(fp);
	}
}
*/

	start = clock();
	printf("\n verify recovered data ...\n");

	//lost_count = source_count;	// If you want to check all source data, enable this line.
	max = each_size / 4;
	for (i = 0; i < lost_count; i++){
		esi = order_buf[i + source_count - lost_count];	// index of lost source data
		int_p = (int *)(work_buf + (size_t)each_size * esi);
		pcg32_srandom(0, esi);
		for (j = 0; j < max; j++){
			if (int_p[j] != pcg32_random())
				break;
		}
		if (each_size & 3){	// each_size may not be multiple of 4.
			rv = pcg32_random();
			if (each_size % 4 == 1){
				if (work_buf[(size_t)each_size * esi + (each_size & ~3)] != (char)rv)
					j = -1;
			} else if (each_size % 4 == 2){
				if (work_buf[(size_t)each_size * esi + (each_size & ~3)] != (char)rv){
					j = -1;
				} else if (work_buf[(size_t)each_size * esi + (each_size & ~3) + 1] != (char)(rv >> 8)){
					j = -2;
				}
			} else if (each_size % 4 == 3){
				if (work_buf[(size_t)each_size * esi + (each_size & ~3)] != (char)rv){
					j = -1;
				} else if (work_buf[(size_t)each_size * esi + (each_size & ~3) + 1] != (char)(rv >> 8)){
					j = -2;
				} else if (work_buf[(size_t)each_size * esi + (each_size & ~3) + 2] != (char)(rv >> 16)){
					j = -3;
				}
			}
		}
		if (j < max){
			printf("Failed: recovered data, %d, ESI = %d\n", i, esi);
			break;
		}
	}
	if (i == lost_count)
		printf("All data repaired !\n");

	free(order_buf);	// End comparison
	free(work_buf);

	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf(" ... %.3f seconds\n", duration);

	return 0;
}

