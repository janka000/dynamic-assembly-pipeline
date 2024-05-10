#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include "sdict.h"
#include "paf.h"
#include "kvec.h"
#include "sys.h"
#include "miniasm.h"

#include "ksort.h"
#define ma_hit_key(a) ((a).qns)
KRADIX_SORT_INIT(hit, ma_hit_t, ma_hit_key, 8)

KSORT_INIT_GENERIC(uint32_t)

typedef kvec_t(uint32_t) uint32_v;

void ma_hit_sort(size_t n, ma_hit_t *a)
{
	radix_sort_hit(a, a + n);
}

void ma_hit_mark_unused(sdict_t *d, size_t n, const ma_hit_t *a)
{
	size_t i;
	for (i = 0; i < d->n_seq; ++i)
		d->seq[i].aux = 0;
	for (i = 0; i < n; ++i)
		d->seq[a[i].qns>>32].aux = d->seq[a[i].tn].aux = 1;
	for (i = 0; i < d->n_seq; ++i) {
		sd_seq_t *s = &d->seq[i];
		if (!s->aux) s->del = 1;
		else s->aux = 0;
	}
}

/**
 * find "contained" reads in the PAF file
 * used in step 0 of mianism main
 * @param const char *fn paf file name
 * @param int min_span (cmd param)
 * @param int min_match min match length (cmd param)
 * @param int max_hang	max over hang length (cmd param)
 * @param float int_frac min end-to-end match ratio (cmd param)
 * 
 * @return pointer to sdict_t with reads from PAF file that are "contained"

*/

sdict_t *ma_hit_no_cont(const char *fn, int min_span, int min_match, int max_hang, float int_frac)
{
	paf_file_t *fp;
	paf_rec_t r; //representation of a line from the PAF file
	sdict_t *d;

	fp = paf_open(fn);
	if (!fp) {
		fprintf(stderr, "[E::%s] could not open PAF file %s\n", __func__, fn);
		exit(1);
	}
	d = sd_init(); //initialize hash table
	while (paf_read(fp, &r) >= 0) { //load line from PAF to r
		int l5, l3;
		if (r.qe - r.qs < min_span || r.te - r.ts < min_span || r.ml < min_match) continue;
		l5 = r.rev? r.tl - r.te : r.ts; //5' end pos based on strand orientation
		l3 = r.rev? r.ts : r.tl - r.te; //3' end pos based on strand orientation
		if (r.ql>>1 > r.tl) {
			if (l5 > max_hang>>2 || l3 > max_hang>>2 || r.te - r.ts < r.tl * int_frac) continue; // internal match
			if ((int)r.qs - l5 > max_hang<<1 && (int)(r.ql - r.qe) - l3 > max_hang<<1)
				sd_put(d, r.tn, r.tl); //add the read to the hash table
		} else if (r.ql < r.tl>>1) {
			if (r.qs > max_hang>>2 || r.ql - r.qe > max_hang>>2 || r.qe - r.qs < r.ql * int_frac) continue; // internal
			if (l5 - (int)r.qs > max_hang<<1 && l3 - (int)(r.ql - r.qe) > max_hang<<1)
				sd_put(d, r.qn, r.ql); //add the read to the hash table
		}
	}
	paf_close(fp);
	if (ma_verbose >= 3) fprintf(stderr, "[M::%s::%s] dropped %d contained reads\n", __func__, sys_timestamp(), d->n_seq); //d->n_seq == # of reads in hash table
	return d; //hash table of contained reads
}

/**
 * creates (and returns) an array of MAPPINGS from PAF file, that are NOT CONTAINED, stores the READS form them into the hash table and counts the mappings stored
 * step 1 of minaism: read mappings
 * @param const char *fn   PAF file name
 * @param int min_span   min span (cmd param)
 * @param int min_match  min match length (cmd param)
 * @param sdict_t *d  hash table of READS
 * @param size_t *n  number of hits (pointer, to be set, the value will be used in the other functions)
 * @param int bi_dir   only one/both directions of an arc are present in input (cmd -b/-B)
 * @param const sdict_t *excl   table of contained reads - reads to be excluded
 * 
 * @return ma_hit_t *a - pointer to ma_hit_t* array -> MAPPINGS from PAF are stored there
 * but the results are: array of mappings, hash table of reads (d), hash table of contained reads (excl) and number of elements int the mappings array (n) 
*/
ma_hit_t *ma_hit_read(const char *fn, int min_span, int min_match, sdict_t *d, size_t *n, int bi_dir, const sdict_t *excl)
{
	paf_file_t *fp;
	paf_rec_t r;
	ma_hit_v h = {0,0,0}; //int m, n ; pointer to ma_hit_t (def. v miniasm.h), based on //https://github.com/attractivechaos/klib, kvec.h for use with klib
	size_t i, tot = 0, tot_len = 0;

	fp = paf_open(fn);
	if (!fp) {
		fprintf(stderr, "[E::%s] could not open PAF file %s\n", __func__, fn);
		exit(1);
	}
	while (paf_read(fp, &r) >= 0) { //iterate PAF file, load PAF entry into r
		ma_hit_t *p; // mame+start of query sequence  + other 8 ints from PAF..
		++tot;
		if (r.qe - r.qs < min_span || r.te - r.ts < min_span || r.ml < min_match) continue; //does not satisfy conditions given by params -> skip
		if (excl && (sd_get(excl, r.qn) >= 0 || sd_get(excl, r.tn) >= 0)) continue; //contained? -> skip
		kv_pushp(ma_hit_t, h, &p); //adds p to h //https://github.com/attractivechaos/klib, kvec.h = "generic dynamic array" ..kv_pushp(type,v,p)
		p->qns = (uint64_t)sd_put(d, r.qn, r.ql)<<32 | r.qs; //store the query sequence to the hash table, store query name + start packed into one int? to p->qns
		p->qe = r.qe; //store query end
		p->tn = sd_put(d, r.tn, r.tl); //store the target sequence into the hash table
		//strore the rest of values from paf
		//ts, te = target start, terget end, ml = num matches, rev = relative strand, bl = alignment block length
		p->ts = r.ts, p->te = r.te, p->rev = r.rev, p->ml = r.ml, p->bl = r.bl; 
		if (bi_dir && p->qns>>32 != p->tn) { //both directions of an arc are present in input and target!=query
			//store the mapping, same as above
			kv_pushp(ma_hit_t, h, &p);
			p->qns = (uint64_t)sd_put(d, r.tn, r.tl)<<32 | r.ts;
			p->qe = r.te;
			p->tn = sd_put(d, r.qn, r.ql);
			p->ts = r.qs, p->te = r.qe, p->rev = r.rev, p->ml = r.ml, p->bl = r.bl;
		}
	}
	paf_close(fp);
	for (i = 0; i < d->n_seq; ++i)
		tot_len += d->seq[i].len; //count the total length
	if (ma_verbose >= 3) fprintf(stderr, "[M::%s::%s] read %ld hits; stored %ld hits and %d sequences (%ld bp)\n", __func__, sys_timestamp(), tot, h.n, d->n_seq, tot_len);
	ma_hit_sort(h.n, h.a); //radix sort
	*n = h.n; //number of elements of the "generic dynamic array"
	return h.a; //pointer to elements of the "generic dynamic array" -> in our case array of values of type ma_hit_t (=PAF mappings)
}

/**
 * creates a filtering array sub to filter the mappings 
 * based on min coverage (using dynamic programming - checking existence of non-0 region with coverage at least min), 
 * based on identity (discards the ones that have match below some threshold or read aligned to itself); 
 * 
 * in step 2:
 * takes the array of mappings created by ma_hit_read (containments filtered) 
 * and aplies some simple further filtering of the mappings
 * end_ clip = 0
 * in step 3:
 * more precise filtering - end_clip = min_span/2 (min_span bases clipped when calculation coverage..)
 * 
 * @param int min_dp 	min covergae (cmd param)
 * @param float min_iden	min identity (cmd param) 
 * @param int end_clip  num of bases to be clipped off the end of the mapping when checking the coverage
 * @param size_t n   size of hit array from previous step (# not contained mappings)
 * @param const ma_hit_t *a  hit array from previous step (not contained mappings from paf)
 * @param size_t n_sub size of sub
 * 
 * @return ma_sub_t *sub - array with delete flags and longest regions with sufficient coverage
*/
ma_sub_t *ma_hit_sub(int min_dp, float min_iden, int end_clip, size_t n, const ma_hit_t *a, size_t n_sub)
{
	size_t i, j, last, n_remained = 0;
	
	//dynamic array structure similar to the one in ma_hit_read
	//  size_t n - current number of elements in the array
    //  size_t m - total capacity of the array
    //  int *a  pointer to the array 
	kvec_t(uint32_t) b = {0,0,0}; 
	
	ma_sub_t *sub = 0;

	sub = (ma_sub_t*)calloc(n_sub, sizeof(ma_sub_t)); //allocate the sub array - store regions with sufficent coverage, set delete flag if there is no such region
	for (i = 1, last = 0; i <= n; ++i) { //iterate over the mappings array
		//qns ==> query name + start packed into one int (?)
		// check if the current element is the last one 
		// or the query sequence id changed (we came to a new query sequence)
		if (i == n || a[i].qns>>32 != a[i-1].qns>>32) { 
			size_t start = 0; // track start position of the region
			int dp, qid = a[i-1].qns>>32; //dp for dynamic programming, qid for the query identifier based on the previous element's query sequence
			ma_sub_t max, max2; //structures - type ma_sub_t - used to store information about selected regions
			kv_resize(uint32_t, b, i - last); // resize the dynamic array b to accommodate the difference between the current index (i) and the last index (last)
			b.n = 0; //number of elements in the array is 0
			for (j = last; j < i; ++j) { // collect all starts and ends
				uint32_t qs, qe;
				if (a[j].tn == qid || a[j].ml < a[j].bl * min_iden) continue; // skip self match and low idednitiy mappings (match length < block_length*min identity)
				qs = (uint32_t)a[j].qns + end_clip, qe = a[j].qe - end_clip; // calculate start and end positions considering the end clipping
				if (qe > qs) { // the read survived clipping.. 
					//push start and end positions to array b
					kv_push(uint32_t, b, qs<<1);
					kv_push(uint32_t, b, qe<<1|1);
				}
			}
			ks_introsort_uint32_t(b.n, b.a); //sort the array (some kind of quicksort, heapsort, and insert sort mix algorithm) - params: num of elements to be sorted and the array
			max.s = max.e = max.del = max2.s = max2.e = max2.del = 0; //start, end, deletion flag for 2 longest matches
			//dynamic programming to identify contiguous regions in the query sequence where coverage is consistently above the specified threshold
			//dp - coverage at the current position of quary sequence, old_dp - previous position
			for (j = 0, dp = 0; j < b.n; ++j) { //loop over the soreted b array 
				int old_dp = dp;
				if (b.a[j]&1) --dp; else ++dp; //update dp based on parity of the element in b
				if (old_dp < min_dp && dp >= min_dp) { //there is a new region that has sufficient coverage (transition from insuffiecnt to sufficient coverage)
					start = b.a[j]>>1; //update the start position
				} else if (old_dp >= min_dp && dp < min_dp) { //region with suffiecnt coverage ends (transition from suffiecnt to insufficient coverage)
					int len = (b.a[j]>>1) - start; //calculate length of the region
					if (len > max.e - max.s) max2 = max, max.s = start, max.e = b.a[j]>>1; //if the region is longer than max, max2 = max and max = current
					else if (len > max2.e - max2.s) max2.s = start, max2.e = b.a[j]>>1; // the region is longer than max2 (but not max) -> max2 = current
				}
			}
			if (max.e - max.s > 0) { //longest region is longer than 0
				assert(qid < n_sub);
				sub[qid].s = max.s - end_clip;
				sub[qid].e = max.e + end_clip;
				sub[qid].del = 0;
				++n_remained;
			} else {
				sub[qid].del = 1; //if the max region is shorther than 0, set the deletion flag
			}
			last = i; //move the index and continue to the next mapping.. 
		}
	}
	free(b.a);
	if (ma_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] %ld query sequences remain after sub\n", __func__, sys_timestamp(), n_remained);
	return sub; //return the array with delete flags and longest regions with sufficient coverage
}

/**
 * using the array from ma_hit_sub, 
 * skip hits associated with deleted regions.
 * ensure that the coordinates (qs, qe, ts, te) are relative to their respective regions in the query and target sequences, considering the orientation of the alignment on the strands
 * check if the alignment satisfies the minimum span criteria
 * update the hit positions and store the hit in the result array if it meets the criteria
 * @param const ma_sub_t *reg  "sub" array - the regions with sufficient coverage
 * @param int min_span min span (cmd param)
 * @param size_t n num of mappings
 * @param ma_hit_t *a  array of mappings
 * 
 * @returns num of hits
 * modifies array of hits ma_hit_t *a
*/
size_t ma_hit_cut(const ma_sub_t *reg, int min_span, size_t n, ma_hit_t *a)
{
	size_t i, m;
	for (i = m = 0; i < n; ++i) {
		ma_hit_t *p = &a[i]; //store current mapping to p
		const ma_sub_t *rq = &reg[p->qns>>32], *rt = &reg[p->tn]; //correspondign regions for query and target
		int qs, qe, ts, te;
		if (rq->del || rt->del) continue; //skip deleted regions
		//calculate new start and end positions after cutting (based on the regions)
		//adjust the coordinates to reflect the correct positions on the query and target sequences
		//handling cases where the alignment extends beyond the ends of the target or query regions
		if (p->rev) { //if the alignemnt is on reverse strand
			qs = p->te < rt->e? (uint32_t)p->qns : (uint32_t)p->qns + (p->te - rt->e); //query start
			qe = p->ts > rt->s? p->qe : p->qe - (rt->s - p->ts); //query end
			ts = p->qe < rq->e? p->ts : p->ts + (p->qe - rq->e); //target start
			te = (uint32_t)p->qns > rq->s? p->te : p->te - (rq->s - (uint32_t)p->qns); //target end
		} else {
			qs = p->ts > rt->s? (uint32_t)p->qns : (uint32_t)p->qns + (rt->s - p->ts);
			qe = p->te < rt->e? p->qe : p->qe - (p->te - rt->e);
			ts = (uint32_t)p->qns > rq->s? p->ts : p->ts + (rq->s - (uint32_t)p->qns);
			te = p->qe < rq->e? p->te : p->te - (p->qe - rq->e);
		}
		//adjuest the position to be relative to their corresponding regions
		qs = (qs > rq->s? qs : rq->s) - rq->s;
		qe = (qe < rq->e? qe : rq->e) - rq->s;
		ts = (ts > rt->s? ts : rt->s) - rt->s;
		te = (te < rt->e? te : rt->e) - rt->s;
		// check if the alignment satisfies the minimum span criteria
		if (qe - qs >= min_span && te - ts >= min_span) {
			// update the hit positions and store the hit in the result array if it meets the criteria
			p->qns = p->qns>>32<<32 | qs, p->qe = qe, p->ts = ts, p->te = te;
			a[m++] = *p;
		}
	}
	if (ma_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] %ld hits remain after cut\n", __func__, sys_timestamp(), m);
	return m;
}

/**
 * filtering - convert hits to arcs, check for containments, keep overlaps only, exclude hits associated with deleted regions
 * @param const ma_sub_t *sub
 * @param max_hang max over hang length
 * @param min_ovlp min overlap
 * @param size_t n
 * @param ma_hit_t *a
 * @param float *cov pointer where to store coverage..  
 * 
 * @returns m number of reamining hits after filtering
 * coverage in *cov
*/
size_t ma_hit_flt(const ma_sub_t *sub, int max_hang, int min_ovlp, size_t n, ma_hit_t *a, float *cov)
{
	size_t i, m;  // i loop variable, m number of filtered hits
	asg_arc_t t; // struct definned in asg.h (ul, v, ol, del) uints - storing the arcs of assebly graph (ma_hit2arc takes pointer of this type to return arc, unused here.. )
	uint64_t tot_dp = 0, tot_len = 0; // total coverage and total length
	for (i = m = 0; i < n; ++i) { //loop over the input hits.. 
		ma_hit_t *h = &a[i]; //current hit
		const ma_sub_t *sq = &sub[h->qns>>32], *st = &sub[h->tn]; //corresponding query and target subregion
		int r;
		if (sq->del || st->del) continue; //skip hits associated with deleted regions
		// ma_hit2arc defined in miniasm.h
		// if r < 0 => filter: internal/query or target contained/short overlap (retrun value < 0) = MA_HT_INT/MA_HT_QCONT/MA_HT_TCONT/MA_HT_SHORT_OVLP = -1/../-4
		// if r >= 0 => it is the length of the overlap
		r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, max_hang, .5, min_ovlp, &t);
		if (r >= 0 || r == MA_HT_QCONT || r == MA_HT_TCONT)
			a[m++] = *h, tot_dp += r >= 0? r : r == MA_HT_QCONT? sq->e - sq->s : st->e - st->s;
	}
	//calculate the total length of the regions that contributed to the coverage
	for (i = 1; i <= m; ++i) 
		if (i == m || a[i].qns>>32 != a[i-1].qns>>32) //new sequence
			tot_len += sub[a[i-1].qns>>32].e - sub[a[i-1].qns>>32].s;
	//calculate crude coverage after filtering
	*cov = (double)tot_dp / tot_len;
	if (ma_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] %ld hits remain after filtering; crude coverage after filtering: %.2f\n", __func__, sys_timestamp(), m, *cov);
	return m;
}

/**
 * merges two sub arrays (results from ma_hit_sub.. ) 
 * probably want to merge two intervals in a way so that the result is inside both of them??
 * HOW DOES THIS WORK??!!
 * 
*/
void ma_sub_merge(size_t n_sub, ma_sub_t *a, const ma_sub_t *b)
{
	size_t i;
	for (i = 0; i < n_sub; ++i){
		a[i].e = a[i].s + b[i].e; //this cannot be inisde b[i], unless a[i].s<=0 ??!!! 
		a[i].s += b[i].s; //this way, the start will be on the right from a[i] and b[i] 
	}
}


/**
 * remove containments
 * 
 * @returns m number of hits remaining after containments removal
 * modifies ma_sub_t sub, ma_hit_t a, sdict_t *d
*/
size_t ma_hit_contained(const ma_opt_t *opt, sdict_t *d, ma_sub_t *sub, size_t n, ma_hit_t *a)
{
	int32_t *map, r;
	size_t i, m, old_n_seq = d->n_seq;
	asg_arc_t t;
	for (i = m = 0; i < n; ++i) { //loop over all hits
		ma_hit_t *h = &a[i];
		ma_sub_t *sq = &sub[h->qns>>32], *st = &sub[h->tn];
		//convert hit to arc
		r = ma_hit2arc(h, sq->e - sq->s, st->e - st->s, opt->max_hang, opt->int_frac, opt->min_ovlp, &t);
		//mark contained for deletion
		if (r == MA_HT_QCONT) sq->del = 1;
		else if (r == MA_HT_TCONT) st->del = 1;
	}
	//mark deleted sequences
	for (i = 0; i < d->n_seq; ++i)
		if (sub[i].del) d->seq[i].del = 1;

	ma_hit_mark_unused(d, n, a); //mark unused hits	
	map = sd_squeeze(d); //squeeze the sequence dictionary (in sdict.c) - removes the elements marked for deletion, returns map to new indices from old..
	for (i = 0; i < old_n_seq; ++i) //update the sub array using map
		if (map[i] >= 0) sub[map[i]] = sub[i];
	for (i = m = 0; i < n; ++i) { //loop over all hits
		ma_hit_t *h = &a[i];
		int32_t qn = map[h->qns>>32], tn = map[h->tn]; 
		if (qn >= 0 && tn >= 0) { //update, if it was not deleted (i.e. not contained)
			a[i].qns = (uint64_t)qn<<32 | (uint32_t)a[i].qns;
			a[i].tn = tn;
			a[m++] = a[i];
		}
	}
	free(map);
	if (ma_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] %d sequences and %ld hits remain after containment removal\n", __func__, sys_timestamp(), d->n_seq, m);
	return m;
}
