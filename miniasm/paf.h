#ifndef PAF_PAF_H
#define PAF_PAF_H

#include <stdint.h>
#include <sys/types.h>

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif

typedef struct { //PAF file pointer + buffer?
	void *fp;
	kstring_t buf;
} paf_file_t;

typedef struct {
	const char *qn, *tn; // these point to the input string; NOT allocated (query name, target name)
	uint32_t ql, qs, qe, tl, ts, te; //query length,start,end; target len,start,end  
	uint32_t ml:31, rev:1, bl; //num matches, relative strand, alignment block length
} paf_rec_t; //representation of a line from the PAF file

#ifdef __cplusplus
extern "C" {
#endif

paf_file_t *paf_open(const char *fn);
int paf_close(paf_file_t *pf);
int paf_read(paf_file_t *pf, paf_rec_t *r);

#ifdef __cplusplus
}
#endif

#endif
