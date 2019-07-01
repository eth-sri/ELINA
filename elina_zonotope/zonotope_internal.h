/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2019 Department of Computer Science, ETH Zurich
 *  This software is distributed under GNU Lesser General Public License Version 3.0.
 *  For more information, see the ELINA project website at:
 *  http://elina.ethz.ch
 *
 *  THE SOFTWARE IS PROVIDED "AS-IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER
 *  EXPRESS, IMPLIED OR STATUTORY, INCLUDING BUT NOT LIMITED TO ANY WARRANTY
 *  THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS OR BE ERROR-FREE AND ANY
 *  IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
 *  TITLE, OR NON-INFRINGEMENT.  IN NO EVENT SHALL ETH ZURICH BE LIABLE FOR ANY     
 *  DAMAGES, INCLUDING BUT NOT LIMITED TO DIRECT, INDIRECT,
 *  SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN
 *  ANY WAY CONNECTED WITH THIS SOFTWARE (WHETHER OR NOT BASED UPON WARRANTY,
 *  CONTRACT, TORT OR OTHERWISE).
 *
 */


#ifndef _ZONOTOPE_INTERNAL_H_
#define _ZONOTOPE_INTERNAL_H_

#include <time.h>
#include <limits.h>
#include <unistd.h>
#include <string.h>
#include  <fenv.h>
#include <errno.h>
#include "zonotope.h"
#include "rdtsc.h"
#include "elina_box_meetjoin.h"



#ifdef __cplusplus
extern "C" {
#endif

#define start_timing()				\
tsc_counter start, end;				\
double cycles;					\
CPUID();					\
RDTSC(start)
    
#define record_timing(counter)		\
RDTSC(end);				\
CPUID();				\
cycles = (double)(COUNTER_DIFF(end, start));	\
counter += cycles
    
    extern double zonotope_copy_time;
    extern double zonotope_is_equal_time;
    extern double zonotope_is_lequal_time;
    extern double zonotope_permute_dimension_time;
    extern double zonotope_add_dimension_time;
    extern double zonotope_remove_dimension_time;
    extern double zonotope_top_time;
    extern double zonotope_bottom_time;
    //extern double meet_time;
    extern double zonotope_join_time;
    //extern double widening_time;
    extern double zonotope_free_time;
    extern double zonotope_forget_array_time;
    extern double zonotope_meet_lincons_time;
    extern double zonotope_to_box_time;
    extern double zonotope_of_box_time;
    extern double zonotope_is_top_time;
    extern double zonotope_is_bottom_time;
    //extern double expand_time;
    //extern double fold_time;
    //extern double sat_lincons_time;
    extern double zonotope_assign_linexpr_time;
    //extern double substitute_linexpr_time;
    //extern double bound_dimension_time;
    //extern double opt_conversion_time;
    //extern long int join_count;
    //extern double poly_is_unconstrained_time;
    
typedef unsigned  int uint_t;

/**************************************************************************************************/
/* INTERNAL DATA TYPE */
/**************************************************************************************************/

/*****************/
/* Noise symbol */
/****************/
typedef enum noise_symbol_t {
    IN,		/* Classical noise symbol */
    UN,		/* Union noise symbol */
} noise_symbol_t;

/* Noise symbol type */
typedef struct zonotope_noise_symbol_t {
    noise_symbol_t	type;		/* type of noise symbol */	/* obsolete */
    uint_t	index;		/* global index, noise symbols of the same index are shared */
//    bool	constrained;	/* true if the noise symbol is constrained */
} zonotope_noise_symbol_t;

/*****************/
/* Zonotope term */ 
/*****************/

/* Zonotope affine arithmetic term */
typedef struct _zonotope_aaterm_t {
    struct _zonotope_aaterm_t*	n;	/* next element */
    zonotope_noise_symbol_t*	pnsym;	/* index of the noise symbol */
    //elina_interval_t*	coeff;	/* coeff, encoded as interval */
   double sup;
   double inf;
} zonotope_aaterm_t;

/************************/
/* Zonotope affine form */
/************************/
struct _zonotope_aff_t {
    double c_inf;	/* center */
    double c_sup;
    zonotope_aaterm_t*	q;	/* first center term (epsilons) aaterm */
    zonotope_aaterm_t*	end;	/* quick jump to the last center term : to add a new term for instance */
    unsigned long long int		l;	/* number of noise symbols */
    unsigned long long int		pby;	/* # pointers to this affine form */
    double 	itv_inf;	/* best known interval concretisation */
    double 	itv_sup;
};
typedef struct _zonotope_aff_t zonotope_aff_t;

    static inline zonotope_aff_t * zonotope_aff_copy(zonotope_aff_t *src){
        zonotope_aff_t *res = (zonotope_aff_t *)malloc(sizeof(zonotope_aff_t));
	res->c_inf = src->c_inf;
        res->c_sup = src->c_sup;
        res->q = src->q==NULL? NULL: (zonotope_aaterm_t *)malloc(sizeof(zonotope_aaterm_t *));
        res->end = src->end==NULL? NULL:(zonotope_aaterm_t *)malloc(sizeof(zonotope_aaterm_t *));
        res->l = src->l;
        res->itv_inf = src->itv_inf; 
	res->itv_sup = src->itv_sup; 
        return res;
    }

typedef struct _zonotope_internal_t {
    //zonotope_internal_t* itv;		/* interval internal representation */
    uint_t		dim;		/* nb of noise symbol used */
    zonotope_noise_symbol_t**	epsilon;	/* array of size index of epsilons */
    elina_funid_t	funid;		/* current function */
    elina_manager_t *	man;		/* back-pointer */
    elina_manager_t *	manNS;		/* abstract domain of noise symbols */
    elina_manager_t *	box;		/* box abstract domain used to compute constraints meet with an hypercube */
    elina_lincons0_array_t moo;		/* array of constraints -1 <= eps_i <= 1; size = 2index */
    elina_interval_t*	muu;		/* [-1,1] (elina_interval_t) */
    elina_interval_t*	ap_muu;		/* [-1,1] (type elina_interval) */
    zonotope_aff_t	*top;		/* top interval */
    zonotope_aff_t	*bot;		/* bottom interval */
    elina_dim_t	*dimtoremove;	/* array to store dimensions to remove after a join */
    elina_dimchange_t*	dimchange;
    elina_abstract0_t*	nsymhypercube;
    uint_t* inputns;
    uint_t epssize;
    uint_t it;
} zonotope_internal_t;

/***********/
/*** T1+ ***/
/***********/
typedef struct _zonotope_t {
    zonotope_aff_t**		paf;            /* array of pointers to Taylor1+ expressions of size dims */
    double *box_inf;  /* reduced product with boxes */
    double *box_sup;
    uint_t		intdim;         /* nb of integer variables */
    uint_t		dims;           /* intdim + realdim */
    elina_abstract0_t* 	abs;        	/* nsym abstract object (=contraints over noise symbols)*/
    elina_dim_t*	nsymcons;       /* array of index of constrained noise symbols */
    elina_interval_t**	gamma;		/* pointer to an array which contains the concretisations of constrained noise symbols if any */
    unsigned long long int		size;		/* size of nsymcons and gamma */
    bool		hypercube;	/* true if no constrained nsym */
    //elina_interval_t**	g;	/* array of the generators of the zonotope - a oublier */
    uint_t		gn;		/* size of generators - a oublier */
} zonotope_t;

/* special object to store and compute meet with lincons */
typedef struct _obj {
  //  uint_t index;
    elina_interval_t* itv;
    elina_interval_t* coeff;
} obj;

static inline void zonotope_aaterm_free(zonotope_internal_t* pr, zonotope_aaterm_t* term);

static inline void zonotope_aaterm_list_free(zonotope_internal_t *pr, zonotope_aaterm_t* head);

static inline void zonotope_aff_free(zonotope_internal_t *pr, zonotope_aff_t *a);

/* memory allocation of affine form term: \alpha_i^x\epsilon_i */
static inline zonotope_aaterm_t* zonotope_aaterm_alloc_init(void)
{
    zonotope_aaterm_t* res = (zonotope_aaterm_t*)malloc(sizeof(zonotope_aaterm_t));
    res->n = NULL;
    res->pnsym = NULL;
    res->inf = 0;
    res->sup = 0;
    return res;
}

static zonotope_aff_t* zonotope_aff_alloc_init(zonotope_internal_t *pr)
{
    zonotope_aff_t* a = (zonotope_aff_t*)malloc(sizeof(zonotope_aff_t));
    a->c_inf = 0.0;
    a->c_sup = 0.0;
    a->q = NULL;
    a->end = NULL;
    a->l = 0;
    a->pby = 0;
    a->itv_inf = INFINITY;
    a->itv_sup = INFINITY;
    return a;
}

static inline zonotope_aff_t * zonotope_aff_top(zonotope_internal_t* pr)
{
    zonotope_aff_t* res = zonotope_aff_alloc_init(pr);
    res->c_inf = INFINITY;
    res->c_sup = INFINITY; 
    res->itv_inf = INFINITY;
    res->itv_sup = INFINITY;
    return res;
}
static inline zonotope_aff_t * zonotope_aff_bottom(zonotope_internal_t* pr)
{
    zonotope_aff_t* res = zonotope_aff_alloc_init(pr);
    res->c_inf = -1;
    res->c_sup = -1;
    res->itv_inf = -1;
    res->itv_sup = -1;
    return res;
}

static inline zonotope_noise_symbol_t* zonotope_noise_symbol_add(zonotope_internal_t *pr, noise_symbol_t type)
    /* increment the global index of used noise symbols and add the noise symbol in pr->eps */
{
    uint_t dim = pr->dim;
    zonotope_noise_symbol_t* res;
	
    /* resize epsilon array */
    if (dim % 1024 == 0){ pr->epsilon = (zonotope_noise_symbol_t**)realloc(pr->epsilon, (dim+1024)*sizeof(zonotope_noise_symbol_t*));
}
    res = pr->epsilon[dim] = (zonotope_noise_symbol_t*)malloc(sizeof(zonotope_noise_symbol_t));
    if (type == IN) {
	if (pr->epssize % 1024 == 0) pr->inputns = (uint_t*)realloc(pr->inputns, (pr->epssize+1024)*sizeof(uint_t));
	pr->inputns[pr->epssize] = dim; pr->epssize++;
    }
    res->index = dim;
    res->type = type;
    pr->dim++;
    return res;
}

    
    static inline void zonotope_noise_symbol_fprint(FILE* stream, zonotope_noise_symbol_t *eps)
    {
        switch (eps->type) {
            case IN: fprintf(stream,"(eps%u)",eps->index); break;
            case UN: fprintf(stream,"(eta%u)",eps->index); break;
            default: fprintf(stderr,"error: unknown type of noise symbol, aborting...\n"); abort(); break;
        }
    }
    
    
    /* Pretty print an affine term */
    static inline void zonotope_aaterm_fprint(zonotope_internal_t *pr, FILE* stream, zonotope_aaterm_t *ptr)
    {
      if ((ptr->inf != INFINITY) && (ptr->inf == ptr->sup)) {
        // elina_scalar_fprint(stream,ptr->sup);
        fprintf(stream, "%.*g", elina_scalar_print_prec, ptr->sup + 0.0);
      } else {
        printf("[%.20f, %.20f]\n", -ptr->inf, ptr->sup);
      }
        fprintf(stream,".");
        zonotope_noise_symbol_fprint(stream,ptr->pnsym);
    }
    

    /* pretty print an Affine expression: expr/top/bottom */
    static inline void zonotope_aff_fprint(zonotope_internal_t* pr, FILE* stream, zonotope_aff_t *expr)
    {
        if(expr==NULL){
            return;
        }
        
        zonotope_aaterm_t* p;
        if ((expr->c_inf != INFINITY) && (expr->c_inf == expr->c_sup)) {
          fprintf(stream, "%.20f", expr->c_sup);
        } else {
          fprintf(stream,"[%.20f,%.20f]",-expr->c_inf,expr->c_sup);
        }
        /* Print values */
        for (p=expr->q; p; p=p->n) {
            fprintf(stream," + ");
            zonotope_aaterm_fprint(pr, stream, p);
        }
        fprintf(stream,"\t;");
        fprintf(stream,"[%.20f,%.20f]",-expr->itv_inf,expr->itv_sup);
        fprintf(stream,"\t");
        fflush(stream);
    }


/*******************************/
/* Zonotope internal structure */
/*******************************/
/* managing internal representation of zonotope abstract value */
static inline zonotope_internal_t* zonotope_internal_alloc(elina_manager_t* manNS);
/* free internal */
static inline void zonotope_internal_free(zonotope_internal_t* a);
/* Initializes some fields of internal t1p structure from manager */
zonotope_internal_t* zonotope_init_from_manager(elina_manager_t* man, elina_funid_t funid);

/* alloc Zonotope noise symbols manager */
elina_manager_t* zonotope_manager_alloc(void);

zonotope_aff_t * zonotope_aff_from_linexpr0(zonotope_internal_t* pr, elina_linexpr0_t * expr, zonotope_t *z);

/* Free memory used by one aaterm */
static inline void zonotope_aaterm_free(zonotope_internal_t* pr, zonotope_aaterm_t* term)
{
   // arg_assert(term, abort(););
    term->n = NULL;
    term->pnsym = NULL;
    free(term);
}
/* free memory used by a chained list starting from one term */
static inline void zonotope_aaterm_list_free(zonotope_internal_t *pr, zonotope_aaterm_t* head)
{
    zonotope_aaterm_t *p,*q;
    p = q = NULL;
    for (p = head; p; p = q) {
	q = p->n;
	zonotope_aaterm_free(pr, p);
    }
}

static inline void zonotope_aff_free(zonotope_internal_t *pr, zonotope_aff_t *a)
{
    if (a->pby) {
	
	//fatal("You are about to free a used affine form\n");
    } else {
	a->pby = 0;
	a->c_inf = 0.0;
	a->c_sup = 0.0;
	if (a->q) zonotope_aaterm_list_free(pr, a->q);
	a->q = NULL;
	a->end = NULL;
	a->l = (uint_t)0;
	a->itv_inf = 0;
	a->itv_sup = 0;
	//elina_interval_free(a->itv);
	free(a);
	a = NULL;
    }
}

static inline void zonotope_aff_check_free(zonotope_internal_t *pr, zonotope_aff_t *a)
{
    if(a==NULL){return;}
    if (a->pby>0) a->pby--;
    if (a->pby == 0) {
      if ((a != pr->top) && (a != pr->bot)) zonotope_aff_free(pr, a);
    }
  // if a is NULL, do nothing
}

static inline void zonotope_aff_noise_symbol_create(zonotope_internal_t *pr, zonotope_aff_t *expr, double coeff_inf, double coeff_sup, noise_symbol_t type)
{
    elina_interval_t *zero = elina_interval_alloc();
    elina_interval_t * coeff = elina_interval_alloc();
    elina_interval_set_double(coeff,coeff_inf,coeff_sup);
    if (elina_interval_cmp(coeff,zero)>=0) {
	zonotope_aaterm_t* ptr = zonotope_aaterm_alloc_init();
	ptr->inf = coeff_inf;
	ptr->sup = coeff_sup;
	ptr->pnsym = zonotope_noise_symbol_add(pr, type);
	if (expr->end) expr->end->n = ptr;
	else expr->q = ptr;
	expr->end = ptr;
	expr->l++;
    }
    elina_interval_free(zero);
    elina_interval_free(coeff);
}



static inline void elina_interval_middev(double *mid_inf, double *mid_sup, double *dev_inf, double *dev_sup, double a_inf, double a_sup)
{
    double tmp[4];
    for(int i=0; i < 4; i++){
	tmp[i] = 0.0;
    }
    
    if((a_inf!=INFINITY) && (-a_inf==a_sup)){
	*mid_inf = a_inf;
	*mid_sup = a_sup;
	*dev_inf = 0.0;
	*dev_sup = 0.0;
    } else if ((a_sup==INFINITY) || (a_inf==INFINITY) || (-a_inf>a_sup)) {
	*mid_inf = INFINITY;
	*mid_sup = INFINITY;
        *dev_inf = INFINITY;
	*dev_sup = INFINITY;
    } else {
	/* a = [x,y], x < y,
	 * tmp[0] = x+y */
	/* Rounding to Nearest, only x86 Unix like */
	tmp[0] = a_sup - a_inf;
	if (!fesetround(FE_TONEAREST)) {
	    /* tmp[1] = (x+y)/2 -- using ldexp if double */
	    tmp[1] = tmp[0]/2;
	    /* mid = [tmp[1], tmp[1]] */
	    *mid_sup = tmp[1];
	    *mid_inf = -tmp[1];
	    
	}
	else {
	    fprintf(stderr,"fesetround: %s\n", strerror (errno));
	    abort();
	}

	/* Come back to default rounding mode */
	if (fesetround(FE_UPWARD)) {
	    fprintf(stderr,"fesetround: %s\n", strerror (errno));
	    abort();
	} else {
	    /* tmp[0] = (x+y)/2 - x */
	    tmp[0] = tmp[1] + a_inf;
	    /* tmp[2] = y - (x+y)/2 */
	    tmp[2] = a_sup - tmp[1];
	    /* tmp[3] = max(tmp[0],tmp[2]) -- fmax if double */
	    tmp[3] = fmax(tmp[0],tmp[2]);
	    /* dev = [tmp[3], tmp[3]] */
	    *dev_sup = tmp[3];
	    *dev_inf = -tmp[3];
	}
    }
    
}



/* convert an itv to an affine form with a fresh noise symbol then add this form to the affine form expr */
static inline void zonotope_aff_add_itv(zonotope_internal_t* pr, zonotope_aff_t *expr, double itv_inf, double itv_sup, noise_symbol_t type)
{
	
    /* itv is a non point interval with finite bounds */
    //elina_interval_t *mid, *dev;
    //mid = elina_interval_alloc();
    //dev = elina_interval_alloc();
    double mid_inf = 0.0;
    double mid_sup = 0.0;
    double dev_inf = 0.0;
    double dev_sup = 0.0;
    if(isfinite(itv_inf)&&(-itv_inf!=itv_sup)){
	elina_interval_middev(&mid_inf, &mid_sup, &dev_inf, &dev_sup, itv_inf,itv_sup);
	//printf("mid\n");
	//elina_scalar_print(mid);
	//printf("\n");
	expr->c_inf = expr->c_inf + mid_inf;
	expr->c_sup = expr->c_sup + mid_sup;
	zonotope_aff_noise_symbol_create(pr, expr, dev_inf, dev_sup, type);
	
    } else{
	 expr->c_inf = expr->c_inf + itv_inf;
	 expr->c_sup = expr->c_sup + itv_sup;
    }
    
	
}

static inline bool findKR(uint_t *res, uint_t x, uint_t* tab, uint_t size)
{
    if(size<=0){
	return false;
    }
    
    int low = 0;
    int high = size - 1;
    int mid = -1;
    while (low <= high) {
	mid = (high + low) / 2;
	if (x < tab[mid]) high = mid - 1;
	else if (x > tab[mid]) low = mid + 1;
	else {
	    *res = mid;
	    return true;
	}
    }
    *res = low;
    return false;
}




static inline uint_t zonotope_noise_symbol_cons_get_dimension(zonotope_internal_t *pr, zonotope_t *z)
{
    elina_dimension_t dimension = elina_abstract0_dimension(pr->manNS, z->abs);
    return (dimension.intdim + dimension.realdim);
}

static inline bool zonotope_noise_symbol_cons_get_dimpos(zonotope_internal_t* pr, elina_dim_t* dim, uint_t nsymIndex, zonotope_t* z)
{
    uint_t size = zonotope_noise_symbol_cons_get_dimension(pr, z);
    if (size) return (elina_dim_t)findKR((uint_t*)dim, (uint_t)nsymIndex, (uint_t*)z->nsymcons, size);
    else return false;
}

static inline double min_double(double a, double b){
	return a < b ? a : b;
}

static inline double max_double(double a, double b){
	return a > b ? a : b;
}


static inline void zonotope_aff_bound(zonotope_internal_t* pr, double *res_inf, double *res_sup, zonotope_aff_t *expr, zonotope_t* z)
{
    if ((expr->c_inf==INFINITY) && (expr->c_sup==INFINITY)) {
        *res_inf = INFINITY;
        *res_sup = INFINITY;
	return;
    } else {
	elina_dim_t dim;
	double tmp_inf = 0.0;
	double tmp_sup = 0.0;
	elina_interval_t *eps_itv = elina_interval_alloc();
	zonotope_aaterm_t* p;
	*res_inf = expr->c_inf;
        *res_sup = expr->c_sup;
	if (z->hypercube) {
	    for (p=expr->q; p; p=p->n) {
		elina_double_interval_mul(&tmp_inf, &tmp_sup,p->inf, p->sup, pr->muu->inf->val.dbl,pr->muu->sup->val.dbl);
		*res_inf = *res_inf + tmp_inf;
		*res_sup = *res_sup + tmp_sup;
	    }
		
	} else {
	    elina_linexpr0_t* linexpr0 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE, 0);
	    elina_linexpr0_set_cst_scalar_double(linexpr0, 0.0);
	    
	    linexpr0->p.linterm = (elina_linterm_t*)malloc(expr->l*sizeof(elina_linterm_t));
	    uint_t k = 0;
	    elina_dim_t dim = 0;
	    for (p=expr->q; p; p=p->n) {
		if (zonotope_noise_symbol_cons_get_dimpos(pr, &dim, p->pnsym->index, z)) {
		    elina_coeff_init(&linexpr0->p.linterm[k].coeff, ELINA_COEFF_INTERVAL);
		    elina_coeff_set_interval_double(&linexpr0->p.linterm[k].coeff, -p->inf,p->sup);
		    linexpr0->p.linterm[k].dim = dim;
		    k++;
		} else {
		    elina_double_interval_mul(&tmp_inf, &tmp_sup, p->inf, p->sup, pr->muu->inf->val.dbl, pr->muu->sup->val.dbl);
		    *res_inf = *res_inf + tmp_inf;
		    *res_sup = *res_sup + tmp_sup;
		}
	    }
		
	    linexpr0->size = k;
	    elina_interval_t* elina_itv = elina_abstract0_bound_linexpr(pr->manNS, z->abs, linexpr0);
	    tmp_inf = -elina_itv->inf->val.dbl;
	    tmp_sup = elina_itv->sup->val.dbl;
	    *res_inf = *res_inf + tmp_inf;
	    *res_sup = *res_sup + tmp_sup;
		
	    linexpr0->p.linterm = (elina_linterm_t*)realloc((void* )linexpr0->p.linterm, k*sizeof(elina_linterm_t));
	    elina_linexpr0_free(linexpr0);
	    elina_interval_free(elina_itv);
	}
	elina_interval_free (eps_itv);
    }
}

static inline bool zonotope_aff_is_known_to_be_zero(zonotope_internal_t *pr, zonotope_aff_t *a)
{
    if ((!a->itv_inf)&& (!a->itv_sup)) return true;
    else return false;
}

static inline bool zonotope_aff_is_top(zonotope_internal_t* pr, zonotope_aff_t *a)
{
    if (a == pr->top) return true;
    else if ((a->c_inf!=INFINITY) || (a->c_sup!=INFINITY)) return false;
    else if ((a->itv_inf!=INFINITY)||(a->itv_sup!=INFINITY)) return false;
    else if (a->q != NULL) return false;
    else return true;
}
static inline bool zonotope_aff_is_bottom(zonotope_internal_t* pr, zonotope_aff_t *a)
{
    if (a == pr->bot) return true;
    else if ((-a->c_inf<=a->c_sup)) return false;
    else if (-a->itv_inf<=a->itv_sup) return false;
    else if (a->q != NULL) return false;
    else return true;
}



zonotope_aff_t* zonotope_aff_mul_itv(zonotope_internal_t* pr, zonotope_aff_t* src, elina_interval_t *lambda);
static inline zonotope_aff_t* zonotope_aff_mul_scalar(zonotope_internal_t* pr, zonotope_aff_t* src, elina_scalar_t *lambda)
    {
        if ((!elina_scalar_sgn(lambda) )|| zonotope_aff_is_known_to_be_zero(pr, src)) {
		zonotope_aff_t *res = zonotope_aff_alloc_init(pr);
		res->itv_inf = 0;
		res->itv_sup = 0;
            	return res;
        } else if (zonotope_aff_is_bottom(pr, src)) {
		
            return zonotope_aff_bottom(pr);
        } else if (zonotope_aff_is_top(pr, src)) {
		
            return zonotope_aff_top(pr);
        } else {
            int sgn = elina_scalar_sgn(lambda);
            zonotope_aff_t* dst = NULL;
            zonotope_aaterm_t *p,*q;
            q = NULL;
            dst = zonotope_aff_alloc_init(pr);
            if(sgn>=0){
		dst->c_inf = src->c_inf*lambda->val.dbl;
		dst->c_sup = src->c_sup*lambda->val.dbl;
            }
            else{
                //elina_scalar_t * add = elina_scalar_alloc_set(tmp->sup);
		dst->c_sup = src->c_inf*-lambda->val.dbl;
		dst->c_inf = src->c_sup*-lambda->val.dbl;
                //elina_scalar_free(add);
            }
            
            if (src->q) {
                dst->q = q = zonotope_aaterm_alloc_init();
                for (p=src->q; p; p=p->n) {
                    if(sgn>=0){
			q->inf = p->inf * lambda->val.dbl;
			q->sup = p->sup * lambda->val.dbl;
                        //elina_scalar_mul(q->inf,p->inf,lambda,ELINA_SCALAR_DOUBLE);
                        //elina_scalar_mul(q->sup,p->sup,lambda,ELINA_SCALAR_DOUBLE);
                    }
                    else{
			q->sup = p->inf * -lambda->val.dbl;
			q->inf = p->sup * -lambda->val.dbl;
                        //elina_scalar_t * add = elina_scalar_alloc_set(tmp->sup);
                        //elina_scalar_mul(q->sup,p->inf,lambda,ELINA_SCALAR_DOUBLE);
                        //elina_scalar_mul(q->inf,p->sup,lambda,ELINA_SCALAR_DOUBLE);
                        //elina_scalar_free(add);
                    }
                    
                    q->pnsym = p->pnsym;
                    if (p->n) {
                        /* continue */
                        q->n = zonotope_aaterm_alloc_init();
                        q = q->n;
                    } else {
                        /* the last iteration */
                        dst->end = q; 
                    }
                }
                
            }
            
            dst->l = src->l;
            if(sgn>=0){
		dst->itv_inf = src->itv_inf*lambda->val.dbl;
                dst->itv_sup = src->itv_sup*lambda->val.dbl;
		
            }
            else{
                //elina_scalar_t * add = elina_scalar_alloc_set(tmp->sup);
		dst->itv_sup = src->itv_inf*-lambda->val.dbl;
                dst->itv_inf = src->itv_sup*-lambda->val.dbl;
                
                //elina_scalar_free(add);
            }
            
            
            return dst;
        }
    }

    
static inline void elina_interval_abs_double(elina_interval_t *a, double inf, double sup)
{
  if (inf<=0){
    /* positive interval */
    elina_interval_set_double(a,-inf, sup);
  }
  else if (sup<=0){
    /* negative interval */
    elina_interval_set_double(a,-sup, inf);
  }
  else {
    a->sup->val.dbl = fmax(inf, sup);
    a->inf->val.dbl = 0;
   // elina_scalar_max(a->sup,b->inf,b->sup);
    //elina_scalar_set_to_int(a->inf,0,discr);
  }
}    

static inline void zonotope_aff_cons_eq_lambda(zonotope_internal_t* pr, elina_interval_t** res, zonotope_aff_t* x, zonotope_aff_t* cons, zonotope_t *z)
{
    zonotope_aaterm_t *p, *q;
    p = q = NULL;
   
    obj** array = NULL;
    if (cons->l + x->l > pr->dim)  array = (obj**)calloc(pr->dim,sizeof(obj*)); 
    else array = (obj**)calloc((cons->l + x->l),sizeof(obj*)); 
    uint_t max = 0;
    uint_t i = 0;
    elina_interval_t *tmp = elina_interval_alloc(); 
    //elina_interval_t *mid = elina_interval_alloc();
    elina_interval_t *dev = elina_interval_alloc();
    double mid_inf = 0.0;
    double mid_sup = 0.0;
    double dev_inf = 0.0;
    double dev_sup = 0.0;
    elina_interval_t *dim_itv = elina_interval_alloc();
    elina_dim_t dim;
    
    for (p=cons->q, q=x->q; p || q;) {
	if (p && q) {
	    if (p->pnsym->index == q->pnsym->index) {
		
		if((p->inf!=0) || (p->sup!=0)) {
		    array[i] = (obj*)calloc(1,sizeof(obj));
		    array[i]->coeff = elina_interval_alloc();
		    array[i]->itv = elina_interval_alloc();
		    if (z->hypercube) {
			elina_interval_abs_double(tmp, p->inf, p->sup);
			elina_interval_set(array[i]->coeff,tmp);
		    } else {
			if (zonotope_noise_symbol_cons_get_dimpos(pr, &dim, p->pnsym->index, z)) {
			    elina_interval_set(dim_itv, z->gamma[dim]);
			    elina_interval_middev(&mid_inf, &mid_sup, &dev_inf,&dev_sup, dim_itv->inf->val.dbl,dim_itv->sup->val.dbl );
			    elina_interval_abs_double(tmp, p->inf, p->sup);
			    elina_interval_set_double(dev,dev_inf,dev_sup);
			    elina_interval_mul(array[i]->coeff, tmp, dev,ELINA_SCALAR_DOUBLE);
			} else {
			    elina_interval_abs_double(tmp, p->inf, p->sup);
			    elina_interval_set(array[i]->coeff,tmp);
			}
		    }
		    elina_interval_t * qitv = elina_interval_alloc();
		    elina_interval_t * pitv = elina_interval_alloc();
		    elina_interval_set_double(qitv,-q->inf,q->sup);
		    elina_interval_set_double(pitv,-p->inf,p->sup);
		    elina_interval_div(tmp, qitv, pitv,ELINA_SCALAR_DOUBLE);
		    elina_interval_free(qitv);
		    elina_interval_free(pitv);
		    if (!elina_scalar_infty(tmp->inf) && elina_scalar_equal(tmp->inf,tmp->sup)) elina_interval_neg(array[i]->itv, tmp);
		   // else {
			//bound_neg(tmp->inf, tmp->sup);
			//elina_interval_set_num(array[i]->itv, bound_numref(tmp->inf));
		    //}
		    i++;
		}
		p = p->n ;
		q = q->n ;
	    } else if (p->pnsym->index < q->pnsym->index) {
		
		array[i] = (obj*)calloc(1,sizeof(obj));
		array[i]->coeff = elina_interval_alloc();
		array[i]->itv = elina_interval_alloc();
		if (z->hypercube) {
		    elina_interval_abs_double(tmp, p->inf, p->sup);
		    elina_interval_set(array[i]->coeff,tmp);
		} else {
		    if (zonotope_noise_symbol_cons_get_dimpos(pr, &dim, p->pnsym->index, z)) {
			elina_interval_set(dim_itv, z->gamma[dim]);
			elina_interval_middev(&mid_inf, &mid_sup, &dev_inf, &dev_sup, dim_itv->inf->val.dbl, dim_itv->sup->val.dbl);
			elina_interval_abs_double(tmp, p->inf, p->sup);
			elina_interval_set_double(dev,dev_inf,dev_sup);
			elina_interval_mul(array[i]->coeff, tmp, dev, ELINA_SCALAR_DOUBLE);
		    } else {
			elina_interval_abs_double(tmp, p->inf, p->sup);
			elina_interval_set(array[i]->coeff,tmp);
		    }
		}
		elina_interval_set_double(array[i]->itv,0,0);
		i++;
		p = p->n ;
	    } else {
		q = q->n ;
	    }
	} else if (p) {
	    
	    array[i] = (obj*)calloc(1,sizeof(obj));
	    array[i]->coeff = elina_interval_alloc();
	    array[i]->itv = elina_interval_alloc();
	    if (z->hypercube) {
		elina_interval_abs_double(tmp, p->inf, p->sup);
		elina_interval_set(array[i]->coeff,tmp);
	    } else {
		if (zonotope_noise_symbol_cons_get_dimpos(pr, &dim, p->pnsym->index, z)) {
		    elina_interval_set(dim_itv, z->gamma[dim]);
		    elina_interval_middev(&mid_inf, &mid_sup, &dev_inf, &dev_sup, dim_itv->inf->val.dbl, dim_itv->sup->val.dbl);
		    elina_interval_abs_double(tmp, p->inf, p->sup);
		    elina_interval_set_double(dev,dev_inf,dev_sup);
		    elina_interval_mul(array[i]->coeff, tmp, dev, ELINA_SCALAR_DOUBLE);
		} else {
		    elina_interval_abs_double(tmp, p->inf, p->sup);
		    elina_interval_set(array[i]->coeff,tmp);
		}
	    }
	    elina_interval_set_double(array[i]->itv,0,0);
	    i++;
	    p = p->n ;
	} else {
	    
	    break;
	}
    }
    max = i; 
    
    i = 0;
    obj* obj;
    if (max > 0) {
	while (i<max-1) {
	    if (elina_interval_cmp(array[i]->itv, array[i+1]->itv)) {
		i++;
	    } else {
		obj = array[i];
		array[i] = array[i+1];
		array[i+1] = obj;
		if (i) i--;
	    }
	}
	
	elina_interval_t *sum = elina_interval_alloc();
	elina_interval_t *q1 = elina_interval_alloc();
	elina_interval_t *q2 = elina_interval_alloc();
	
	elina_interval_set(sum, array[0]->coeff);
	for (i=1; i<max; i++) elina_interval_add(sum, sum, array[i]->coeff, ELINA_SCALAR_DOUBLE);
	elina_interval_set(q2, sum);
	i = 0;
	if (max == 1) {
	    elina_interval_set(*res, array[i]->itv);
	} else {
	    do {
		elina_interval_add(q1, q1, array[i]->coeff,ELINA_SCALAR_DOUBLE);
		elina_interval_sub(q2,q2,array[i]->coeff, ELINA_SCALAR_DOUBLE);
		i++;
	    } while (elina_interval_cmp(q1,q2));
	    /* we shall have q1 >= q2 */
	    if (elina_interval_equal(q1,q2)) {
		/* all lambda between [array[i-1].itv, array[i].itv] are solutions */
		
		elina_interval_set(*res, array[i-1]->itv);
	    }
	    else {
		/* array[i-1].itv is the only possible solution */
		elina_interval_set(*res, array[i-1]->itv);
	    }
	}
	elina_interval_free(sum); 
	elina_interval_free(q1); 
	elina_interval_free(q2);
    }
    elina_interval_free(dev);
    elina_interval_free(tmp);
    elina_interval_free(dim_itv);
    for (i=0;i<max;i++) {
	elina_interval_free(array[i]->coeff);
	elina_interval_free(array[i]->itv);
	free(array[i]);
    }
    free(array);
}


static inline void zonotope_noise_symbol_cons_get_gamma(zonotope_internal_t * pr, double *res_inf, double *res_sup, uint_t nsymIndex, zonotope_t* z)
{
    if (z->hypercube) {
	*res_inf = 1.0;
        *res_sup = 1.0;
    } else {
	elina_dim_t dim;
	if (zonotope_noise_symbol_cons_get_dimpos(pr, &dim, nsymIndex, z)) {
		*res_inf = -z->gamma[dim]->inf->val.dbl;
		*res_sup = z->gamma[dim]->sup->val.dbl;
		
	}
	else {
		*res_inf = 1.0;
        	*res_sup = 1.0;
	}
    }
}

zonotope_aff_t* zonotope_aff_add(zonotope_internal_t* pr, zonotope_aff_t* exprA, zonotope_aff_t* exprB, zonotope_t* abs);

static inline void zonotope_set_lincons_dim(zonotope_internal_t* pr, elina_dim_t dim) {
    pr->moo.p[0].linexpr0->p.linterm[0].dim = dim;
    pr->moo.p[1].linexpr0->p.linterm[0].dim = dim;
}

static inline bool zonotope_insert_constrained_noise_symbol(zonotope_internal_t *pr, elina_dim_t* res, uint_t nsymIndex, zonotope_t *z)
{
    uint_t j = 0;
    void* dst = NULL;
    uint_t size = zonotope_noise_symbol_cons_get_dimension(pr, z);
    uint_t dim = 0;
    bool addconsnsym = false;
   
    /* resize nsymcons array if needed */
    if (size>=z->size) {
	z->size += 128;
	z->nsymcons = (elina_dim_t*)realloc(z->nsymcons, (z->size)*sizeof(elina_dim_t));
	z->gamma = (elina_interval_t**)realloc(z->gamma, (z->size)*sizeof(elina_interval_t*));
    }
    if (size == 0) {
	z->nsymcons[size] = nsymIndex;
	z->gamma[size] = pr->ap_muu;
	*res = size;
	pr->dimchange->dim[0] = size;
	elina_abstract0_add_dimensions(pr->manNS, true, z->abs, pr->dimchange, false);
	addconsnsym = true;
    } else if (nsymIndex > z->nsymcons[size-1]) {
	z->nsymcons[size] = nsymIndex;
	z->gamma[size] = pr->ap_muu;
	*res = size;
	pr->dimchange->dim[0] = size;
	elina_abstract0_add_dimensions(pr->manNS, true, z->abs, pr->dimchange, false);
        
	addconsnsym= true;
    } else {
	if (!findKR((uint_t *)&dim, nsymIndex, (uint_t*)z->nsymcons, size)) {
	    dst = (void *)(&z->nsymcons[dim+1]);
	    memmove(dst,(void*)(&z->nsymcons[dim]),(size-dim)*sizeof(elina_dim_t));
	    z->nsymcons[dim] = nsymIndex;
	    dst = (void *)(&z->gamma[dim+1]);
	    memmove(dst,(void*)(&z->gamma[dim]),(size-dim)*sizeof(elina_interval_t*));
	    z->gamma[dim] = pr->ap_muu;
	    *res = dim;
	    pr->dimchange->dim[0] = dim;
	    elina_abstract0_add_dimensions(pr->manNS, true, z->abs, pr->dimchange, false);
	    addconsnsym = true;
	} else {
	    *res = dim;
	    addconsnsym = false;
	}
    }
    
    if (addconsnsym) {
	zonotope_set_lincons_dim(pr, pr->dimchange->dim[0]);
        //elina_abstract0_fprint(stdout,pr->manNS,z->abs,NULL);
	elina_abstract0_meet_lincons_array(pr->manNS, true, z->abs, &pr->moo);
         //elina_abstract0_fprint(stdout,pr->manNS,z->abs,NULL);
    }
    //printf("addconsnsym: %d\n",addconsnsym);
    return addconsnsym;
}

zonotope_aff_t * zonotopeaff_from_linexpr0(zonotope_internal_t* pr, elina_linexpr0_t * expr, zonotope_t *z);
    
static inline elina_linexpr0_t * elina_linexpr0_from_zonotope(zonotope_internal_t* pr, zonotope_aff_t * aff, zonotope_t *z){
	elina_linexpr0_t *res = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,aff->l);
	elina_linexpr0_set_cst_interval_double(res, -aff->c_inf,aff->c_sup);		
	uint_t k = 0;
	elina_dim_t dim;
	zonotope_aaterm_t *p;
	for(p=aff->q; p; p=p->n){
		elina_coeff_init(&res->p.linterm[k].coeff, ELINA_COEFF_SCALAR);
		elina_coeff_set_interval_double(&res->p.linterm[k].coeff, -p->inf,p->sup);
		/* update a->abs with new constrained noise symbols */
		zonotope_insert_constrained_noise_symbol(pr, &dim, p->pnsym->index, z);
		res->p.linterm[k].dim = dim;
        //printf("k: %d dim: %d\n",k,dim);
		k++;
	}
	return res;
}



/* reduce the center and coefficients of the central part (C) to smaller intervals and add a new noise symbol */
static inline bool zonotope_aff_reduce(zonotope_internal_t* pr, zonotope_aff_t *expr)
{
    zonotope_aaterm_t *p,*t;
    double eps = 0.0;
    double err = 0.0;
    //elina_scalar_t *eps = elina_scalar_alloc();
    //elina_scalar_t *err = elina_scalar_alloc();
    //elina_interval_t *mid = elina_interval_alloc();
    //elina_interval_t *dev = elina_interval_alloc();
    double mid_inf = 0.0;
    double mid_sup = 0.0;
    double dev_inf = 0.0;
    double dev_sup = 0.0;
    double sum_inf = 0.0;
    double sum_sup = 0.0;
    double itv_inf = 0.0;
    double itv_sup = 0.0;
    //elina_interval_t *sum = elina_interval_alloc();
    //elina_interval_t *itv = elina_interval_alloc();
    bool ok = false;
    //printf("aff\n");
    //zonotope_aff_fprint(pr,stdout,expr);
    if ((expr->c_inf==INFINITY) || (expr->c_sup==INFINITY)){
        ok = false;
        //printf("ok2\n");
    }
    else {
        eps = 5*1.11022302462515654042e-16;  /* threshold : last bit in the mantissa in double precision */
        err = expr->c_sup + expr->c_inf; 
	if (err> eps) {
	    elina_interval_middev(&mid_inf, &mid_sup,&dev_inf,&dev_sup, expr->c_inf, expr->c_sup);
	    expr->c_inf = mid_inf;
	    expr->c_sup = mid_sup;
	    sum_inf = sum_inf + dev_inf;
	    sum_sup = sum_sup + dev_sup;
	}
	for(p = expr->q; p; p=p->n) {
	    if ((p->inf==INFINITY) || (p->sup==INFINITY)) {
		if ((p->inf==INFINITY) && (p->sup==INFINITY)) {
		    /* reduce to top */
		    if (expr->q) zonotope_aaterm_list_free(pr, expr->q);
		    expr->q = NULL;
		    expr->end = NULL;
		    expr->l = 0;
		    expr->c_inf = INFINITY;
		    expr->c_sup = INFINITY;
		    expr->itv_inf = INFINITY;
		    expr->itv_sup = INFINITY; 
		}
		ok = false;
		sum_inf = 0.0;
		sum_sup = 0.0;
		break;
	    } else {
		double err1 =p->sup + p->inf;
		
		err = err1;
		if ((err> eps)) {
			
		    elina_interval_middev(&mid_inf,&mid_sup, &dev_inf,&dev_sup, p->inf,p->sup);
		    p->inf = mid_inf;
		    p->sup = mid_sup;
		    sum_inf = sum_inf + dev_inf;
	    	    sum_sup = sum_sup + dev_sup;
		}
	    }
	}
        
        //elina_interval_fprint(stdout,sum);
	if ((!sum_inf) && (!sum_sup)) ok = false;
	else {
		
	    zonotope_aff_noise_symbol_create(pr, expr, sum_sup, sum_sup, UN);
	    ok = true; /* adding new symbol */
	}
    }
    return ok;
}

/* z1->itv and z2->itv may be different */
static inline bool zonotope_aff_is_eq(zonotope_internal_t* pr, zonotope_aff_t *z1, zonotope_aff_t *z2)
{
    if (z1 == z2) return true;
    else if (z1->l != z2->l) return false;
    else if ((z1->c_inf!= z2->c_inf) || (z1->c_sup!=z2->c_sup)) return false;
    else {
	zonotope_aaterm_t *p, *q;
	for (p=z1->q, q=z2->q; p && q;) {
	    if (p->pnsym != q->pnsym) return false;
	    else if ((p->inf!=q->inf) || (p->sup!=q->sup)) return false;
	    else {
		p = p->n;
		q = q->n;
	    }
	}
	return true;
    }
}

/* extension of the argmin operator to intervals:
 * - if "a" or "b" contains 0, then the argmin is 0
 * - if "a" is a pos interval and "b" is a neg interval, then the argmin is 0
 * - if "a" is a neg interval and "b" is a pos interval, then the argmin is 0
 * - if "a" and "b" are pos, then choose the one with the minimum inf
 * - if "a" and "b" are neg, then choose the one with the maximum sup
 * RETURN: -1 if argmin = a; 1 if argmin = b; 0 otherwise
 */


static inline int argmin(zonotope_internal_t* pr, elina_interval_t *res, elina_interval_t *a, elina_interval_t *b)
    /* IN: a, b */
    /* OUT: int */
{
    elina_interval_t *zero = elina_interval_alloc();
    int dir = 0;
    if (((elina_scalar_cmp(zero->sup,a->sup)<=0)&&(elina_scalar_cmp(zero->inf,a->inf)<=0)) || ((elina_scalar_cmp(zero->sup,b->sup)<=0)&&(elina_scalar_cmp(zero->inf,b->inf)<=0))) {
	elina_interval_set_double(res,0,0);
	dir = 0;
    } else if (elina_scalar_sgn(a->inf)<=0 && elina_scalar_sgn(b->sup)<=0) {
	elina_interval_set_double(res,0,0);
	dir = 0;
    } else if (elina_scalar_sgn(a->sup)<=0 && elina_scalar_sgn(b->inf)<=0) {
	elina_interval_set_double(res,0,0);
	dir = 0;
    } else {
	/* a and b have the same sign */
	if (elina_scalar_sgn(a->inf)<=0) {
	    if (elina_scalar_cmp(a->inf, b->inf) >= 0) {
		elina_interval_set(res, a);
		dir = -1;
	    } else {
		elina_interval_set(res, b);
		dir = 1;
	    }
	} else {
	    if (elina_scalar_cmp(a->sup, b->sup) >= 0) {
		elina_interval_set(res, a);
		dir = -1;
	    } else {
		elina_interval_set(res, b);
		dir = 1;
	    }
	}
    }
    elina_interval_free(zero);
	
    return dir;
}


/* uses the internal structure pr->dimtoremove */
/* starts with a memory zone setted to zero  (memset) */
static inline void zonotope_delete_constrained_noise_symbol(zonotope_internal_t *pr, uint_t nsymIndex, zonotope_t *z)
{
    uint_t size = zonotope_noise_symbol_cons_get_dimension(pr, z);
    if (size == 0) {
	/* the constrained noise symbol abstract object is already empty */
    } else {
	uint_t dim;
	/* worst case : all constrained nsym have to be removed */
	if (findKR(&dim, nsymIndex, (uint_t*)z->nsymcons, size)) {
	    pr->dimtoremove[dim]++;
	    if (pr->dimtoremove[dim] == z->dims) {
		void* dst = NULL;
		//fprintf(stdout, "remove %d\n", dim);
		pr->dimchange->dim[0] = dim;
		elina_abstract0_remove_dimensions(pr->manNS, true, z->abs, pr->dimchange);
		if (size-1 == 0) z->hypercube = true;
		/* UPDATE z->nsymcons and z->gamma */
		dst = (void *)(&z->nsymcons[dim]);
		memmove(dst,(void*)(&z->nsymcons[dim+1]),(size-dim-1)*sizeof(elina_dim_t));
		//fprintf(stdout, "toto %d #################### \n",size-1-1);
		if (size >= 1){
			 z->nsymcons[size-1] = 0;
		} 
		else {
			z->nsymcons[0] = 0;
		}
		if (z->gamma[dim]) {
		    if (z->gamma[dim] != pr->ap_muu) {elina_interval_free(z->gamma[dim]);z->gamma[dim] = NULL;}
		}
		dst = (void *)(&z->gamma[dim]);
		memmove(dst,(void*)(&z->gamma[dim+1]),(size-dim-1)*sizeof(elina_interval_t*));
	    }
	}
    }
}


static inline void elina_interval_sub_double(elina_interval_t* a, double inf1, double sup1, double inf2, double sup2)
{
  a->inf->val.dbl = inf1 + sup2;
  a->sup->val.dbl = sup1 + inf2;
}

static inline zonotope_aff_t * zonotope_aff_join_constrained6(zonotope_internal_t* pr, zonotope_aff_t *exp1, zonotope_aff_t *exp2, zonotope_t* z1, zonotope_t* z2, zonotope_t* z3)
{
    if(!exp1 || !exp2){
	return NULL;
    }
    zonotope_aff_t * res = zonotope_aff_alloc_init(pr);
    elina_interval_t *tmp = elina_interval_alloc(); 
    elina_interval_t *tmp1 = elina_interval_alloc(); 
    elina_interval_t *tmp2 = elina_interval_alloc();
    elina_interval_t *betaA = elina_interval_alloc(); 
    elina_interval_t *betaB = elina_interval_alloc();
    elina_interval_t *c = elina_interval_alloc();
    elina_interval_t *c1 = elina_interval_alloc();
    elina_interval_t *c2 = elina_interval_alloc();
    elina_interval_t *d = elina_interval_alloc();
    elina_interval_t *d1 = elina_interval_alloc();
    elina_interval_t *d2 = elina_interval_alloc();
    elina_interval_t *mid = elina_interval_alloc();
    elina_interval_t *dev = elina_interval_alloc();
    double mid_inf = 0.0;
    double mid_sup = 0.0;
    double dev_inf = 0.0;
    double dev_sup = 0.0;
    elina_interval_t *nsymItv1 = elina_interval_alloc();
    elina_interval_t *nsymItv2 = elina_interval_alloc();
    elina_interval_t *pmptr = elina_interval_alloc();
    elina_interval_t *qmptr = elina_interval_alloc();
    elina_interval_t *argminpq = elina_interval_alloc();

    zonotope_aaterm_t *p, *q, *ptr;
    res->itv_inf = fmax(exp1->itv_inf,exp2->itv_inf);
    res->itv_sup = fmax(exp1->itv_sup,exp2->itv_sup);
    //elina_scalar_min(res->itv->inf,exp1->itv->inf,exp2->itv->inf);
    //elina_scalar_max(res->itv->sup,exp1->itv->sup,exp2->itv->sup);
    ptr = NULL;
    int s = 0;

    if (exp1->q || exp2->q) {
	elina_interval_set_double(c1, exp1->c_inf,exp1->c_sup);
	elina_interval_set_double(c2, exp2->c_inf, exp2->c_sup);
	ptr = zonotope_aaterm_alloc_init();
	for(p = exp1->q, q = exp2->q; p || q;) {
	    if (p && q) {
		if (p->pnsym->index == q->pnsym->index) {
		    zonotope_noise_symbol_cons_get_gamma(pr, &nsymItv1->inf->val.dbl, &nsymItv1->sup->val.dbl, p->pnsym->index, z1);
		    zonotope_noise_symbol_cons_get_gamma(pr, &nsymItv2->inf->val.dbl, &nsymItv2->sup->val.dbl, p->pnsym->index, z2);
		    if (p->pnsym->type == UN) {
			
			elina_double_interval_mul(&tmp->inf->val.dbl,&tmp->sup->val.dbl , nsymItv1->inf->val.dbl, nsymItv1->sup->val.dbl, p->inf, p->sup);
			elina_interval_add(betaA, betaA, tmp, ELINA_SCALAR_DOUBLE);
			elina_double_interval_mul(&tmp->inf->val.dbl, &tmp->sup->val.dbl, nsymItv2->inf->val.dbl, nsymItv2->sup->val.dbl, q->inf, q->sup);
			elina_interval_add(betaB, betaB, tmp, ELINA_SCALAR_DOUBLE);
		    } else if (p->pnsym->type == IN) {
			
			ptr->pnsym = p->pnsym;
			elina_interval_t *pitv = elina_interval_alloc();
			elina_interval_t *qitv = elina_interval_alloc();
			elina_interval_t *ptritv = elina_interval_alloc();
			elina_interval_set_double(pitv,p->inf,p->sup);
			elina_interval_set_double(qitv,q->inf,q->sup);		
			s = argmin(pr, ptritv, pitv, qitv);
			
			ptr->inf = ptritv->inf->val.dbl;
			ptr->sup = ptritv->sup->val.dbl;
			elina_interval_free(pitv);
			elina_interval_free(qitv);
			elina_interval_free(ptritv);
			if (s == -1) {
			    elina_interval_sub_double(qmptr, q->inf, q->sup, ptr->inf, ptr->sup);
			    //elina_interval_set_int(pmptr, 0,0);
			} else if (s == 1) {
			    //elina_interval_set_int(qmptr, 0,0);
			    elina_interval_sub_double(pmptr, p->inf, p->sup, ptr->inf, ptr->sup);
			} else {
			    elina_interval_set_double(qmptr,q->inf, q->sup);
			    elina_interval_set_double(pmptr,p->inf, p->sup);
			}
		    }
		    p = p->n ;
		    q = q->n ;
		} else if (p->pnsym->index < q->pnsym->index) {
			
		    zonotope_noise_symbol_cons_get_gamma(pr, &nsymItv1->inf->val.dbl, &nsymItv1->sup->val.dbl, p->pnsym->index, z1);
		    zonotope_noise_symbol_cons_get_gamma(pr, &nsymItv2->inf->val.dbl, &nsymItv2->sup->val.dbl, p->pnsym->index, z2);
		    if (p->pnsym->type == UN) {
			elina_double_interval_mul(&tmp->inf->val.dbl, &tmp->sup->val.dbl, nsymItv1->inf->val.dbl, nsymItv1->sup->val.dbl, p->inf, p->sup);
			elina_interval_add(betaA, betaA, tmp, ELINA_SCALAR_DOUBLE);
		    } else if (p->pnsym->type == IN) {
			elina_interval_set_double(pmptr,p->inf, p->sup);
			elina_interval_set_double(qmptr,0,0);
		    }
		    zonotope_delete_constrained_noise_symbol(pr, p->pnsym->index, z3);
		    p = p->n;
		} else {
			
		    zonotope_noise_symbol_cons_get_gamma(pr, &nsymItv1->inf->val.dbl, &nsymItv1->sup->val.dbl, q->pnsym->index, z1);
		    zonotope_noise_symbol_cons_get_gamma(pr, &nsymItv2->inf->val.dbl, &nsymItv2->sup->val.dbl, q->pnsym->index, z2);
		    if (q->pnsym->type == UN) {
			elina_double_interval_mul(&tmp->inf->val.dbl, &tmp->sup->val.dbl, nsymItv2->inf->val.dbl, nsymItv2->sup->val.dbl, q->inf, q->sup);
			elina_interval_add(betaB, betaB, tmp, ELINA_SCALAR_DOUBLE);
		    } else if (q->pnsym->type == IN) {
			elina_interval_set_double(qmptr,q->inf, q->sup);
			elina_interval_set_double(pmptr,0,0);
		    }
		    zonotope_delete_constrained_noise_symbol(pr, q->pnsym->index, z3);
		    q = q->n;
		}
	    } else if (p) {
		
		zonotope_noise_symbol_cons_get_gamma(pr, &nsymItv1->inf->val.dbl, &nsymItv1->sup->val.dbl, p->pnsym->index, z1);
		zonotope_noise_symbol_cons_get_gamma(pr, &nsymItv2->inf->val.dbl, &nsymItv2->sup->val.dbl, p->pnsym->index, z2);
		if (p->pnsym->type == UN) {
		    elina_double_interval_mul(&tmp->inf->val.dbl, &tmp->sup->val.dbl, nsymItv1->inf->val.dbl, nsymItv1->sup->val.dbl, p->inf, p->sup);
		    elina_interval_add(betaA, betaA, tmp, ELINA_SCALAR_DOUBLE);
			//printf("nysmitv %g %g\n",p->inf,p->sup);
			//elina_interval_fprint(stdout,nsymItv1);
			//elina_interval_fprint(stdout,tmp);
			//elina_interval_fprint(stdout,betaA);
			//fflush(stdout);
		} else if (p->pnsym->type == IN) {
		    elina_interval_set_double(pmptr,p->inf, p->sup);
		    elina_interval_set_double(qmptr,0,0);
		}
		zonotope_delete_constrained_noise_symbol(pr, p->pnsym->index, z3);
		p = p->n;
	    } else {
		
		zonotope_noise_symbol_cons_get_gamma(pr, &nsymItv1->inf->val.dbl, &nsymItv1->sup->val.dbl, q->pnsym->index, z1);
		zonotope_noise_symbol_cons_get_gamma(pr, &nsymItv2->inf->val.dbl, &nsymItv2->sup->val.dbl, q->pnsym->index, z2);
		if (q->pnsym->type == UN) {
		    elina_double_interval_mul(&tmp->inf->val.dbl, &tmp->sup->val.dbl, nsymItv2->inf->val.dbl, nsymItv2->sup->val.dbl, q->inf, q->sup);
		    elina_interval_add(betaB, betaB, tmp, ELINA_SCALAR_DOUBLE);
		} else if (q->pnsym->type == IN) {
		    elina_interval_set_double(qmptr,q->inf, q->sup);
		    elina_interval_set_double(pmptr,0,0);
		}
		zonotope_delete_constrained_noise_symbol(pr, q->pnsym->index, z3);
		q = q->n;
	    }
	    elina_double_interval_mul(&tmp1->inf->val.dbl, &tmp1->sup->val.dbl, nsymItv1->inf->val.dbl, nsymItv1->sup->val.dbl, pmptr->inf->val.dbl, pmptr->sup->val.dbl);
	    elina_interval_add(c1, c1, tmp1, ELINA_SCALAR_DOUBLE);
	    elina_double_interval_mul(&tmp2->inf->val.dbl, &tmp2->sup->val.dbl, nsymItv2->inf->val.dbl,  nsymItv2->sup->val.dbl, qmptr->inf->val.dbl, qmptr->sup->val.dbl);
	    elina_interval_add(c2, c2, tmp2, ELINA_SCALAR_DOUBLE);
	    elina_interval_set_double(pmptr,0,0);
	    elina_interval_set_double(qmptr,0,0);
	    if ((ptr->inf==0) && (ptr->sup==0)) {
		if (!(p||q)) {
		    /* the last iteration */
		    zonotope_aaterm_free(pr, ptr);
		    if (res->end) res->end->n = NULL;
		}
	    } else {
		/* keep this term */
		if (!res->q) res->q = ptr;
		res->end = ptr;
		res->l++;
		if (p||q) {
		    /* continuing */
		    ptr->n = zonotope_aaterm_alloc_init();
		    ptr=ptr->n;
		} else {
		    /* the last iteration */
		}
	    }
	}
	//printf("betaA\n");
	//elina_interval_print(betaA);
	//printf("betaB\n");
	//elina_interval_print(betaB);
	//fflush(stdout);
	elina_interval_middev(&mid->inf->val.dbl, &mid->sup->val.dbl, &dev->inf->val.dbl,&dev->sup->val.dbl, betaA->inf->val.dbl, betaA->sup->val.dbl);
	elina_interval_set(d1, dev);
	elina_interval_add(c1, c1, mid, ELINA_SCALAR_DOUBLE);

	elina_interval_middev(&mid->inf->val.dbl, &mid->sup->val.dbl, &dev->inf->val.dbl,&dev->sup->val.dbl, betaB->inf->val.dbl,betaB->sup->val.dbl);
	elina_interval_set(d2, dev);
	elina_interval_add(c2, c2, mid, ELINA_SCALAR_DOUBLE);

	
	elina_scalar_t *c0 = elina_scalar_alloc(); 
	elina_scalar_t *min = elina_scalar_alloc(); 
	elina_scalar_t *max = elina_scalar_alloc(); 
	elina_scalar_t *dmin = elina_scalar_alloc(); 
	elina_scalar_t *dmax = elina_scalar_alloc();
	elina_interval_t *beta = elina_interval_alloc();

	elina_scalar_t **array = (elina_scalar_t **)malloc(4*sizeof(elina_scalar_t *));
	array[0] = elina_scalar_alloc();
	array[1] = elina_scalar_alloc();
	array[2] = elina_scalar_alloc();
	array[3] = elina_scalar_alloc();

	elina_scalar_set(array[0],c1->inf);
	elina_scalar_neg(array[0],array[0]);
	elina_scalar_add(array[0],array[0],d1->sup, ELINA_SCALAR_DOUBLE);
	elina_scalar_set(array[1],c2->inf);
	elina_scalar_neg(array[1],array[1]);
	elina_scalar_add(array[1],array[1],d2->sup, ELINA_SCALAR_DOUBLE);
	elina_scalar_set(array[2],c1->sup);
	elina_scalar_add(array[2],array[2],d1->sup, ELINA_SCALAR_DOUBLE);
	elina_scalar_set(array[3],c2->sup);
	elina_scalar_add(array[3],array[3],d2->sup, ELINA_SCALAR_DOUBLE);

	int i = 0;
	elina_scalar_set(dmax, array[0]);
	for (i=1;i<4;i++)
	    if (elina_scalar_cmp(dmax,array[i]) < 0) elina_scalar_set(dmax, array[i]);

	elina_scalar_set(array[0],c1->inf);
    	//elina_scalar_neg(array[0],array[0]);
	elina_scalar_add(array[0],array[0],d1->sup, ELINA_SCALAR_DOUBLE);
	elina_scalar_set(array[1],c2->inf);
    	//elina_scalar_neg(array[1],array[1]);
	elina_scalar_add(array[1],array[1],d2->sup, ELINA_SCALAR_DOUBLE);
	elina_scalar_set(array[2],c1->sup);
	elina_scalar_neg(array[2],array[2]);
	elina_scalar_add(array[2],array[2],d1->sup, ELINA_SCALAR_DOUBLE);
	elina_scalar_set(array[3],c2->sup);
	elina_scalar_neg(array[3],array[3]);
	elina_scalar_add(array[3],array[3],d2->sup, ELINA_SCALAR_DOUBLE);

	elina_scalar_set(dmin, array[0]);
	for (i=1;i<4;i++)
	    if (elina_scalar_cmp(dmin,array[i]) < 0) elina_scalar_set(dmin, array[i]);

	elina_scalar_add(c0,dmax,dmin, ELINA_SCALAR_DOUBLE);
	elina_scalar_div_2(c0, c0);
	
	elina_interval_set_double(beta,-c0->val.dbl,c0->val.dbl);

	elina_scalar_sub(c0,dmax,dmin, ELINA_SCALAR_DOUBLE);
	elina_scalar_div_2(c0, c0);
	res->c_inf = -c0->val.dbl;
	res->c_sup = c0->val.dbl; 

	if (elina_scalar_cmp(c1->inf,c2->inf) < 0) {
	    elina_scalar_set(min, c2->inf);
	} else {
	    elina_scalar_set(min, c1->inf);
	}
	elina_scalar_neg(min,min);

	if (elina_scalar_cmp(c1->sup,c2->sup) < 0) {
	    elina_scalar_set(max, c2->sup);
	} else {
	    elina_scalar_set(max, c1->sup);
	}
	//printf("beta\n");
	//zonotope_aff_fprint(pr,stdout,exp1); 
	//zonotope_aff_fprint(pr,stdout,exp2);
	//zonotope_aff_fprint(pr,stdout,res);
	//elina_interval_fprint(stdout,betaA);
	//elina_interval_fprint(stdout,betaB);
	//elina_interval_fprint(stdout,beta);
	//printf("\n");
	//elina_scalar_fprint(stdout,c0);
	//printf("\n");
	//elina_scalar_fprint(stdout,min);
	//printf("\n");
	//elina_scalar_fprint(stdout,max);
	//printf("\n");
	//fflush(stdout);
	if (elina_scalar_cmp(c0,min) <= 0 || elina_scalar_cmp(c0,max) >= 0) {
	    zonotope_aff_free(pr, res);
	    elina_interval_sub_double(tmp, d2->inf->val.dbl, d2->sup->val.dbl, d1->inf->val.dbl, d1->sup->val.dbl);
	    //elina_interval_sub(tmp,d2,d1, ELINA_SCALAR_DOUBLE);
	    if (elina_scalar_sgn(tmp->inf)<=0) res = exp2;
	    else res = exp1;
	} else{
		
		 zonotope_aff_noise_symbol_create(pr, res, beta->inf->val.dbl, beta->sup->val.dbl, UN);
	}
	elina_scalar_free(array[0]);
	elina_scalar_free(array[1]);
	elina_scalar_free(array[2]);
	elina_scalar_free(array[3]);
	free(array);
	elina_interval_free(beta);
	elina_scalar_free(c0); 
	elina_scalar_free(min); 
	elina_scalar_free(max); 
	elina_scalar_free(dmin); 
	elina_scalar_free(dmax);
	
    } else {
	zonotope_aff_add_itv(pr, res, res->itv_inf, res->itv_sup, UN);
    }
    elina_interval_free(tmp); elina_interval_free(tmp1); elina_interval_free(tmp2);
    elina_interval_free(betaA); elina_interval_free(betaB);
    elina_interval_free(c);
    elina_interval_free(c1);
    elina_interval_free(c2);
    elina_interval_free(d);
    elina_interval_free(d1);
    elina_interval_free(d2);
    elina_interval_free(nsymItv1);
    elina_interval_free(nsymItv2);
    elina_interval_free(argminpq);
    elina_interval_free(pmptr);
    elina_interval_free(qmptr);
    elina_interval_free(mid);
    elina_interval_free(dev);
  
    return res;
}

static inline void zonotope_noise_symbol_cons_set_hypercube(zonotope_internal_t* pr, zonotope_t* z)
{
    uint_t size = zonotope_noise_symbol_cons_get_dimension(pr, z);
    elina_interval_t** tinterval = (elina_interval_t**)malloc(size*sizeof(elina_interval_t*));
    uint_t i = 0;
    for (i=0; i<size; i++) tinterval[i] = pr->ap_muu;

    elina_abstract0_free(pr->manNS, pr->nsymhypercube);
    pr->nsymhypercube = elina_abstract0_of_box(pr->manNS, 0, size, tinterval);
    free(tinterval);
}




static inline void zonotope_update_noise_symbol_cons_gamma(zonotope_internal_t* pr, zonotope_t* z)
{
    uint_t i = 0;
    uint_t nsymcons_size = zonotope_noise_symbol_cons_get_dimension(pr, z);
    elina_interval_t* bound = NULL;
    zonotope_noise_symbol_cons_set_hypercube(pr, z);
    if (elina_abstract0_is_leq(pr->manNS, pr->nsymhypercube, z->abs)) {
        //printf("hypercube %d\n",nsymcons_size);
        //fflush(stdout);
	z->hypercube = true;
	elina_abstract0_free(pr->manNS, z->abs);
	z->abs = elina_abstract0_top(pr->manNS, 0,0);
	for (i=0; i<nsymcons_size; i++) {
	    if (z->gamma[i] != pr->ap_muu && z->gamma[i]) elina_interval_free(z->gamma[i]);
	    z->gamma[i] = NULL;
	}
	memset((void*)z->nsymcons, 0, nsymcons_size*sizeof(elina_dim_t));
    } else {
        //printf("not heypercube %d\n",nsymcons_size);
        //fflush(stdout);
	z->hypercube = false;
	for (i=0; i<nsymcons_size; i++) {
		
	    bound = elina_abstract0_bound_dimension(pr->manNS, z->abs, i);
		
	    if (z->gamma[i] == NULL) {
		z->gamma[i] = bound;
	    } else if (elina_interval_is_leq(z->gamma[i], bound)) {
		elina_interval_free(bound);
	    } else {
		if (z->gamma[i] != pr->ap_muu) elina_interval_free(z->gamma[i]);
		z->gamma[i] = bound;
	    }
	    bound = NULL;
	}
    }
}

static inline void zonotope_internal_free(zonotope_internal_t* pr)
{
    uint_t i = 0;
    if (pr) {
	//elina_internal_free(pr->itv);
	zonotope_aff_free(pr, pr->top);
	zonotope_aff_free(pr, pr->bot);
	
	for (i=0;i<pr->dim;i++) free(pr->epsilon[i]);
	pr->dim = (uint_t)0;
	free(pr->epsilon);
	pr->epsilon= NULL;
	pr->funid = ELINA_FUNID_UNKNOWN;
	pr->man = NULL;
	elina_abstract0_free(pr->manNS, pr->nsymhypercube);
	pr->nsymhypercube = NULL;
	elina_manager_free(pr->manNS);
	pr->manNS = NULL;
	elina_manager_free(pr->box);
	pr->box = NULL;
	elina_interval_free(pr->muu);
	elina_interval_free(pr->ap_muu);
	pr->ap_muu = NULL;
	elina_lincons0_array_clear(&(pr->moo));
	free(pr->dimtoremove);
	elina_dimchange_free(pr->dimchange);
	pr->dimchange = NULL;
	pr->it = 0;
	free(pr->inputns);
	free(pr);
    }
}


static inline zonotope_internal_t* zonotope_internal_alloc(elina_manager_t* manNS)
{
    //CALL();
    zonotope_internal_t* pr = (zonotope_internal_t*)malloc(sizeof(zonotope_internal_t));
    //pr->itv = elina_internal_alloc();
    pr->dim = 0;
    pr->funid = ELINA_FUNID_UNKNOWN;
    pr->man = NULL;
    pr->manNS = manNS;
    pr->box = elina_box_manager_alloc();
    pr->muu = elina_interval_alloc();
    elina_scalar_set_double(pr->muu->inf, (double)1.0);
    elina_scalar_set_double(pr->muu->sup, (double)1.0);
    pr->ap_muu = elina_interval_alloc();
    elina_scalar_set_double(pr->ap_muu->inf, (double)-1.0);
    elina_scalar_set_double(pr->ap_muu->sup, (double)1.0);
    //ap_interval_set_itv(pr->itv, pr->ap_muu, pr->muu);
    pr->moo = elina_lincons0_array_make(2);
    pr->epsilon = (zonotope_noise_symbol_t**)calloc(1024, sizeof(zonotope_noise_symbol_t*));	/* starts with a limit of 1024 noise symbols */
    pr->top = zonotope_aff_top(pr);
    pr->bot = zonotope_aff_bottom(pr);
    pr->nsymhypercube = elina_abstract0_top(pr->manNS, 0,0);
    pr->dimtoremove = (elina_dim_t*)calloc(128, sizeof(elina_dim_t));
    pr->dimchange = elina_dimchange_alloc(0,1);
    /* -eps + 1 */
    elina_linexpr0_t* mnspone = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE, (uint_t)1);	/* mnspone : (m)inus (n)oise (s)ymbol (p)lus (o)ne */
    elina_linexpr0_set_cst_scalar_double(mnspone, 1.0);
    elina_linexpr0_set_coeff_scalar_double(mnspone, (elina_dim_t)0, (double)-1.0);
    /* eps + 1 */
    elina_linexpr0_t* nspone = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE, (uint_t)1);
    elina_linexpr0_set_cst_scalar_double(nspone, 1.0);
    elina_linexpr0_set_coeff_scalar_double(nspone, (elina_dim_t)0, (double)1.0);

    /* 0 <= -eps + 1 */
    pr->moo.p[0] = elina_lincons0_make(ELINA_CONS_SUPEQ, mnspone, NULL);
    /* 0 <= eps + 1 */
    pr->moo.p[1] = elina_lincons0_make(ELINA_CONS_SUPEQ, nspone, NULL);

    pr->inputns = (uint_t*)calloc(1024, sizeof(uint_t));	/* starts with a limit of 1024 noise symbols */
    pr->epssize = 0;
    pr->it = 0;
    return pr;
}

#ifdef __cplusplus
}
#endif

#endif
