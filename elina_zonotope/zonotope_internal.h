/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2017 Department of Computer Science, ETH Zurich
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




#ifdef __cplusplus
extern "C" {
#endif

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
    elina_interval_t*	coeff;	/* coeff, encoded as interval */
} zonotope_aaterm_t;

/************************/
/* Zonotope affine form */
/************************/
struct _zonotope_aff_t {
    elina_interval_t *	c;	/* center */
    zonotope_aaterm_t*	q;	/* first center term (epsilons) aaterm */
    zonotope_aaterm_t*	end;	/* quick jump to the last center term : to add a new term for instance */
    unsigned long long int		l;	/* number of noise symbols */
    unsigned long long int		pby;	/* # pointers to this affine form */
    elina_interval_t*	itv;	/* best known interval concretisation */
};
typedef struct _zonotope_aff_t zonotope_aff_t;

    static inline zonotope_aff_t * zonotope_aff_copy(zonotope_aff_t *src){
        zonotope_aff_t *res = (zonotope_aff_t *)malloc(sizeof(zonotope_aff_t));
        res->c  = elina_interval_alloc_set(src->c);
        res->q = src->q==NULL? NULL: (zonotope_aaterm_t *)malloc(sizeof(zonotope_aaterm_t *));
        res->end = src->end==NULL? NULL:(zonotope_aaterm_t *)malloc(sizeof(zonotope_aaterm_t *));
        res->l = src->l;
        res->itv = elina_interval_alloc_set(src->itv);
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
    elina_interval_t**	box;		/* reduced product with boxes */
    uint_t		intdim;         /* nb of integer variables */
    uint_t		dims;           /* intdim + realdim */
    elina_abstract0_t* 	abs;        	/* nsym abstract object (=contraints over noise symbols)*/
    elina_dim_t*	nsymcons;       /* array of index of constrained noise symbols */
    elina_interval_t**	gamma;		/* pointer to an array which contains the concretisations of constrained noise symbols if any */
    unsigned long long int		size;		/* size of nsymcons and gamma */
    bool		hypercube;	/* true if no constrained nsym */
    elina_interval_t**	g;	/* array of the generators of the zonotope - a oublier */
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
    res->coeff = elina_interval_alloc();
    return res;
}

static zonotope_aff_t* zonotope_aff_alloc_init(zonotope_internal_t *pr)
{
    zonotope_aff_t* a = (zonotope_aff_t*)malloc(sizeof(zonotope_aff_t));
    a->c = elina_interval_alloc();
    a->q = NULL;
    a->end = NULL;
    a->l = 0;
    a->pby = 0;
    a->itv = elina_interval_alloc();
    return a;
}

static inline zonotope_aff_t * zonotope_aff_top(zonotope_internal_t* pr)
{
    zonotope_aff_t* res = zonotope_aff_alloc_init(pr);
    elina_interval_set_top(res->c);
    elina_interval_set_top(res->itv);
    return res;
}
static inline zonotope_aff_t * zonotope_aff_bottom(zonotope_internal_t* pr)
{
    zonotope_aff_t* res = zonotope_aff_alloc_init(pr);
    elina_interval_set_bottom(res->c);
    elina_interval_set_bottom(res->itv);
    return res;
}

static inline zonotope_noise_symbol_t* zonotope_noise_symbol_add(zonotope_internal_t *pr, noise_symbol_t type)
    /* increment the global index of used noise symbols and add the noise symbol in pr->eps */
{
    uint_t dim = pr->dim;
    zonotope_noise_symbol_t* res;
    /* resize epsilon array */
    if ((dim+1) % 1024 == 0) pr->epsilon = (zonotope_noise_symbol_t**)realloc(pr->epsilon, (dim+1024)*sizeof(zonotope_noise_symbol_t*));
    res = pr->epsilon[dim] = (zonotope_noise_symbol_t*)malloc(sizeof(zonotope_noise_symbol_t));
    if (type == IN) {pr->inputns[pr->epssize] = dim; pr->epssize++;}
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
        if (!elina_scalar_infty(ptr->coeff->inf) && elina_scalar_equal(ptr->coeff->inf,ptr->coeff->sup)){
            elina_scalar_fprint(stream,ptr->coeff->sup);
        }
        else elina_interval_fprint(stream, ptr->coeff);
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
        if (!elina_scalar_infty(expr->c->inf) && elina_scalar_equal(expr->c->inf,expr->c->sup)){
            elina_scalar_fprint(stream,expr->c->sup);
        }
        else elina_interval_fprint(stream, expr->c);
        /* Print values */
        for (p=expr->q; p; p=p->n) {
            fprintf(stream," + ");
            zonotope_aaterm_fprint(pr, stream, p);
        }
        fprintf(stream,"\t;");
        elina_interval_fprint(stream, expr->itv);
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


/* Free memory used by one aaterm */
static inline void zonotope_aaterm_free(zonotope_internal_t* pr, zonotope_aaterm_t* term)
{
   // arg_assert(term, abort(););
    term->n = NULL;
    term->pnsym = NULL;
    elina_interval_free(term->coeff);
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
	elina_interval_free(a->c);
	if (a->q) zonotope_aaterm_list_free(pr, a->q);
	a->q = NULL;
	a->end = NULL;
	a->l = (uint_t)0;
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

static inline void zonotope_aff_noise_symbol_create(zonotope_internal_t *pr, zonotope_aff_t *expr, elina_interval_t *coeff, noise_symbol_t type)
{
    elina_interval_t *zero = elina_interval_alloc();
    if (elina_interval_cmp(coeff,zero)>=0) {
	zonotope_aaterm_t* ptr = zonotope_aaterm_alloc_init();
	elina_interval_set(ptr->coeff, coeff);
	ptr->pnsym = zonotope_noise_symbol_add(pr, type);
	if (expr->end) expr->end->n = ptr;
	else expr->q = ptr;
	expr->end = ptr;
	expr->l++;
    }
    elina_interval_free(zero);
}


static inline void elina_interval_middev(elina_interval_t *mid, elina_interval_t *dev, elina_interval_t *a)
{
    elina_scalar_t ** tmp = (elina_scalar_t **)malloc(4*sizeof(elina_scalar_t*));
    for(int i=0; i < 4; i++){
	tmp[i] = elina_scalar_alloc();
    }
    
    if(!elina_scalar_infty(a->inf) && elina_scalar_equal(a->inf,a->sup)){
	elina_interval_set(mid,a);
	elina_interval_set_int(dev,0,0);
    } else if (elina_scalar_infty(a->sup) || elina_scalar_infty(a->inf) || elina_interval_is_bottom(a)) {
	elina_interval_set_top(mid);
	elina_interval_set_top(dev);
    } else {
	/* a = [x,y], x < y,
	 * tmp[0] = x+y */
	/* Rounding to Nearest, only x86 Unix like */
	elina_scalar_add(tmp[0], a->sup, a->inf,ELINA_SCALAR_DOUBLE);
	if (!fesetround(FE_TONEAREST)) {
	    /* tmp[1] = (x+y)/2 -- using ldexp if double */
	    elina_scalar_div_2(tmp[1], tmp[0]);
	    /* mid = [tmp[1], tmp[1]] */
	    elina_scalar_set(mid->sup, tmp[1]);
	    elina_scalar_set(mid->inf, tmp[1]);
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
	    elina_scalar_sub(tmp[0], tmp[1], a->inf,ELINA_SCALAR_DOUBLE);
	    /* tmp[2] = y - (x+y)/2 */
	    elina_scalar_sub(tmp[2], a->sup, tmp[1],ELINA_SCALAR_DOUBLE);
	    /* tmp[3] = max(tmp[0],tmp[2]) -- fmax if double */
	    elina_scalar_max(tmp[3], tmp[0], tmp[2]);
	    /* dev = [tmp[3], tmp[3]] */
	    elina_scalar_set(dev->sup, tmp[3]);
	    elina_scalar_set(dev->inf, tmp[3]);
	}
    }
    for(int i=0; i < 4; i++){
	elina_scalar_free(tmp[i]);
    }
    free(tmp);
}

/* convert an itv to an affine form with a fresh noise symbol then add this form to the affine form expr */
static inline void zonotope_aff_add_itv(zonotope_internal_t* pr, zonotope_aff_t *expr, elina_interval_t* itv, noise_symbol_t type)
{
    /* itv is a non point interval with finite bounds */
    elina_interval_t *mid, *dev;
    mid = elina_interval_alloc();
    dev = elina_interval_alloc();
    if(!elina_scalar_equal(itv->inf,itv->sup)){
	elina_interval_middev(mid, dev, itv);
	//printf("mid\n");
	//elina_scalar_print(mid);
	//printf("\n");
	elina_interval_add(expr->c, expr->c, mid,ELINA_SCALAR_DOUBLE);
	zonotope_aff_noise_symbol_create(pr, expr, dev, type);
    } else elina_interval_add(expr->c, expr->c, itv, ELINA_SCALAR_DOUBLE);
    elina_interval_free(mid); 
    elina_interval_free(dev);
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

static inline void zonotope_aff_bound(zonotope_internal_t* pr, elina_interval_t *res, zonotope_aff_t *expr, zonotope_t* z)
{
    if (elina_interval_is_top(expr->c)) {
	elina_interval_set_top(res);
	return;
    } else {
	elina_dim_t dim;
	elina_interval_t *tmp = elina_interval_alloc();
	elina_interval_t *eps_itv = elina_interval_alloc();
	zonotope_aaterm_t* p;
	elina_interval_set(res,expr->c);
	if (z->hypercube) {
	    for (p=expr->q; p; p=p->n) {
		elina_interval_mul(tmp,p->coeff, pr->muu,ELINA_SCALAR_DOUBLE);
		elina_interval_add(res, res, tmp,ELINA_SCALAR_DOUBLE);
	    }
	} else {
	    elina_linexpr0_t* linexpr0 = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE, 0);
	    elina_linexpr0_set_cst_scalar_int(linexpr0, (int)0);
	    
	    linexpr0->p.linterm = (elina_linterm_t*)malloc(expr->l*sizeof(elina_linterm_t));
	    uint_t k = 0;
	    elina_dim_t dim = 0;
	    for (p=expr->q; p; p=p->n) {
		if (zonotope_noise_symbol_cons_get_dimpos(pr, &dim, p->pnsym->index, z)) {
		    elina_coeff_init(&linexpr0->p.linterm[k].coeff, ELINA_COEFF_INTERVAL);
		    elina_coeff_set_interval(&linexpr0->p.linterm[k].coeff, p->coeff);
		    linexpr0->p.linterm[k].dim = dim;
		    k++;
		} else {
		    elina_interval_mul(tmp, p->coeff, pr->muu,ELINA_SCALAR_DOUBLE);
		    elina_interval_add(res, res, tmp,ELINA_SCALAR_DOUBLE);
		}
	    }
	    linexpr0->size = k;
	    elina_interval_t* elina_itv = elina_abstract0_bound_linexpr(pr->manNS, z->abs, linexpr0);
	    elina_interval_set(tmp, elina_itv);
	    elina_interval_add(res, res, tmp,ELINA_SCALAR_DOUBLE);
	    linexpr0->p.linterm = (elina_linterm_t*)realloc((void* )linexpr0->p.linterm, k*sizeof(elina_linterm_t));
	    elina_linexpr0_free(linexpr0);
	    elina_interval_free(elina_itv);
	}
	elina_interval_free(tmp); 
	elina_interval_free (eps_itv);
    }
}

static inline bool zonotope_aff_is_known_to_be_zero(zonotope_internal_t *pr, zonotope_aff_t *a)
{
    if (!elina_scalar_sgn(a->itv->inf)&& !elina_scalar_sgn(a->itv->sup)) return true;
    else return false;
}

static inline bool zonotope_aff_is_top(zonotope_internal_t* pr, zonotope_aff_t *a)
{
    if (a == pr->top) return true;
    else if (!elina_interval_is_top(a->c)) return false;
    else if (!elina_interval_is_top(a->itv)) return false;
    else if (a->q != NULL) return false;
    else return true;
}
static inline bool zonotope_aff_is_bottom(zonotope_internal_t* pr, zonotope_aff_t *a)
{
    if (a == pr->bot) return true;
    else if (!elina_interval_is_bottom(a->c)) return false;
    else if (!elina_interval_is_bottom(a->itv)) return false;
    else if (a->q != NULL) return false;
    else return true;
}



static inline zonotope_aff_t* zonotope_aff_mul_itv(zonotope_internal_t* pr, zonotope_aff_t* src, elina_interval_t *lambda)
{
    if ((!elina_scalar_sgn(lambda->inf) && !elina_scalar_sgn(lambda->sup) )|| zonotope_aff_is_known_to_be_zero(pr, src)) {
	return zonotope_aff_alloc_init(pr);
    } else if (zonotope_aff_is_bottom(pr, src) || elina_interval_is_bottom(lambda)) {
	return zonotope_aff_bottom(pr);
    } else if (zonotope_aff_is_top(pr, src) || elina_interval_is_top(lambda)) {
	return zonotope_aff_top(pr);
    } else {
	zonotope_aff_t* dst = NULL;
	zonotope_aaterm_t *p,*q;
	q = NULL;
	dst = zonotope_aff_alloc_init(pr);
	elina_interval_mul(dst->c, lambda, src->c,ELINA_SCALAR_DOUBLE);
	if (src->q) {
	    dst->q = q = zonotope_aaterm_alloc_init();
	    for (p=src->q; p; p=p->n) {
		elina_interval_mul(q->coeff, lambda, p->coeff,ELINA_SCALAR_DOUBLE);
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
	elina_interval_mul(dst->itv, src->itv, lambda,ELINA_SCALAR_DOUBLE);
	return dst;
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
    elina_interval_t *mid = elina_interval_alloc();
    elina_interval_t *dev = elina_interval_alloc();
    elina_interval_t *dim_itv = elina_interval_alloc();
    elina_dim_t dim;
    
    for (p=cons->q, q=x->q; p || q;) {
	if (p && q) {
	    if (p->pnsym->index == q->pnsym->index) {
		
		if(elina_scalar_sgn(p->coeff->inf) || elina_scalar_sgn(p->coeff->sup)) {
		    array[i] = (obj*)calloc(1,sizeof(obj));
		    array[i]->coeff = elina_interval_alloc();
		    array[i]->itv = elina_interval_alloc();
		    if (z->hypercube) {
			elina_interval_abs(tmp, p->coeff,ELINA_SCALAR_DOUBLE);
			elina_interval_set(array[i]->coeff,tmp);
		    } else {
			if (zonotope_noise_symbol_cons_get_dimpos(pr, &dim, p->pnsym->index, z)) {
			    elina_interval_set(dim_itv, z->gamma[dim]);
			    elina_interval_middev(mid, dev, dim_itv);
			    elina_interval_abs(tmp, p->coeff, ELINA_SCALAR_DOUBLE);
			    elina_interval_mul(array[i]->coeff, tmp, dev,ELINA_SCALAR_DOUBLE);
			} else {
			    elina_interval_abs(tmp, p->coeff, ELINA_SCALAR_DOUBLE);
			    elina_interval_set(array[i]->coeff,tmp);
			}
		    }
		    elina_interval_div(tmp, q->coeff, p->coeff,ELINA_SCALAR_DOUBLE);
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
		    elina_interval_abs(tmp, p->coeff, ELINA_SCALAR_DOUBLE);
		    elina_interval_set(array[i]->coeff,tmp);
		} else {
		    if (zonotope_noise_symbol_cons_get_dimpos(pr, &dim, p->pnsym->index, z)) {
			elina_interval_set(dim_itv, z->gamma[dim]);
			elina_interval_middev(mid, dev, dim_itv);
			elina_interval_abs(tmp, p->coeff, ELINA_SCALAR_DOUBLE);
			elina_interval_mul(array[i]->coeff, tmp, dev, ELINA_SCALAR_DOUBLE);
		    } else {
			elina_interval_abs(tmp, p->coeff, ELINA_SCALAR_DOUBLE);
			elina_interval_set(array[i]->coeff,tmp);
		    }
		}
		elina_interval_set_int(array[i]->itv,0,0);
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
		elina_interval_abs(tmp, p->coeff, ELINA_SCALAR_DOUBLE);
		elina_interval_set(array[i]->coeff,tmp);
	    } else {
		if (zonotope_noise_symbol_cons_get_dimpos(pr, &dim, p->pnsym->index, z)) {
		    elina_interval_set(dim_itv, z->gamma[dim]);
		    elina_interval_middev(mid, dev, dim_itv);
		    elina_interval_abs(tmp, p->coeff, ELINA_SCALAR_DOUBLE);
		    elina_interval_mul(array[i]->coeff, tmp, dev, ELINA_SCALAR_DOUBLE);
		} else {
		    elina_interval_abs(tmp, p->coeff, ELINA_SCALAR_DOUBLE);
		    elina_interval_set(array[i]->coeff,tmp);
		}
	    }
	    elina_interval_set_int(array[i]->itv,0,0);
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
    elina_interval_free(mid);
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


static inline void zonotope_noise_symbol_cons_get_gamma(zonotope_internal_t * pr, elina_interval_t *res, uint_t nsymIndex, zonotope_t* z)
{
    if (z->hypercube) {
	elina_interval_set(res, pr->muu);
    } else {
	elina_dim_t dim;
	if (zonotope_noise_symbol_cons_get_dimpos(pr, &dim, nsymIndex, z)) {
		
		elina_interval_set(res, z->gamma[dim]);
		
	}
	else elina_interval_set(res, pr->muu);
    }
}

static inline zonotope_aff_t* zonotope_aff_add(zonotope_internal_t* pr, zonotope_aff_t* exprA, zonotope_aff_t* exprB, zonotope_t* abs)
{
    elina_interval_t *box = elina_interval_alloc();
    elina_interval_t *tmp = elina_interval_alloc();
    
    zonotope_aff_t* res = zonotope_aff_alloc_init(pr);
    zonotope_aaterm_t *p, *q, *ptr;
    elina_interval_add(res->c, exprA->c, exprB->c, ELINA_SCALAR_DOUBLE);
    elina_interval_set(box, res->c);
    if (exprA->q || exprB->q) {
	ptr = zonotope_aaterm_alloc_init();
	for(p = exprA->q, q = exprB->q; p || q;) {
	    if (p && q) {
		if (p->pnsym->index == q->pnsym->index) {
		    elina_interval_add(ptr->coeff, p->coeff, q->coeff, ELINA_SCALAR_DOUBLE);
		    ptr->pnsym = p->pnsym;
		    p = p->n ;
		    q = q->n ;
		} else if (p->pnsym->index < q->pnsym->index) {
		    elina_interval_set(ptr->coeff, p->coeff);
		    ptr->pnsym = p->pnsym;
		    p = p->n ;
		} else {
		    elina_interval_set(ptr->coeff, q->coeff);
		    ptr->pnsym = q->pnsym;
		    q = q->n ;
		}
	    } else if (p) {
		elina_interval_set(ptr->coeff, p->coeff);
		ptr->pnsym = p->pnsym;
		p = p->n ;
	    } else {
		elina_interval_set(ptr->coeff, q->coeff);
		ptr->pnsym = q->pnsym;
		q = q->n ;
	    }
	    if (!elina_scalar_sgn(ptr->coeff->inf) && !elina_scalar_sgn(ptr->coeff->sup)) {
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
		zonotope_noise_symbol_cons_get_gamma(pr, tmp, ptr->pnsym->index, abs);
		elina_interval_mul(tmp, tmp, ptr->coeff, ELINA_SCALAR_DOUBLE);
		elina_interval_add(box, box, tmp, ELINA_SCALAR_DOUBLE);
		if (p||q) {
		    /* continuing */
		    ptr->n = zonotope_aaterm_alloc_init();
		    ptr=ptr->n;
		}
	    }
	}
    }
	
    elina_interval_add(res->itv, exprA->itv, exprB->itv, ELINA_SCALAR_DOUBLE);
	
    elina_scalar_max(res->itv->inf, res->itv->inf, box->inf);
    elina_scalar_min(res->itv->sup, res->itv->sup, box->sup);
    elina_interval_free(box);
    elina_interval_free(tmp);
	
    return res;
}

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
    if ((size + 1) % 128 == 0) {
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

static inline zonotope_aff_t * zonotope_aff_from_linexpr0(zonotope_internal_t* pr, elina_linexpr0_t * expr, zonotope_t *z){
	size_t i;
	elina_dim_t dim;
	elina_coeff_t *coeff;
	zonotope_aff_t *res = zonotope_aff_alloc_init(pr);
	elina_coeff_t * cst = &(expr->cst);
	if(cst->discr==ELINA_COEFF_SCALAR){
		elina_interval_set_scalar(res->c, cst->val.scalar,cst->val.scalar);
        elina_interval_set_scalar(res->itv, cst->val.scalar,cst->val.scalar);
	}
	else{
		elina_interval_set(res->c, cst->val.interval);
        elina_interval_set(res->itv, cst->val.interval);
	}
	elina_linexpr0_ForeachLinterm(expr,i,dim,coeff) {
		zonotope_aff_t *aff = z->paf[dim];
    		if(coeff->discr==ELINA_COEFF_SCALAR){
			elina_scalar_t * scalar = elina_scalar_alloc();
			elina_scalar_set(scalar,coeff->val.scalar);
			elina_coeff_reinit(coeff,ELINA_COEFF_INTERVAL,ELINA_SCALAR_DOUBLE);
			elina_coeff_set_interval_scalar(coeff,scalar,scalar);
			elina_scalar_free(scalar);
    		}
    		elina_interval_t *interval = coeff->val.interval;
		zonotope_aff_t *tmp = zonotope_aff_mul_itv(pr,aff,interval);
		zonotope_aff_t *tmp1 = res;
      
		res = zonotope_aff_add(pr,tmp1,tmp,z);		
       
		zonotope_aff_free(pr,tmp);
		zonotope_aff_free(pr,tmp1);
  	}
	return res;
} 

static inline elina_linexpr0_t * elina_linexpr0_from_zonotope(zonotope_internal_t* pr, zonotope_aff_t * aff, zonotope_t *z){
	elina_linexpr0_t *res = elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE,aff->l);
	elina_linexpr0_set_cst_interval(res, aff->c);		
	uint_t k = 0;
	elina_dim_t dim;
	zonotope_aaterm_t *p;
	for(p=aff->q; p; p=p->n){
		elina_coeff_init(&res->p.linterm[k].coeff, ELINA_COEFF_INTERVAL);
		elina_coeff_set_interval(&res->p.linterm[k].coeff, p->coeff);
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
    elina_scalar_t *eps = elina_scalar_alloc();
    elina_scalar_t *err = elina_scalar_alloc();
    elina_interval_t *mid = elina_interval_alloc();
    elina_interval_t *dev = elina_interval_alloc();
    elina_interval_t *sum = elina_interval_alloc();
    elina_interval_t *itv = elina_interval_alloc();
    bool ok = false;
    //printf("aff\n");
    //zonotope_aff_fprint(pr,stdout,expr);
    if (elina_scalar_infty(expr->c->inf) || elina_scalar_infty(expr->c->sup)){
        ok = false;
        //printf("ok2\n");
    }
    else {
        elina_scalar_set_double(eps,5*1.11022302462515654042e-16);  /* threshold : last bit in the mantissa in double precision */
	elina_scalar_sub(err,expr->c->sup,expr->c->inf,ELINA_SCALAR_DOUBLE);
	if (elina_scalar_cmp(err, eps) > 0) {
	    elina_interval_middev(mid, dev, expr->c);
	    elina_interval_set(expr->c, mid);
	    elina_interval_add(sum, sum, dev, ELINA_SCALAR_DOUBLE);
	}
	for(p = expr->q; p; p=p->n) {
	    if (elina_scalar_infty(p->coeff->inf) || elina_scalar_infty(p->coeff->sup)) {
		if (elina_interval_is_top(p->coeff)) {
		    /* reduce to top */
		    if (expr->q) zonotope_aaterm_list_free(pr, expr->q);
		    expr->q = NULL;
		    expr->end = NULL;
		    expr->l = 0;
		    elina_interval_set_top(expr->c);
		    elina_interval_set_top(expr->itv);
		}
            //printf("ok1\n");
		ok = false;
		elina_interval_set_int(sum,0, ELINA_SCALAR_DOUBLE);
		break;
	    } else {
            //printf("ok3\n");
            //elina_interval_fprint(stdout,p->coeff);
		elina_scalar_sub(err,p->coeff->sup,p->coeff->inf,ELINA_SCALAR_DOUBLE);
		if (elina_scalar_cmp(err, eps) > 0) {
		    elina_interval_middev(mid, dev, p->coeff);
		    elina_interval_set(p->coeff, mid);
		    elina_interval_add(sum, sum, dev,ELINA_SCALAR_DOUBLE);
		}
	    }
	}
        //printf("sum\n");
        //elina_interval_fprint(stdout,sum);
	if (!elina_scalar_sgn(sum->inf) && !elina_scalar_sgn(sum->sup)) ok = false;
	else {
	    elina_scalar_set(itv->sup, sum->sup);
	    elina_scalar_set(itv->inf, itv->sup);
	    zonotope_aff_noise_symbol_create(pr, expr, itv, UN);
	    ok = true; /* adding new symbol */
	}
    }
   
    elina_interval_free(mid);
    elina_interval_free(dev);
    elina_interval_free(sum);
    elina_interval_free(itv);
    elina_scalar_free(eps);
    elina_scalar_free(err);
    return ok;
}

/* z1->itv and z2->itv may be different */
static inline bool zonotope_aff_is_eq(zonotope_internal_t* pr, zonotope_aff_t *z1, zonotope_aff_t *z2)
{
    if (z1 == z2) return true;
    else if (z1->l != z2->l) return false;
    else if (!elina_interval_equal(z1->c, z2->c)) return false;
    else {
	zonotope_aaterm_t *p, *q;
	for (p=z1->q, q=z2->q; p && q;) {
	    if (p->pnsym != q->pnsym) return false;
	    else if (!elina_interval_equal(p->coeff, q->coeff)) return false;
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
    if (elina_interval_cmp(zero,a) <= 0 || elina_interval_cmp(zero,b)<=0) {
	elina_interval_set_int(res,0,0);
	dir = 0;
    } else if (elina_scalar_sgn(a->inf)>=0 && elina_scalar_sgn(b->sup)<=0) {
	elina_interval_set_int(res,0,0);
	dir = 0;
    } else if (elina_scalar_sgn(a->sup)<=0 && elina_scalar_sgn(b->inf)>=0) {
	elina_interval_set_int(res,0,0);
	dir = 0;
    } else {
	/* a and b have the same sign */
	if (elina_scalar_sgn(a->inf)>=0) {
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
		if (size > 1){
			 z->nsymcons[size-1-1] = 0;
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
    elina_interval_t *nsymItv1 = elina_interval_alloc();
    elina_interval_t *nsymItv2 = elina_interval_alloc();
    elina_interval_t *pmptr = elina_interval_alloc();
    elina_interval_t *qmptr = elina_interval_alloc();
    elina_interval_t *argminpq = elina_interval_alloc();

    zonotope_aaterm_t *p, *q, *ptr;

    elina_scalar_min(res->itv->inf,exp1->itv->inf,exp2->itv->inf);
    elina_scalar_max(res->itv->sup,exp1->itv->sup,exp2->itv->sup);
    ptr = NULL;
    int s = 0;

    if (exp1->q || exp2->q) {
	elina_interval_set(c1, exp1->c);
	elina_interval_set(c2, exp2->c);
	ptr = zonotope_aaterm_alloc_init();
	for(p = exp1->q, q = exp2->q; p || q;) {
	    if (p && q) {
		if (p->pnsym->index == q->pnsym->index) {
		    zonotope_noise_symbol_cons_get_gamma(pr, nsymItv1, p->pnsym->index, z1);
		    zonotope_noise_symbol_cons_get_gamma(pr, nsymItv2, p->pnsym->index, z2);
		    if (p->pnsym->type == UN) {
			elina_interval_mul(tmp, nsymItv1, p->coeff, ELINA_SCALAR_DOUBLE);
			elina_interval_add(betaA, betaA, tmp, ELINA_SCALAR_DOUBLE);
			elina_interval_mul(tmp, nsymItv2, q->coeff, ELINA_SCALAR_DOUBLE);
			elina_interval_add(betaB, betaB, tmp, ELINA_SCALAR_DOUBLE);
		    } else if (p->pnsym->type == IN) {
			ptr->pnsym = p->pnsym;
			s = argmin(pr, ptr->coeff, p->coeff, q->coeff);
			if (s == -1) {
			    elina_interval_sub(qmptr, q->coeff, ptr->coeff, ELINA_SCALAR_DOUBLE);
			    //elina_interval_set_int(pmptr, 0,0);
			} else if (s == 1) {
			    //elina_interval_set_int(qmptr, 0,0);
			    elina_interval_sub(pmptr, p->coeff, ptr->coeff, ELINA_SCALAR_DOUBLE);
			} else {
			    elina_interval_set(qmptr,q->coeff);
			    elina_interval_set(pmptr,p->coeff);
			}
		    }
		    p = p->n ;
		    q = q->n ;
		} else if (p->pnsym->index < q->pnsym->index) {
		    zonotope_noise_symbol_cons_get_gamma(pr, nsymItv1, p->pnsym->index, z1);
		    zonotope_noise_symbol_cons_get_gamma(pr, nsymItv2, p->pnsym->index, z2);
		    if (p->pnsym->type == UN) {
			elina_interval_mul(tmp, nsymItv1, p->coeff,ELINA_SCALAR_DOUBLE);
			elina_interval_add(betaA, betaA, tmp, ELINA_SCALAR_DOUBLE);
		    } else if (p->pnsym->type == IN) {
			elina_interval_set(pmptr,p->coeff);
			elina_interval_set_int(qmptr,0,0);
		    }
		    zonotope_delete_constrained_noise_symbol(pr, p->pnsym->index, z3);
		    p = p->n;
		} else {
		    zonotope_noise_symbol_cons_get_gamma(pr, nsymItv1, q->pnsym->index, z1);
		    zonotope_noise_symbol_cons_get_gamma(pr, nsymItv2, q->pnsym->index, z2);
		    if (q->pnsym->type == UN) {
			elina_interval_mul(tmp, nsymItv2, q->coeff,ELINA_SCALAR_DOUBLE);
			elina_interval_add(betaB, betaB, tmp, ELINA_SCALAR_DOUBLE);
		    } else if (q->pnsym->type == IN) {
			elina_interval_set(qmptr,q->coeff);
			elina_interval_set_int(pmptr,0,0);
		    }
		    zonotope_delete_constrained_noise_symbol(pr, q->pnsym->index, z3);
		    q = q->n;
		}
	    } else if (p) {
		zonotope_noise_symbol_cons_get_gamma(pr, nsymItv1, p->pnsym->index, z1);
		zonotope_noise_symbol_cons_get_gamma(pr, nsymItv2, p->pnsym->index, z2);
		if (p->pnsym->type == UN) {
		    elina_interval_mul(tmp, nsymItv1, p->coeff,ELINA_SCALAR_DOUBLE);
		    elina_interval_add(betaA, betaA, tmp, ELINA_SCALAR_DOUBLE);
		} else if (p->pnsym->type == IN) {
		    elina_interval_set(pmptr,p->coeff);
		    elina_interval_set_int(qmptr,0,0);
		}
		zonotope_delete_constrained_noise_symbol(pr, p->pnsym->index, z3);
		p = p->n;
	    } else {
		zonotope_noise_symbol_cons_get_gamma(pr, nsymItv1, q->pnsym->index, z1);
		zonotope_noise_symbol_cons_get_gamma(pr, nsymItv2, q->pnsym->index, z2);
		if (q->pnsym->type == UN) {
		    elina_interval_mul(tmp, nsymItv2, q->coeff, ELINA_SCALAR_DOUBLE);
		    elina_interval_add(betaB, betaB, tmp, ELINA_SCALAR_DOUBLE);
		} else if (q->pnsym->type == IN) {
		    elina_interval_set(qmptr,q->coeff);
		    elina_interval_set_int(pmptr,0,0);
		}
		zonotope_delete_constrained_noise_symbol(pr, q->pnsym->index, z3);
		q = q->n;
	    }
	    elina_interval_mul(tmp1, nsymItv1, pmptr, ELINA_SCALAR_DOUBLE);
	    elina_interval_add(c1, c1, tmp1, ELINA_SCALAR_DOUBLE);
	    elina_interval_mul(tmp2, nsymItv2, qmptr, ELINA_SCALAR_DOUBLE);
	    elina_interval_add(c2, c2, tmp2, ELINA_SCALAR_DOUBLE);
	    elina_interval_set_int(pmptr,0,0);
	    elina_interval_set_int(qmptr,0,0);
	    if (!elina_scalar_sgn(ptr->coeff->inf) && !elina_scalar_sgn(ptr->coeff->sup)) {
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

	elina_interval_middev(mid, dev, betaA);
	elina_interval_set(d1, dev);
	elina_interval_add(c1, c1, mid, ELINA_SCALAR_DOUBLE);

	elina_interval_middev(mid, dev, betaB);
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
	//elina_scalar_neg(array[0],array[0]);
	elina_scalar_add(array[0],array[0],d1->sup, ELINA_SCALAR_DOUBLE);
	elina_scalar_set(array[1],c2->inf);
	//elina_scalar_neg(array[1],array[1]);
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
    elina_scalar_neg(array[0],array[0]);
	elina_scalar_add(array[0],array[0],d1->sup, ELINA_SCALAR_DOUBLE);
	elina_scalar_set(array[1],c2->inf);
    elina_scalar_neg(array[1],array[1]);
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
	elina_interval_set_scalar(beta,c0,c0);

	elina_scalar_sub(c0,dmax,dmin, ELINA_SCALAR_DOUBLE);
	elina_scalar_div_2(c0, c0);
	elina_interval_set_scalar(res->c,c0,c0);

	if (elina_scalar_cmp(c1->inf,c2->inf) < 0) {
	    elina_scalar_set(min, c2->inf);
	} else {
	    elina_scalar_set(min, c1->inf);
	}
	//elina_scalar_neg(min,min);

	if (elina_scalar_cmp(c1->sup,c2->sup) < 0) {
	    elina_scalar_set(max, c2->sup);
	} else {
	    elina_scalar_set(max, c1->sup);
	}

	if (elina_scalar_cmp(c0,min) <= 0 || elina_scalar_cmp(c0,max) >= 0) {
	    zonotope_aff_free(pr, res);
	    elina_interval_sub(tmp,d2,d1, ELINA_SCALAR_DOUBLE);
	    if (elina_scalar_sgn(tmp->inf)>=0) res = exp2;
	    else res = exp1;
	} else zonotope_aff_noise_symbol_create(pr, res, beta, UN);

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
	zonotope_aff_add_itv(pr, res, res->itv, UN);
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
	z->hypercube = true;
	elina_abstract0_free(pr->manNS, z->abs);
	z->abs = elina_abstract0_top(pr->manNS, 0,0);
	for (i=0; i<nsymcons_size; i++) {
	    if (z->gamma[i] != pr->ap_muu && z->gamma[i]) elina_interval_free(z->gamma[i]);
	    z->gamma[i] = NULL;
	}
	memset((void*)z->nsymcons, 0, nsymcons_size*sizeof(elina_dim_t));
    } else {
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
    elina_scalar_set_int(pr->muu->inf, (long int)-1);
    elina_scalar_set_int(pr->muu->sup, (long int)1);
    pr->ap_muu = elina_interval_alloc();
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
