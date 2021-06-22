/*
 *
 *  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
 *  ELINA is Copyright © 2021 Department of Computer Science, ETH Zurich
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


/* ********************************************************************** */
/* opt_pk_lait.c: the learned Lait transformer */
/* ********************************************************************** */


#include "opt_pk_config.h"
#include "opt_pk_vector.h"
#include "opt_pk_matrix.h"
#include "opt_pk.h"
#include "opt_pk_representation.h"
#include "opt_pk_user.h"
#include "opt_pk_constructor.h"
#include "opt_pk_test.h"
#include "opt_pk_meetjoin.h"
#include "opt_pk_project.h"
#include "opt_pk_cherni.h"
#include "opt_pk_lait.h"
#include "python3.6/Python.h"


bool constraint_removal(elina_manager_t* man, opt_pk_t* poly, int* rmap) {
	opt_pk_internal_t* opk = opt_pk_init_from_manager(man, ELINA_FUNID_JOIN);
	opt_matrix_t* C = poly->C;

	if (!C) return false;

	for (size_t i=0; i<C->nbrows; i++)
		if (rmap[i] && C->p[i][0] == 0)
			--poly->nbeq;

	size_t nbrows = C->nbrows;
	size_t j=nbrows - 1;
	bool removed = false;
	for (size_t i=0; i<C->nbrows; i++) {
		if (rmap[i]) {
			removed = true;
			--nbrows;
			while(j>i && rmap[j]) --j;
			if (j>i) {
				opt_matrix_exch_rows(C, i, j);
			}
			--j;
		}
	}
	C->nbrows = nbrows;

	if (removed) {
		bool has_positivity_constraint = false;
		for (size_t i=0; i<C->nbrows; i++)
			if (opt_vector_is_positivity_constraint(opk, C->p[i], C->nbcolumns))
				has_positivity_constraint = true;
		if (!has_positivity_constraint) {
			++C->nbrows;
			opt_matrix_fill_constraint_top(opk, C, C->nbrows-1);
		}

		if (poly->F != NULL) {
			opt_matrix_free(poly->F);
			poly->F = NULL;
		}
		// if (poly->satF != NULL) {
		// 	opt_satmat_free(poly->satF);
		// 	poly->satF = NULL;
		// }
		// if (poly->satC != NULL) {
		// 	opt_satmat_free(poly->satC);
		// 	poly->satC = NULL;
		// }
		opt_poly_chernikova(man, poly, NULL);
	}

	return removed;
}

bool is_common_constraint(opt_pk_internal_t* opk, unsigned short int cons_maxcols, opt_numint_t *cons, comp_list_t *cons_cl, opt_pk_array_t *op){
	size_t i;
	array_comp_list_t *acl = op->acl;
	unsigned short int maxcols = op->maxcols;
	comp_list_t *cl = acl->head;
	unsigned short int k=0;
	while(cl!=NULL){
		if(is_included(cons_cl,cl,maxcols)) {
			unsigned short int *ca_cons_cl = to_sorted_array(cons_cl,cons_maxcols);
			unsigned short int *ca = to_sorted_array(cl,maxcols);
			opt_pk_t *poly = op->poly[k];
			opt_matrix_t *C = poly->C;
			size_t nbrows = C->nbrows;
			unsigned short int nbcolumns = C->nbcolumns;
			for(i=0; i < nbrows; i++){
				if((cons[0]!=C->p[i][0]) || (cons[1]!=C->p[i][1])){
					continue;
				}
				unsigned short int j,l=2;
				bool flag = true;
				for(j=2; j < nbcolumns; j++){
					if(ca_cons_cl[l-2]==ca[j-2]){
						if(cons[l]!=C->p[i][j]){
							flag = false;
							break;
						}
						l++;
					}
					else{
						if(C->p[i][j]!=0){
							flag = false;
							break;
						}
					}
					
				}
				if(flag){
					free(ca_cons_cl);
					free(ca);
					return true;
				}
			}
			free(ca_cons_cl);
			free(ca);
		}
		cl = cl->next;
		k++;
	}
	
	return false;
}

// python_path: the path contatining the python file opt_pk_lait.py
// model_path: the path to opt_pk_lait_model.pt
void opt_pk_lait_init(char* python_path, char* model_path) {
	char* original_python_path = getenv("PYTHONPATH");
	if (original_python_path == NULL) {
		setenv("PYTHONPATH", python_path, 1);
	} else {
		char new_python_path[200] = "";
		strcat(new_python_path, original_python_path);
		strcat(new_python_path, ":");
		strcat(new_python_path, python_path);
		setenv("PYTHONPATH", new_python_path, 1);
	}

    Py_Initialize();
    PyEval_InitThreads();

	PyGILState_STATE gstate = PyGILState_Ensure();
	PyObject* pName = PyUnicode_DecodeFSDefault("opt_pk_lait");
	PyObject* pModule = PyImport_Import(pName);
	PyObject* pFunc = PyObject_GetAttrString(pModule, "init");
	PyObject* pArgs = PyTuple_New(1);
	PyTuple_SetItem(pArgs, 0, PyBytes_FromString(model_path));
	PyObject* res = PyObject_CallObject(pFunc, pArgs);
	PyGILState_Release(gstate);
}

// oa: the first join input
// ob: the second join input
// res: the join output
// head: the abstract element at the loop head
// loop_iter: the loop iteration index
opt_pk_array_t* opt_pk_lait(elina_manager_t* man, bool destructive, opt_pk_array_t* oa, opt_pk_array_t* ob, opt_pk_array_t* res, opt_pk_array_t* head, int loop_iter) {
	res = destructive ? res : opt_pk_copy(man, res);

    opt_pk_internal_t* opk = opt_pk_init_from_manager(man, ELINA_FUNID_JOIN);
	comp_list_t* cl;

	opt_pk_t** poly = res->poly;
	array_comp_list_t* acl = res->acl;
	unsigned short int num_comp = acl->size;

	int maxcols = oa->maxcols > ob->maxcols? oa->maxcols : ob->maxcols;
	int blocks_a[maxcols];
	int blocks_b[maxcols];
	for (int i=0; i<maxcols; i++) blocks_a[i] = -1;
	for (int i=0; i<maxcols; i++) blocks_b[i] = -1;

    if (oa->acl != NULL) {
		cl = oa->acl->head;
		for (int i=0; i<oa->acl->size; i++) {
			comp_t* c = cl->head;
			while(c!=NULL) {
				blocks_a[c->num-opk->dec] = i;
				c = c->next;
			}
			cl = cl->next;
		}
    }

	if (ob->acl != NULL) {
		cl = ob->acl->head;
		for (int i=0; i<ob->acl->size; i++) {
			comp_t* c = cl->head;
			while(c!=NULL) {
				blocks_b[c->num - opk->dec] = i;
				c = c->next;
			}
			cl = cl->next;
		}
	}

    int num_cons = 0;
	for (int i=0; i<num_comp; i++) {
		opt_matrix_t* C = poly[i]->C;
		if (C)
			for (size_t j=0; j<C->nbrows; j++)
				if (!opt_vector_is_positivity_constraint(opk, C->p[j], C->nbcolumns))
					++num_cons;
	}

    unsigned int vars_head[head->maxcols-opk->dec];
    int num_vars_head = 0;
    array_comp_list_t* head_acl = head->acl;
    if (head_acl != NULL) {
        cl = head_acl->head;
        while (cl != NULL) {
            comp_t* c = cl->head;
            while (c != NULL) {
                vars_head[num_vars_head++] = c->num-opk->dec;
                c = c->next;
            }
            cl = cl->next;
        }
    }

	int num_features = 12;
	int features[num_cons][num_features];

	int index_cons = 0;
	cl = acl->head;

	for (int i=0; i<num_comp; i++) {
		opt_matrix_t* C = poly[i]->C;

		if (C) {
			unsigned short int* ca = to_sorted_array(cl, res->maxcols);

			int num_vars = 0;
			int is_var_head[cl->size];
			int var_cons[cl->size];
			for (int j=0; j<cl->size; j++) {
				is_var_head[j] = 0;
				var_cons[j] = 0;
				for (size_t k=0; k<C->nbrows; k++) {
					if (C->p[k][j+opk->dec] != 0) {
						++num_vars;
						for (int l=0; l<num_vars_head; l++) {
							if (vars_head[l] == ca[j] - opk->dec) {
								is_var_head[j] = 1;
								break;
							}
						}
						break;
					}
				}
				for (size_t k=0; k<C->nbrows; k++)
					if (C->p[k][j+opk->dec] != 0)
						++var_cons[j];
				if (var_cons[j] > 0) --var_cons[j];
			}

			bool has_positivity = false;
			for (size_t j=0; j<C->nbrows; j++)
				if (opt_vector_is_positivity_constraint(opk, C->p[j], C->nbcolumns))
					has_positivity = true;

			for(size_t j=0; j<C->nbrows; j++) {
				if (opt_vector_is_positivity_constraint(opk, C->p[j], C->nbcolumns)) continue;

				features[index_cons][0] = num_vars;
				features[index_cons][1] = has_positivity? C->nbrows - 1 : C->nbrows;
				features[index_cons][2] = poly[i]->F ? poly[i]->F->nbrows : 0;
				features[index_cons][3] = C->p[j][0] == 0;
				features[index_cons][4] = 0; 
				features[index_cons][5] = 0; 
				features[index_cons][6] = 0; 
				features[index_cons][7] = is_common_constraint(opk, res->maxcols, C->p[j], cl, oa) && is_common_constraint(opk, res->maxcols, C->p[j], cl, ob);
				features[index_cons][8] = 0;
				features[index_cons][9] = 0;
				features[index_cons][10] = 0;
				features[index_cons][11] = is_common_constraint(opk, res->maxcols, C->p[j], cl, head);

				int block_a = -1;
				int block_b = -1;
				bool block_a_diff = false;
				bool block_b_diff = false;
				for (int k=0; k<cl->size; k++) {
					unsigned short var = ca[k] - opk->dec;
					opt_numint_t coeff = C->p[j][k+opk->dec];
					if (coeff == 0) continue;

					features[index_cons][4]++;
					if (is_var_head[k]) features[index_cons][10]++;
					if (coeff > 100 || coeff <= -100) features[index_cons][5] += 1;
					if (block_a != -1 && blocks_a[var] != -1 && block_a != blocks_a[var]) {
						block_a_diff = true;
					}
					if (block_b != -1 && blocks_b[var] != -1 && block_b != blocks_b[var]) {
						block_b_diff = true;
					}
					block_a = blocks_a[var];
					block_b = blocks_b[var];
					features[index_cons][8] += var_cons[k];
				}

				if (block_a_diff && block_b_diff) {
					features[index_cons][6] = 1;
				}

				features[index_cons][9] = features[index_cons][4] == features[index_cons][10] ? 1 : 0;

				++index_cons;
			}

			free(ca);
		}

		cl = cl->next;
	}

	PyGILState_STATE gstate = PyGILState_Ensure();

	PyObject* pName = PyUnicode_DecodeFSDefault("opt_pk_lait");
	PyObject* pModule = PyImport_Import(pName);
	PyObject* pFunc = PyObject_GetAttrString(pModule, "lait");
	PyObject* pArgs = PyTuple_New(3);

	PyObject* pBlockLens = PyList_New(0);
	for (int i=0; i<num_comp; i++) {
		opt_matrix_t* C = poly[i]->C;
		if (C) {
			int num = 0;
			for (size_t j=0; j<C->nbrows; j++)
				if (!opt_vector_is_positivity_constraint(opk, C->p[j], C->nbcolumns))
					++num;
			PyList_Append(pBlockLens, PyLong_FromLong(num));
		}
	}
	PyTuple_SetItem(pArgs, 0, pBlockLens);

	PyObject* pFeatures = PyList_New(num_cons);
	for (int i=0; i<num_cons; i++) {
		PyObject* pFeature = PyList_New(num_features);
		for (int j=0; j<num_features; j++)
			PyList_SetItem(pFeature, j, PyLong_FromLong(features[i][j]));
		PyList_SetItem(pFeatures, i, pFeature);
	}
	PyTuple_SetItem(pArgs, 1, pFeatures);

	PyObject* pEdges = PyList_New(0);
	index_cons = 0;
	for (int i=0; i<num_comp; i++) {
		opt_matrix_t* C = poly[i]->C;
		if (C) {
			int indices[C->nbrows];
			int vars_in_cons[C->nbrows][C->nbcolumns-opk->dec];
			for (size_t j=0; j<C->nbrows; j++) {
				if (!opt_vector_is_positivity_constraint(opk, C->p[j], C->nbcolumns)) {
					indices[j] = index_cons;
					for (size_t k=0; k<C->nbcolumns-opk->dec; k++) {
						if (C->p[j][k+opk->dec] != 0) {
							vars_in_cons[j][k] = 1;
						} else {
							vars_in_cons[j][k] = 0;
						}
					}
					++index_cons;
				} else {
					indices[j] = -1;
				}
			}

			for (size_t j=0; j<C->nbrows; j++) {
				if (indices[j] == -1) continue;
				for (size_t k=0; k<C->nbrows; k++) {
					if (indices[k] == -1) continue;
					int edge_type = 0;
					for (size_t l=0; l<C->nbcolumns-opk->dec; l++)
						if (vars_in_cons[j][l] && vars_in_cons[k][l])
							++edge_type;
					if (edge_type > 0) {
						PyObject* pEdge = PyList_New(3);
						PyList_SetItem(pEdge, 0, PyLong_FromLong(indices[j]));
						PyList_SetItem(pEdge, 1, PyLong_FromLong(indices[k]));
						PyList_SetItem(pEdge, 2, PyLong_FromLong(edge_type));
						PyList_Append(pEdges, pEdge);
					}
				}
			}
		}
	}
	PyTuple_SetItem(pArgs, 2, pEdges);

    PyObject* pRemove = PyObject_CallObject(pFunc, pArgs);

	index_cons = 0;
	cl = acl->head;
	for (int i=0; i<num_comp; i++) {
		opt_matrix_t* C = poly[i]->C;

		if (C) {
			int rmap[C->nbrows];
			for (size_t j=0; j<C->nbrows; j++) {
				if (opt_vector_is_positivity_constraint(opk, C->p[j], C->nbcolumns)) {
					rmap[j] = 0;
					continue;
				}

				rmap[j] = PyLong_AsLong(PyList_GetItem(pRemove, index_cons));
				index_cons++;
			}

			if (!constraint_removal(man, poly[i], rmap)) {
				cl = cl->next;
				continue;
			}

			if(opk->exn){
				for (size_t j=0; j<C->nbrows; j++) {
					if (opt_vector_is_positivity_constraint(opk, C->p[j], C->nbcolumns)) {
						rmap[j] = 0;
						continue;
					}

					int num_vars = 0;
					for (int k=0; k<cl->size; k++)
						if (C->p[j][k+opk->dec] != 0)
							++num_vars;

					if (num_vars > 3) {
						rmap[j] = 1;
					}
				}
				opk->exn = ELINA_EXC_NONE;
				if (!constraint_removal(man, poly[i], rmap) || opk->exn) {
					poly[i]->C->nbrows = 1;
					opt_matrix_fill_constraint_top(opk, poly[i]->C, 0);
					opt_poly_chernikova(man, poly[i], NULL);
					opk->exn = ELINA_EXC_NONE;
				}
			}
		}

		cl = cl->next;
	}

	PyGILState_Release(gstate);

	return res;
}