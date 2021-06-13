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

#include "opt_oct_hmat.h"
#include "opt_oct_lait.h"
#include "python3.6/Python.h"

void compute_partition(opt_oct_mat_t* oo, int dim, int *num_comp, int *num_cons_arg, int **map, int **cm, int **cons_count, int **var_count, int **var_cons) {
  *map = (int*)malloc(dim * sizeof(int));
  memset(*map, 0, dim * sizeof(int));
  *cm = (int*)malloc(dim * sizeof(int));
  memset(*cm, 0, dim * sizeof(int));
  *var_cons = (int*)malloc(dim * sizeof(int));
  memset(*var_cons, 0, dim * sizeof(int));
  int num_cons = 0;

  if (!oo->is_dense) {
    array_comp_list_t* acl = oo->acl;
    *num_comp = acl->size;
    *cons_count = (int*)malloc((*num_comp+1) * sizeof(int));
    memset(*cons_count, 0, (*num_comp+1) * sizeof(int));

    int part = 1;
    comp_list_t *cl = acl->head;
    while (cl != NULL) {
      comp_t *c = cl->head;
      while (c != NULL) {
        unsigned short num = c->num;
        (*map)[num] = 1;
        (*cm)[num] = part;
        c = c->next;
      }
      ++part;
      cl = cl->next;
    }

    for (int i=0; i<2*dim; i++) {
      if (!(*map)[i/2]) continue;
      for (int j=0; j<=(i|1); j++) {
        if (!(*map)[j/2]) continue;
        if ((*cm)[j/2] != (*cm)[i/2]) continue;
        if (i == j || oo->mat[opt_matpos2(i, j)] == INFINITY) continue;
        ++num_cons;
        (*cons_count)[(*cm)[j/2]] += 1;
        (*var_cons)[i/2] += 1;
        (*var_cons)[j/2] += 1;
      }
    }
  } else {
    *num_comp = 1;
    *cons_count = (int*)malloc((*num_comp+1) * sizeof(int));
    memset(*cons_count, 0, (*num_comp+1) * sizeof(int));

    for (int i=0; i<2*dim; i++) {
      for (int j=0; j<=(i|1); j++) {
        if (i == j || oo->mat[opt_matpos2(i, j)] == INFINITY) continue;
        (*map)[i/2] = 1;
        (*map)[j/2] = 1;
        (*cm)[i/2] = 1;
        (*cm)[j/2] = 1;
        (*var_cons)[i/2] += 1;
        (*var_cons)[j/2] += 1;
        ++num_cons;
      }
    }
    (*cons_count)[1] = num_cons;
  }

  *var_count = (int*)malloc((*num_comp+1) * sizeof(int));
  memset(*var_count, 0, (*num_comp+1) * sizeof(int));
  for (int i=0; i<dim; i++) {
    if ((*map)[i]) {
      (*var_count)[(*cm)[i]] += 1;
    }
  }

  *num_cons_arg = num_cons;
}

array_comp_list_t * compute_finest(opt_oct_mat_t * oo, unsigned short int dim){
	array_comp_list_t *res = create_array_comp_list();
	comp_list_t * cl = oo->acl->head;
	double * m = oo->mat;
	while(cl!=NULL){
		unsigned short int comp_size = cl->size;
		unsigned short int *ca = to_sorted_array(cl,dim);
		unsigned short int i,j;
		for(i=0; i < comp_size; i++){
			unsigned short int i1 = ca[i];
			if(m[opt_matpos2(2*i1,2*i1+1)]!=INFINITY || m[opt_matpos2(2*i1+1,2*i1)]!=INFINITY){
				comp_list_t * cl = create_comp_list();
				insert_comp(cl,i1);
				insert_comp_list(res,cl);
			}
		}
		for(i=0; i < comp_size; i++){
			unsigned short int i1 = ca[i];
			for(j=0; j < i; j++){
				unsigned short int j1 = ca[j];
				if(m[opt_matpos2(2*i1,2*j1+1)]!=INFINITY || m[opt_matpos2(2*i1+1,2*j1)]!=INFINITY || m[opt_matpos2(2*i1,2*j1)]!=INFINITY || m[opt_matpos2(2*i1+1,2*j1+1)]!=INFINITY){
					comp_list_t * cl = create_comp_list();
					insert_comp(cl,i1);
					insert_comp(cl,j1);
					insert_comp_list_with_union(res,cl,dim);
				}
			}
		}
		free(ca);
		cl = cl->next;
	}
	return res;
}

// python_path: the path contatining the python file opt_oct_lait.py
// model_path: the path to opt_oct_lait_model.pt
void opt_oct_lait_init(char* python_path, char* model_path) {
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
	PyObject* pName = PyUnicode_DecodeFSDefault("opt_oct_lait");
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
opt_oct_t* opt_oct_lait(elina_manager_t* man, bool destructive, opt_oct_t* o1, opt_oct_t* o2, opt_oct_t* o_res, opt_oct_t* o_head, int loop_iter) {
    o_res = destructive ? o_res : opt_oct_copy(man, o_res);

    if (!o1->closed && !o1->m) return o_res;
    if (!o2->closed && !o2->m) return o_res;

    opt_oct_mat_t* oo1 = o1->closed ? o1->closed : o1->m;
    opt_oct_mat_t* oo2 = o2->closed ? o2->closed : o2->m;
    opt_oct_mat_t* oo_res = o_res->closed ? o_res->closed : o_res->m;
    int dim = o1->dim;

    int *map, *cm, *cons_count, *var_count, *var_cons;
    int num_cons, num_comp;
    compute_partition(oo_res, dim, &num_comp, &num_cons, &map, &cm, &cons_count, &var_count, &var_cons);

    int *map1, *cm1, *cons_count1, *var_count1, *var_cons1;
    int num_cons1, num_comp1;
    compute_partition(oo1, dim, &num_comp1, &num_cons1, &map1, &cm1, &cons_count1, &var_count1, &var_cons1);

    int *map2, *cm2, *cons_count2, *var_count2, *var_cons2;
    int num_cons2, num_comp2;
    compute_partition(oo2, dim, &num_comp2, &num_cons2, &map2, &cm2, &cons_count2, &var_count2, &var_cons2);

    PyGILState_STATE gstate = PyGILState_Ensure();

    PyObject* pName = PyUnicode_DecodeFSDefault("opt_oct_lait");
    PyObject* pModule = PyImport_Import(pName);
    PyObject* pFunc = PyObject_GetAttrString(pModule, "lait");
	PyObject* pArgs = PyTuple_New(3);

    PyObject* pBlockLens = PyList_New(0);
    for (int i=0; i<num_comp; i++) {
        PyList_Append(pBlockLens, PyLong_FromLong(cons_count[i+1]));
    }
    PyTuple_SetItem(pArgs, 0, pBlockLens);

    int num_cons_predict = 0;
    int cons_dim[num_cons][2];
    if (!oo_res->is_dense) {
        for (int i=0; i<2*dim; i++) {
            if (!map[i/2]) continue;
            for (int j=0; j<=(i|1); j++) {
                if (!map[j/2]) continue;
                if (cm[j/2] != cm[i/2]) continue;
                if (i == j || oo_res->mat[opt_matpos2(i, j)] == INFINITY) continue;
                if (i/2 == j/2) continue;
                if (cm1[i/2] == cm1[j/2] || cm2[i/2] == cm2[j/2]) continue;
    
                cons_dim[num_cons_predict][0] = i;
                cons_dim[num_cons_predict][1] = j;
                num_cons_predict++;
            }
        }
    } else {
        for (int i=0; i<2*dim; i++) {
            for (int j=0; j<=(i|1); j++) {
                if (i == j || oo_res->mat[opt_matpos2(i, j)] == INFINITY) continue;
                if (i/2 == j/2) continue;
                if (cm1[i/2] == cm1[j/2] || cm2[i/2] == cm2[j/2]) continue;

                cons_dim[num_cons_predict][0] = i;
                cons_dim[num_cons_predict][1] = j;
                num_cons_predict++;
            }
        }
    }

    opt_oct_mat_t* oo_head = o_head->closed? o_head->closed : o_head->m;
    int *map_h, *cm_h, *cons_count_h, *var_count_h, *var_cons_h, num_comp_h, num_cons_h;
    compute_partition(oo_head, o_head->dim, &num_comp_h, &num_cons_h, &map_h, &cm_h, &cons_count_h, &var_count_h, &var_cons_h);
    int* is_var_head = (int*)calloc(dim, sizeof(int));
    for (int i=0; i<o_head->dim; i++) {
        if (map_h[i]) is_var_head[i] = 1;
    }

    int num_features = 12;
    int features[num_cons_predict][num_features];

    for (int k=0; k<num_cons_predict; k++) {
        int i = cons_dim[k][0];
        int j = cons_dim[k][1];

        int pos = opt_matpos2(i, j);
        features[k][0] = cons_count[cm[i/2]];
        features[k][1] = var_count[cm[i/2]];
        features[k][2] = num_cons_predict;
        features[k][3] = fabs(oo_res->mat[pos]);
        features[k][4] = 0;
        features[k][5] = var_cons[i/2];
        features[k][6] = var_cons[j/2];
        if (oo_res->mat[opt_matpos2(i/2*2+1, i/2*2)] != INFINITY) {
            if (oo_res->mat[opt_matpos2(i/2*2, i/2*2+1)] != INFINITY) {
                features[k][7] = 2;
            } else {
                features[k][7] = 1;
            }
        } else {
            if (oo_res->mat[opt_matpos2(i/2*2, i/2*2+1)] != INFINITY) {
                features[k][7] = 1;
            } else {
                features[k][7] = 0;
            }
        }
        if (oo_res->mat[opt_matpos2(j/2*2+1, j/2*2)] != INFINITY) {
            if (oo_res->mat[opt_matpos2(j/2*2, j/2*2+1)] != INFINITY) {
                features[k][8] = 2;
            } else {
                features[k][8] = 1;
            }
        } else {
            if (oo_res->mat[opt_matpos2(j/2*2, j/2*2+1)] != INFINITY) {
                features[k][8] = 1;
            } else {
                features[k][8] = 0;
            }
        }
        features[k][9] = (int)(oo_res->mat[opt_matpos2(i^1, j^1)] != INFINITY);
        features[k][4] = loop_iter;
        features[k][10] = is_var_head[i/2] + is_var_head[j/2];
        features[k][11] = oo_head->mat[pos] != INFINITY && oo_head->mat[pos] == oo_res->mat[pos];
    }

    PyObject* pFeatures = PyList_New(0);
    for (int i=0; i<num_cons_predict; i++) {
        PyObject* pFeature = PyList_New(num_features);
        for (int j=0; j<num_features; j++) {
            PyList_SetItem(pFeature, j, PyLong_FromLong(features[i][j]));
        }
        PyList_Append(pFeatures, pFeature);
    }
    PyTuple_SetItem(pArgs, 1, pFeatures);

    int var_to_cons[2*dim][num_cons_predict];
    int var_to_cons_count[2*dim];
    for (int i=0; i<2*dim; i++) var_to_cons_count[i] = 0;

    for (int k=0; k<num_cons_predict; k++) {
        int i = cons_dim[k][0];
        int j = cons_dim[k][1];

        var_to_cons[i][var_to_cons_count[i]] = k;
        ++var_to_cons_count[i];
        var_to_cons[j][var_to_cons_count[j]] = k;
        ++var_to_cons_count[j];
    }

    PyObject* pEdges = PyList_New(0);
    for (int i=0; i<2*dim; i++) {
        for (int j=0; j<var_to_cons_count[i]; j++) {
            for (int k=0; k<var_to_cons_count[i]; k++) {
                PyObject* pEdge = PyList_New(3);
                PyList_SetItem(pEdge, 0, PyLong_FromLong(var_to_cons[i][j]));
                PyList_SetItem(pEdge, 1, PyLong_FromLong(var_to_cons[i][k]));
                PyList_SetItem(pEdge, 2, PyLong_FromLong(1));
                PyList_Append(pEdges, pEdge);
            }
        }
    }
    PyTuple_SetItem(pArgs, 2, pEdges);

    PyObject* pRemove = PyObject_CallObject(pFunc, pArgs);

    for (int k=0; k<num_cons_predict; k++) {
        if (PyLong_AsLong(PyList_GetItem(pRemove, k)) == 1) {
            int i = cons_dim[k][0];
            int j = cons_dim[k][1];
            oo_res->mat[opt_matpos2(i, j)] = INFINITY;
        }
    }
    if (!oo_res->is_dense) {
        oo_res->acl = compute_finest(oo_res, dim);
    }

    free(is_var_head);
    free(map);
    free(cm);
    free(cons_count);
    free(var_count);
    free(var_cons);
    free(map1);
    free(cm1);
    free(cons_count1);
    free(var_count1);
    free(var_cons1);
    free(map2);
    free(cm2);
    free(cons_count2);
    free(var_count2);
    free(var_cons2);
    free(map_h);
    free(cm_h);
    free(cons_count_h);
    free(var_count_h);
    free(var_cons_h);

    PyGILState_Release(gstate);

    return o_res;
}