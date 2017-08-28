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



#include "opt_oct_incr_closure_comp_sparse.h"

bool incremental_closure_comp_sparse(opt_oct_mat_t *oo,int dim, int v, bool is_int){
	double *m = oo->mat;
	array_comp_list_t *acl = oo->acl;
	int n = 2*dim;
	int ii = 2*v + 1;
	int j;
	double *temp1, *temp2;
  	temp1 = (double *)calloc(2*dim, sizeof(double));
	temp2 = (double *)calloc(2*dim, sizeof(double));
	unsigned short int *index1, *index2;
	int count = oo->nni;
	index1 = (unsigned short int *)calloc(2*(2*dim + 1),sizeof(unsigned short int));
	index2 = (unsigned short int *)calloc(2*(2*dim + 1),sizeof(unsigned short int));
	int size = 2*dim*(dim+1);
  	comp_list_t * cn = find(acl,v);
        /******
		Find the component set containing v, if it is null
		then we do not need to perform incremental closure.
	*******/
	if(cn!=NULL){
		//comp_list_t *cl = acl->head;
		//int l = 0;
		//while(cl!=cn){
		//	cl = cl->next;
			//l++;
		//}
		
		unsigned short int comp_size = cn->size;
		comp_t *ck = cn->head;
		/******
			incremental Floyd-Warshall : v in end-point position 
		******/
		for(unsigned k = 0; k < comp_size; k++){
			int k1 = 2*ck->num;
			ck = ck->next;
			int v1 = 2*v;
			int v2 = 2*v + 1;
			int v1v2 = v2 + (((v1 + 1)*(v1 + 1))/2);
			int v2v1 = v1 + (((v2 + 1)*(v2 + 1))/2);
			int kk1 = (k1^1);
			int br1 = k1 < v1 ? k1 : v1;
			int br2 = kk1 < v1 ? kk1 : v1;
			for(unsigned i = 2*v; i < 2*v + 2; i++){
				//double ik = m[n*i + k];
				int ind_ik, ind_ikk;
				if(k1 <=i){
					ind_ik = k1 + (((i + 1)*(i + 1))/2);
				}
				else{
					ind_ik = (i^1) + ((((k1^1) + 1)*((k1^1) + 1))/2);
				}
			
				if(kk1 <=i){
					ind_ikk = kk1 + (((i + 1)*(i + 1))/2);
				}
				else{
					ind_ikk = (i^1) + ((((kk1^1) + 1)*((kk1^1) + 1))/2);
				}
				double ik = m[ind_ik];
				double ikk = m[ind_ikk];
				//double ki = m[n*k + i];
				int ind_ki, ind_kki;
				if ( k1 <= i){
					ind_ki = (k1^1) + ((((i^1) + 1)*((i^1) + 1))/2);
				}
				else{
					ind_ki = i + (((k1 + 1)*(k1 + 1))/2);
				}

				if ( kk1 <= i){
					ind_kki = (kk1^1) + ((((i^1) + 1)*((i^1) + 1))/2);
				}
				else{
					ind_kki = i + (((kk1 + 1)*(kk1 + 1))/2);
				}
			
				//int ind_ki = i + (((k + 1)*(k + 1))/2);
				double ki = m[ind_ki];	
				double kki = m[ind_kki];
				/* v in first end-point position */
				if(ik != INFINITY){	
					comp_t * cj = cn->head;	
					for(j = 0; j <2*comp_size; j++){
						//double kj = m[n*k + j];
						int j1 = (j%2==0) ? 2*cj->num : 2*cj->num+1;
						if(j1 < v1){
							int ind_kj = opt_matpos2(k1,j1);
							double kj = m[ind_kj];
							//double jk = m[n*j + k];
							//int ind_jk = k + (((j + 1)*(j + 1))/2);
							//double jk = m[ind_jk];
							//m[n*i + j] = min(m[n*i + j], ik + kj);
							int ind_ij = j1 + (((i + 1)*(i + 1))/2);
							//if(m[ind_ij]!=INFINITY){
							m[ind_ij] = min(m[ind_ij], ik + kj);
						}
						//}
						//else{
						//m[ind_ij] = ik + kj;
							//count++;
						//}
						//m[n*j + i] = min(m[n*j + i], jk + ki);
						if(j%2==1){
							cj = cj->next;
						}
					}
				}

				/*if(ik != INFINITY){
					for(; j < v1; j++){
						//double kj = m[n*k + j];
						int ind_kj = (k1^1) + ((((j^1) + 1)*((j^1) + 1))/2);
						double kj = m[ind_kj];
						//double jk = m[n*j + k];
						int ind_ij = j + (((i + 1)*(i + 1))/2);
						//m[n*i + j] = min(m[n*i + j], ik + kj);
						//if(m[ind_ij]!=INFINITY){
						m[ind_ij] = min(m[ind_ij], ik + kj);
						//}
						//else{
						//	m[ind_ij] = ik + kj;
							//count++;
						//}
						//m[n*j + i] = min(m[n*j + i], jk + ki);
					}
				}*/
				/* v in second end-point position */
				if(ki != INFINITY){
					comp_t * cj = cn->head;	
					for(j= 0; j < 2*comp_size; j++ ){
						int j1 = (j%2==0) ? 2*cj->num : 2*cj->num+1;
						if(j1>=(2*v+2)){
							int ind_jk = opt_matpos2(j1,k1);
							double jk = m[ind_jk];
							int ind_ji = i + (((j1 + 1)*(j1 + 1))/2);
							//if(m[ind_ji]!=INFINITY){
							m[ind_ji] = min(m[ind_ji], jk + ki);
						}
						//}
						//else{
						//	m[ind_ji] = jk + ki;
							//count++;
						//}
						if(j%2==1){
							cj = cj->next;
						}
					}
				}

				/*if(ki != INFINITY){
					for(; j < 2*dim; j++){
						int ind_jk = k1 + (((j + 1)*(j + 1))/2);
						double jk = m[ind_jk];
						int ind_ji = i + (((j + 1)*(j + 1))/2);
						//if(m[ind_ji] != INFINITY){
						m[ind_ji] = min(m[ind_ji], jk + ki);
						//}
						//else{
						//	m[ind_ji] = jk + ki;
							//count++;
						//}
					}
				}*/
				/* v in first end-point position */
				if(ikk != INFINITY){
					comp_t * cj = cn->head;	
					for(j = 0; j <2*comp_size; j++){
						//double kj = m[n*k + j];
						int j1 = (j%2==0) ? 2*cj->num : 2*cj->num+1;
						if(j1 < v1){
							int ind_kkj = opt_matpos2(kk1,j1);
							double kkj = m[ind_kkj];
							//double jk = m[n*j + k];
							//int ind_jk = k + (((j + 1)*(j + 1))/2);
							//double jk = m[ind_jk];
							//m[n*i + j] = min(m[n*i + j], ik + kj);
							int ind_ij = j1 + (((i + 1)*(i + 1))/2);
							m[ind_ij] = min(m[ind_ij], ikk + kkj);
							//m[n*j + i] = min(m[n*j + i], jk + ki);
						}
						if(j%2==1){
							cj = cj->next;
						}
					}
				}
				/*if(ikk != INFINITY){
					for(; j < v1; j++){
						//double kj = m[n*k + j];
						int ind_kkj = (kk1^1) + ((((j^1) + 1)*((j^1) + 1))/2);
						double kkj = m[ind_kkj];
						//double jk = m[n*j + k];
						int ind_ij = j + (((i + 1)*(i + 1))/2);
						//m[n*i + j] = min(m[n*i + j], ik + kj);
						m[ind_ij] = min(m[ind_ij], ikk + kkj);
						//m[n*j + i] = min(m[n*j + i], jk + ki);
					}
				}*/
				/* v in second end-point position */
				if(kki != INFINITY){
					comp_t * cj = cn->head;	
					for(j= 0; j < 2*comp_size; j++ ){
						int j1 = (j%2==0) ? 2*cj->num : 2*cj->num+1;
						if(j1 >=(2*v+2)){
							int ind_jkk = opt_matpos2(j1,kk1);
							double jkk = m[ind_jkk];
							int ind_ji = i + (((j1 + 1)*(j1 + 1))/2);
							m[ind_ji] = min(m[ind_ji], jkk + kki);
						}
						if(j%2==1){
							cj = cj->next;
						}
					}
				}

				/*if(kki != INFINITY){
					for(; j < 2*dim; j++){
						int ind_jkk = kk1 + (((j + 1)*(j + 1))/2);
						double jkk = m[ind_jkk];
						int ind_ji = i + (((j + 1)*(j + 1))/2);
						m[ind_ji] = min(m[ind_ji], jkk + kki);
					}
				}*/
			
			}
			//if(m[v1v2]!=INFINITY){
			m[v1v2] = min(m[v1v2],m[opt_matpos2(v1,k1)] + m[opt_matpos2(k1,v2)]);
			//}
			//else{
				//m[v1v2] =  m[opt_matpos2(v1,k)] + m[opt_matpos2(k,v2)];
				//count++;
			//}
			m[v1v2] = min(m[v1v2],m[opt_matpos2(v1,kk1)] + m[opt_matpos2(kk1,v2)]);
			//if(m[v2v1] != INFINITY){
			m[v2v1] = min(m[v2v1],m[opt_matpos2(v2,k1)] + m[opt_matpos2(k1,v1)]);
			//}
			//else{
			//	m[v2v1] = m[opt_matpos2(v2,k)] + m[opt_matpos2(k,v1)];
				//count++;
			//}
			m[v2v1] = min(m[v2v1],m[opt_matpos2(v2,kk1)] + m[opt_matpos2(kk1,v1)]);
		 
		}
		
		int v1 = (2*v);
		int v2 = (2*v)^1;
		int vi = (((v1 + 1)*(v1 + 1))/2);
		int vvi = (((v2 + 1)*(v2 + 1))/2);
		int pos1 = v2 + vi;
		//int pos2 = opt_matpos2((2*k)^1, 2*k);
		int pos2 = v1 + vvi;
		//variable v in pivot position
		//if(m[pos1]!= INFINITY){
		unsigned short int * ca = to_sorted_array(cn,dim);
		//comp_t *ci = cn->head;
			for(int i = 0; i < 2*comp_size;i++){
				int i1 = (i%2==0)? 2*ca[i/2] : 2*ca[i/2]+1;
				if(i1 >= (v1+2)){
					int ind1 = v2 + (((i1+1)*(i1+1))/2);
					int ind2 = v1 + (((i1+1)*(i1+1))/2);
					m[ind1] = min(m[ind1], m[ind2] + m[pos1] );
					temp2[i1] = m[ind1];
				}
				if(temp2[i]!=INFINITY){
					count++;
				}
				
			}
		//}

		//if(m[pos2]!=INFINITY){
		
			for(int i = 0; i < 2*comp_size; i++){
				int i1 = (i%2==0)? 2*ca[i/2] : 2*ca[i/2]+1;
				if(i1 >=(v1+2)){
					int ind1 = v2 + (((i1+1)*(i1+1))/2);
					int ind2 = v1 + (((i1+1)*(i1+1))/2);
					m[ind2] = min(m[ind2], m[ind1] + m[pos2] );
					temp1[i1] = m[ind2];
				}
				if(temp1[i]!=INFINITY){
					count++;
				}
				
			}
		//}

		if(m[pos2]!=INFINITY){
			
			for(int j = 0; j < 2*comp_size; j++){
				//int ind3 = opt_matpos2((2*k)^1,j);
				int j1 = (j%2==0)? 2*ca[j/2] : 2*ca[j/2]+1;
				if(j1 < v1){
					int ind3 = j1 + vvi;
					//int ind4 = opt_matpos2( 2*k,j);
					int ind4 = j1 + vi;
					//result[n*((2*k)^1) + j] = min(result[n*((2*k)^1) + j], result[n*((2*k)^1) + 2*k] + result[n*(2*k) + j]);
					m[ind3] = min(m[ind3], m[pos2] + m[ind4]);
					if(m[ind3]!=INFINITY){
						count++;
					}
				}
				
				
			}
		}

		if(m[pos1] != INFINITY){
			
			for(int j = 0; j < 2*comp_size; j++){
				//int ind3 = opt_matpos2((2*k)^1,j);
				int j1 = (j%2==0) ? 2*ca[j/2] : 2*ca[j/2]+1;
				if(j1 < v1){
					int ind3 = j1 + vvi;
					//int ind4 = opt_matpos2(2*k,j);
					int ind4 = j1 + vi;
					//result[n*2*k + j] = min(result[n*2*k + j], result[n*2*k + ((2*k)^1)] + result[n*((2*k)^1) + j]);
					m[ind4] = min(m[ind4], m[pos1] + m[ind3]);
					if(m[ind4]!=INFINITY){
						count++;
					}
				}
				
				
			}
		}
		
		/******
			Compute the index
		******/
		compute_index_comp_sparse(m,ca,comp_size,index1,index2,v,dim);
		free(ca);
		//for(unsigned k = 2*v; k < 2*v + 2; k++){
		int ind1_k = index1[0];
		int ind2_k = index1[n + 1];
		int ind3_k = index2[0];
		int ind4_k = index2[n + 1];
		/* incremental Floyd-Warshall : v in pivot position */
		//First part	
		for(int i = 0; i < ind1_k; i++){
			int i1 = index1[i + 1];
			int i2 = (i1%2==0) ? (i1 + 1): i1;
			int br = i2 < 2*v ? i2 : 2*v - 1;
			int ind1 = i1 + ((((2*v) + 1)*((2*v) + 1))/2);
			//double t1 = result[n*(2*k) + i1];
			double t1 = m[ind1];
			//double t2 = result[n*((2*k)^1) + (i^1)];
			//int j2 = (j/2)*2;
			for(int j = 0;j < ind2_k ; j++){
				//int ind2 = get_index(k,j);
		        	//int j1 = index1[m*((2*k)^1) + j + 1];
				int j1 = index1[n + j + 2];
				if(j1 > br){
					break;
					//continue;
				}
				int ind2 = j1 + (((((2*v)^1) + 1)*(((2*v)^1) + 1))/2);
				//double op1 = t1 + result[n*((2*k)^1) + j1];
				double op1 = t1 + m[ind2];
				//double op2 = t2 + result[n*(2*k) + j];
				//double op3 = min(op1, op2);
				int ind3 = j1 + ((((i1^1) + 1)*((i1^1) + 1))/2);
				//result[n*(i1^1) + j1] = min(result[n*(i1^1) + j1],op1 );
				if(m[ind3]!=INFINITY){
					m[ind3] = min(m[ind3],op1 );
				}
				else{
					m[ind3] = op1;
					count++;
				}
			}
			//for(int j = 0; j < index2[m*2*k]; j++){
			for(int j = 0;j < ind3_k; j++){
				int j1 = index2[j + 1];
				if(j1>i2){
					 break;
					//continue;
			    	}
				double op1 = t1 + temp1[j1];
				//double op2 = t2 + temp2[j];
				//double op3 = min(op1, op2);
				int ind3 = (j1^1) + ((((i1^1) + 1)*((i1^1) + 1))/2);
				//result[n*(i1^1) + (j1^1)] = min(result[n*(i1^1) + (j1^1)],op1 );
				if(m[ind3] != INFINITY){
					m[ind3] = min(m[ind3],op1 );
				}
				else{
					m[ind3] = op1;
					count++;
				}
			}
			//}
  
		}
		
		//Second Part
		for(int i = 0; i < ind2_k; i++){
		    int i1 = index1[n + i + 2];
		    int i2 = (i1%2==0) ? (i1 + 1): i1;
		    int br = i2 < 2*v ? i2 : 2*v - 1;
		    //double t1 = result[n*(2*k) + i1];
		    int ind1 = i1 + (((((2*v)^1) + 1)*(((2*v)^1) + 1))/2);
		    //double t2 = result[n*((2*k)^1) + i1];
		    double t2 = m[ind1];
		    //int j2 = (j/2)*2;
		    for(int j = 0; j < ind1_k; j++){
			    int j1 = index1[j + 1];
			    if(j1 > br){
				break;
				//continue;
			    }
			    int ind2 = j1 + ((((2*v) + 1)*((2*v) + 1))/2);
		            //double op2 = t2 + result[n*(2*k) + j1];
			    double op2 = t2 + m[ind2];
			    int ind3 = j1 + ((((i1^1) + 1)*((i1^1) + 1))/2);
		            //result[n*(i1^1) + j1] = min(result[n*(i1^1) + j1],op2 );
			    if(m[ind3]!=INFINITY){
			    	m[ind3] = min(m[ind3],op2 );
			    }
			    else{
				m[ind3] = op2;
				count++;
			    }
		        }
		        //for(int j = 0; j < index2[m*((2*k)^1)]; j++){
			  for(int j = 0; j < ind4_k; j++){
			     int j1 = index2[n + j + 2];
			     if(j1>i2){
				break;
				//continue;
			    }
		            double op2 = t2 + temp2[j1];
			    int ind3 = (j1^1) + ((((i1^1) + 1)*((i1^1) + 1))/2);
		            //result[n*(i1^1) + (j1^1)] = min(result[n*(i1^1) + (j1^1)],op2 );
			    if(m[ind3]!=INFINITY){
			    	m[ind3] = min(m[ind3],op2 );
			    }
			    else{
				m[ind3] = op2;
				count++;
			    }
		        }
		    //}
		}
	    
		//Third Part
		for(int i = 0; i < ind4_k; i++){
		    int i1 = index2[n + i + 2];
		    int i2 = (i1%2==0) ? (i1 + 1): i1;
		    int br = i2 < 2*v ? i2 : 2*v - 1;
		    int ind1 = ((2*v)^1) + (((i1 + 1)*(i1 + 1))/2);
		    //double t1 = result[n*i1 + ((2*k)^1)];
		    double t1 = m[ind1];
		    //double t2 = result[n*((2*k)^1) + (i^1)];
		    //int j2 = (j/2)*2;
		   
		    for(int j = 0; j < ind2_k; j++){
		   	    //int j1 = index1[m*((2*k)^1) + j + 1];
			    int j1 = index1[n + j + 2];
			    if(j1 > br){
				break;
				//continue;
			    }
			    int ind2 = j1 + (((((2*v)^1) + 1)*(((2*v)^1) + 1))/2);
		            //double op1 = t1 + result[n*((2*k)^1) + j1];
			    double op1 = t1 + m[ind2];
			    int ind3 = j1 + (((i1 + 1)*(i1 + 1))/2);
		            //result[n*i1 + j1] = min(result[n*i1 + j1],op1 );
			    if(m[ind3]!=INFINITY){
			    	m[ind3] = min(m[ind3],op1 );
			    }
			    else{
				m[ind3] = op1;
				count++;
			    }
		        }
		        
		 	for(int j = 0; j < ind3_k ; j++){
		            
			    int j1 = index2[j + 1];
			    if(j1>i2){
				break;
				//continue;
			     }
			    //cout<<result[n*(j1^1) + 2*k]<<"\t"<<result[n*j1 + 2*k]<<"\n";
		            double op1 = t1 + temp1[j1];
		            //double op1 = t1 + result[n*(j1) + 2*k];
		     	    int ind3 = (j1^1) + (((i1 + 1)*(i1 + 1))/2);
		            //result[n*i1 + (j1^1)] = min(result[n*i1 + (j1^1)],op1 );
			    if(m[ind3]!=INFINITY){
			    	m[ind3] = min(m[ind3],op1 );
			    }
			  else{
				m[ind3] = op1;
				count++;
			    }
		        }
		    //}
		}
		
		//Fourth Part
		for(int i = 0; i < ind3_k; i++){
		    //int i1 = index2[m*(2*k) + i + 1];
		    int i1 = index2[i + 1];
		    int i2 = (i1%2==0) ? (i1 + 1): i1;
		     int br = i2 < 2*v ? i2 : 2*v - 1;
		    //double t1 = result[n*(2*k) + i1];
		    //double t2 = result[n*i1 + (2*k)];
		    int ind1 = (2*v) + (((i1 + 1)*(i1 + 1))/2);
		    double t2 = m[ind1];
		    for(int j = 0; j < ind1_k; j++){
			    int j1 = index1[j + 1];
			    if(j1 > br){
				break;
				//continue;
			    }
		            //double op2 = t2 + result[n*(2*k) + j1];
			    int ind2 = j1 + ((((2*v) + 1)*((2*v) + 1))/2);
			    double op2 = t2 + m[ind2];
		            //double op3 = min(op1, op2);
			  int ind3 = j1 + (((i1 + 1)*(i1 + 1))/2);
		            //result[n*i1 + j1] = min(result[n*i1 + j1],op2 );
			    if(m[ind3]!=INFINITY){
			    	m[ind3] = min(m[ind3],op2 );
			    }
			    else{
				m[ind3] = op2;
				count++;
			    }
		        }
		       
		        
		       for(int j = 0; j < ind4_k ; j++){
			    int j1 = index2[n + j + 2];
			    if(j1>i2){
				break;
				//continue;
			    }
		            double op2 = t2 + temp2[j1];
		            //double op3 = min(op1, op2);
			    int ind3 = (j1^1) + (((i1 + 1)*(i1 + 1))/2);
		            //result[n*i1 + (j1^1)] = min(result[n*i1 + (j1^1)],op2 );
			    if(m[ind3]!=INFINITY){
			    	m[ind3] = min(m[ind3],op2 );
			    }
			    else{
				m[ind3] = op2;
				count++;
			    }
		        }
		   // }
		}
	}
  	//incr_closure_time += cycles;
	int res;
	
	oo->nni = count;
	 
  	
	if(is_int){
		
		res = strengthning_int_comp_sparse(oo,index1,temp1, n);
		
	}
	else{
        	res = strengthning_comp_sparse(oo, index1, temp1, n);
	}
	free(index2);
	free(index1);
        free(temp1);
	free(temp2);
	return res;
}
