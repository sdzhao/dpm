#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "R.h"

void dpm(double * Sigma, double * omg, int *dd, int *diags, double * lambda, int *nnlambda, double * gamma, int * max_ite, double * prec, double * e_i)
{
    int i,j,k,m,dim,junk_a,size_a,size_a_pre,w_idx,rs_idx,nlambda;
    int ite_ext,ite_int1,ite_int2,gap_ext,max_ite1,max_ite2,max_ite3, tmp_m;
    double gap_int,ilambda,tmp1,tmp2,err1,err2,omg_2norm,mu_2norm,omg_dif,eps1,eps2,omg_temp,alp_dif,mu_dif,max_dif;
    double omg_sum1, omg_sum2, omg_dif_sum1, omg_out1;

    dim = *dd;
    //dim_sq = dim*dim;
    nlambda = *nnlambda;
    double *omg0 = (double*) malloc(dim*sizeof(double));
    double *omg1 = (double*) malloc(dim*sizeof(double));
    int *idx_a = (int*) malloc(dim*sizeof(int)); //sizes of active sets
    int *idx_i = (int*) malloc(dim*sizeof(int)); //sizes of inactive sets
    double *alp_tild = (double*) malloc(dim*sizeof(double));
    double *mu_grad = (double*) malloc(dim*sizeof(double));
    double *omg_grad = (double*) malloc(dim*sizeof(double));
    double *omg_pre = (double*) malloc(dim*sizeof(double));
    double *alp = (double*) malloc(dim*sizeof(double));
    double *mu = (double*) malloc(dim*sizeof(double));
    //double *e_i = (double*) malloc(dim*sizeof(double)); // comment out because I'm SETTING e_i in the argument
    double *S_col = (double*) malloc(dim*sizeof(double));
    double *gamma_col = (double*) malloc(dim*sizeof(double));
    double *S_omg = (double*) malloc(dim*sizeof(double));
    double *omg_tild = (double*) malloc(dim*sizeof(double));
    double *r = (double*) malloc(dim*sizeof(double));
    double *y_i = (double*) malloc(dim*sizeof(double));
    double *SS = (double*) malloc(dim*dim*sizeof(double));
    double *Sy = (double*) malloc(dim*sizeof(double));

    max_ite1 = * max_ite;
    max_ite2 = 1e2;
    max_ite3 = 1e2;
    eps1 = * prec;
    eps2 = 1e-3;
    //cnz = 0;
    for(i=0; i<dim; i++){
        S_col[i] = 0;
        for(j=0; j<dim; j++){
            SS[j*dim+i] = 0;
            for(k=0; k<dim; k++)
                SS[j*dim+i] += Sigma[k*dim+i]*Sigma[j*dim+k];
        }
        S_col[i] = SS[i*dim+i];
        gamma_col[i] = *gamma/S_col[i];
    }
    //for(i=0; i<dim; i++){
    i=0; // only want first column
        for(j=0; j<dim; j++) {
  	    //e_i[j] = 0; // comment out because I'm SETTING e_i in the argument
            omg0[j] = 0;
            idx_i[j] = 1;
            alp[j] = 0;
            mu[j] = 0;
        }
        // idx_i[i] = 0;
        // e_i[i] = 1; // comment out because I'm SETTING e_i in the argument
        //for(j=0; j<dim; j++)
	// {
	//   printf("%f\n",e_i[j]); // see if e is correct
	//  }
        size_a = 0;

        for(m=0; m<nlambda; m++) {
            /* for(j=0; j<dim; j++) { */
            /*     //alp[j] = 0; */
            /*     //mu[j] = 0; */
            /*     //omg0[j] = 0; */
            /*     //idx_i[j] = 1; */
            /* } */
            //size_a = 0;
            ite_ext = 0;
            //tmp_m = m*dim_sq+i*dim;
	    tmp_m=m*dim; // i'm only concerned about one column
            max_dif = 1;
            while(max_dif > eps1 && ite_ext < max_ite1){
	      //while(ite_ext < max_ite1){
                // update alpha
                for(j=0; j<dim; j++){
                    S_omg[j]=0;
                    for(k=0; k<size_a; k++){
                        w_idx = idx_a[k];
                        S_omg[j]+=Sigma[w_idx*dim+j]*omg0[w_idx];
                    }
                    alp_tild[j]=e_i[j]-S_omg[j]-mu[j];
                }
                alp_dif = 0;
                for(j=0; j<dim; j++){
		  // change lambda for diagonal elements
		  ilambda = lambda[m];
		  if(diags[j]==1){
		    ilambda = lambda[m]/2;
		  }
                    if (alp_tild[j]<=-ilambda){
                        alp[j]=-ilambda;
                        alp_dif = fabs(alp_tild[j]+ilambda)>alp_dif ? fabs(alp_tild[j]+ilambda) : alp_dif;
                    }
                    else {
                        if (alp_tild[j]>=ilambda){
                            alp[j] = ilambda;
                            alp_dif = fabs(alp_tild[j]-ilambda)>alp_dif ? fabs(alp_tild[j]-ilambda) : alp_dif;
                        }
                        else
                            alp[j] = alp_tild[j];
                    }
                }

                // update omega
                for(j=0; j<dim; j++)
                    y_i[j] = e_i[j]-alp[j]-mu[j];
                for(j=0; j<dim; j++){
                    Sy[j] = 0;
                    omg_pre[j] = omg0[j];
                    for(k=0; k<dim; k++){
                        Sy[j] += Sigma[j*dim+k]*y_i[k];
                    }
                }

                gap_ext = 1;
                ite_int1 = 0;
                while(gap_ext !=0 && ite_int1<max_ite2){
                    size_a_pre = size_a;
                    for(j=0; j<dim; j++){
                        if(idx_i[j] == 1){
                            omg_tild[j] = 0;
                            for(k=0; k<size_a; k++){
                                w_idx = idx_a[k];
                                omg_tild[j] += SS[w_idx*dim+j]*omg0[w_idx];
                            }
                            omg_tild[j] = (Sy[j]-omg_tild[j]+S_col[j]*omg0[j])/S_col[j];
                            if(fabs(omg_tild[j])<=gamma_col[j]) {
                                omg1[j] = 0;
                            }
                            else{
                                if(omg_tild[j]>gamma_col[j])
                                    omg1[j] = omg_tild[j] - gamma_col[j];
                                else
                                    omg1[j] = omg_tild[j] + gamma_col[j];
                                idx_a[size_a] = j;
                                size_a++;
                                idx_i[j] = 0;
                            }
                            omg0[j] = omg1[j];
                        }
                    }
                    gap_ext = size_a - size_a_pre;
                    gap_int = 1;
                    ite_int2 = 0;
                    while(gap_int>eps2 && ite_int2<max_ite3){
                        tmp1 = 0;
                        tmp2 = 0;
                        for(j=0; j<size_a; j++){
                            w_idx = idx_a[j];
                            omg_tild[w_idx] = 0;
                            for(k=0; k<size_a; k++){
                                rs_idx = idx_a[k];
                                omg_tild[w_idx] += SS[rs_idx*dim+w_idx]*omg0[rs_idx];
                            }
                            omg_tild[w_idx] = (Sy[w_idx]-omg_tild[w_idx]+S_col[w_idx]*omg0[w_idx])/S_col[w_idx];
                            if (fabs(omg_tild[w_idx]) <= gamma_col[w_idx]) {
                                omg1[w_idx] = 0;
                            }
                            else {
                                if (omg_tild[w_idx]>gamma_col[w_idx])
                                    omg1[w_idx] = omg_tild[w_idx] - gamma_col[w_idx];
                                else
                                    omg1[w_idx] = omg_tild[w_idx] + gamma_col[w_idx];
                                tmp2 = tmp2+fabs(omg1[w_idx]);
                            }
                            omg_dif = omg1[w_idx]-omg0[w_idx];
                            tmp1 = tmp1+fabs(omg_dif);
                            omg0[w_idx] = omg1[w_idx];
                        }
                        gap_int = tmp1/tmp2;
                        ite_int2++;
                    }
                    //ite_cnt_int2[m*dim+i] += ite_int2;

                    junk_a = 0;
                    for(j=0; j<size_a; j++){
                        w_idx = idx_a[j];
                        if (omg1[w_idx] == 0){
                            junk_a++;
                            idx_i[w_idx] = 1;
                        }
                        else
                            idx_a[j-junk_a] = w_idx;
                    }
                    size_a = size_a - junk_a;
                    ite_int1++;
                }
                //ite_cnt_int1[m*dim+i] += ite_int1;

                omg_sum1 = 0;
                omg_dif_sum1 = 0;
                omg_sum2 = 0;
                for(j=0; j<dim; j++){
                    //omg_dif_sum1 += fabs(omg_pre[j] - omg0[j]);
                    omg_sum1 += fabs(omg0[j]);
                    omg_sum2 += fabs(omg_pre[j]);
                    S_omg[j]=0;
                    for(k=0; k<size_a; k++){
                        w_idx = idx_a[k];
                        S_omg[j]+=Sigma[w_idx*dim+j]*omg1[w_idx];
                    }
                }
                //omg_out1 = omg_dif_sum1/omg_sum1;
                omg_out1 = fabs(omg_sum2-omg_sum1)/omg_sum1;

                // update mu
                mu_dif = 0;
                for(j=0; j<dim; j++){
                    mu_grad[j]=alp[j]+S_omg[j]-e_i[j];
                    mu[j] += mu_grad[j];
                    mu_dif = fabs(mu_grad[j])>mu_dif ? fabs(mu_grad[j]) : mu_dif;
                }
                max_dif = omg_out1>mu_dif ? omg_out1 : mu_dif;
//if(ite_ext%10==0 && i%100==0)
//printf("ite_ext=%d,omg_out1=%f,omg_dif_sum1=%f,omg_sum1=%f,mu_dif=%f \n",ite_ext,omg_out1,omg_dif_sum1,omg_sum1,mu_dif);
//printf("ite_ext=%d,obj_dif=%f,omg_sum1=%f,mu_dif=%f \n",ite_ext,fabs(omg_sum2-omg_sum1),omg_sum1,mu_dif);
                ite_ext++;
            }
            //ite_cnt_ext[m*dim+i] = ite_ext;
	    
            for(j=0;j<size_a;j++){
                w_idx = idx_a[j];
                omg[tmp_m+w_idx] = omg1[w_idx];
                /* if(w_idx != i) { */
		/*     x[cnz] = omg1[w_idx]; */
                /*     row_idx[cnz] = m*dim+w_idx; */
                /*     cnz++; */
                /* } */
            }
        }
        //col_cnz[i+1]=cnz;
    //}
    free(omg0);
    free(omg1);
    free(idx_a);
    free(idx_i);
    free(alp_tild);
    free(mu_grad);
    free(omg_grad);
    free(omg_pre);
    free(alp);
    free(mu);
    //free(e_i); // comment out because I'm SETTING e_i in the argument
    free(S_col);
    free(gamma_col);
    free(S_omg);
    free(omg_tild);
    free(r);
    free(y_i);
    free(SS);
    free(Sy);
}
