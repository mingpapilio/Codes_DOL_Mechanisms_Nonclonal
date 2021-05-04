/* 
 * This file runs time series of evolving proportion of helper/ level of coordination.
 * Average property of populations are printed from the simulations.
 * Please check comments on parameters for their descriptions.
 ************************************************************************************************
 * # This section is example bash code for gcc compiler #
 * # Execution #
 * gcc-9 hsta_rts.c -O2
 * ./a.out
 * # Please make sure the location of dSFMT folder (random number generator) 
 * and include-instruction match each other #
 ************************************************************************************************
 * Output files:
 * out_q: target proportion of helper (genetic trait)
 * out_porp: proportion of helpers in the population (phenotype)
 * out_SD_porp: standard deviation of the proportion of helpers
 * out_s: level of coordination (genetic trait)
 ************************************************************************************************
 * Key switches and parameters of simulation:
 * S1: Whether coordination is evolving
 * S2: Type of coordination cost
 * S3: Whether target proportion of helper is evolving
 * n: size of the social group
 * e: essentiality of helping
 * N_ini_idvl: the number of founders in a group, notated as 'l' in main text
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "../dSFMT-src-2.2.3/dSFMT.c"

// Functions
    double Mean_array(double p[], int length);                                  // Calculating average of a list
    double SD_array(double p[], int length);                                    // Calculating standard deviation of a list
    double sum(double p[], int length);                                         // Sum of a list
    void swap(double *p1, double *p2);                                          // swapping two elements
    int RandFromProb(double p[], int length, double RandNum);                   // Random sampling with weighted probability
    double normal_dist_BM (double mean, double sd, double u1, double u2);       // Normal-distributed random number generator
    double **d_array2d(long size_1, long size_2);                               // Creating a 2-D array in heap mamory
    void free_d_array2d(double **x);                                            // Releasing a memory containing a 2-D array
    double ***d_array3d(long size_1, long size_2, long size_3);                 // Creating a 3-D array in heap mamory
    void free_d_array3d(double ***x);                                           // Releasing a memory containing a 3-D array
    double ****d_array4d(long size_1, long size_2, long size_3, long size_4);   // Creating a 4-D array in heap mamory
    void free_d_array4d(double ****x);                                          // Releasing a memory containing a 4-D array

// Main codes
int main(){
    // Switches
        int s1= 2;                      // 0: No coordination (fixed at 0), 1: full coordinaiton (fixed at 1), 2: evolving coordination
        int s2= 1;                      // 0: Linear cost of coordination, 1: diminishing (nonlinear) cost
        int s3= 0;                      // 0: Evolvable target proportion of helper (q), 1: fixed q= q*_FR, 2: q fixed at any given value

    // General parameters
        int T= 1E5;                     // Number of generations
        int N_rep= 10;                  // Number of repetitions
        int n= 18;                      // Number of individual in each social group
        int N_alloc= 5*n;               // Number of allocating divisions 
        int N_ini_idvl= 2;              // Number of initial group members (Relatedness= its inverse)
        int N_total= 10000;             // Total number of individual (approximate number, actual number is determined by G)
        int G= ceil((double)N_total/n); // Number of groups
        double theta= 0.025;            // Coefficient of coordination cost
        double e= 0.5+ 0.5*(10 -1)/9;   // Trait essentiality, related to fitness return function
        double lambda= (double)(n-1)/n; // Trait sociality (assuming equal share of helping accross group members here)
        double mut_rate= 0.01;          // Mutation rate
        double MutStep= 0.1;            // Size of mutation for traits

    // Random number genertor
        int seed;
        dsfmt_t dsfmt;
        seed= time(NULL);
        if(seed==0)seed= 1;
        dsfmt_init_gen_rand(&dsfmt,seed);

    // Variables
        int i,j,k,i_1,i_2;
        int t, rep, idx, num_div_idvl;
        double pp, pp2, pp3, u1, u2, rho, cost;
        double avg_Q, asd_Q, avg_H1, asd_H1, avg_H2, asd_H2, avg_S, asd_S, avg_porp, asd_porp, sd_avg_Q, sd_avg_H1, sd_avg_H2, sd_asd_Q, sd_asd_H1, sd_asd_H2, sd_avg_S, sd_asd_S, sd_avg_porp, sd_asd_porp;
        double p_i;                // proportion of interacting division 1 individuals
        int n_i;                   // number of interacting neughbours
        int k_i;                   // number of interacting division 1 individuals

    // Parameters for repeats
        int length_param= 10;
        double essential_list[length_param];
        for(i_1=0; i_1< length_param; ++i_1){
            essential_list[i_1]= 0.5+ i_1*(0.5)/(length_param-1);
        }

    // Indexes
        int id_q=   0;                  // Target proportion of helper of focal individual
        int id_h1=  1;                  // Level of helping in division 1
        int id_h2=  2;                  // Level of helping in division 2
        int id_s=   3;                  // Probability og knowing the division of other group members
        int id_type=4;                  // Phenotype of the individual (0: undetermined, 1: div 1, 2: div 2)
    // Temporal space
        double tmp[G];                  // Temporary array to store a property of the group members within the focal social group
        double idvl_tmp[n];             // Temporary array for individual levels of coordination
        double ****Popn= d_array4d(N_rep, G, n, 5);     // The population matrix
        double ****PopnNext= d_array4d(N_rep, G, n, 5); // The population in next generation, sampled from the offsprings
        double **Poptmp= d_array2d(N_rep, G*n);         // Temporal space for taking population average
        double **ProbIdvl= d_array2d(N_rep,G*n);        // Probability of each group to be chosen to survive (global competition)
        double **PorpDiv1= d_array2d(N_rep, G);
        double ****link= d_array4d(N_rep, G, n, n);     // Linkage matrix of each group
        //
        double globalQ[N_rep], globalQ_sd[N_rep], globalH1[N_rep], globalH1_sd[N_rep], globalH2[N_rep], globalH2_sd[N_rep], globalS[N_rep], globalS_sd[N_rep], globalDiv1[N_rep], globalDiv1_sd[N_rep], globalQ_avg[N_rep];
        //
        double log_ini_idvl[N_rep][N_ini_idvl][4];
        int k_num[N_rep];
        int samp, g_samp, n_samp;
        double social_temp[N_rep], social_wofocal[N_rep], mut_size[N_rep], q_temp[N_rep], h1_temp[N_rep], h2_temp[N_rep], s_temp[N_rep];

    // Output file
        FILE *out_q, *out_s, *out_porp, *out_SD_porp;
        out_q= fopen("out_q_10.txt","w");
        fprintf(out_q, "T\tp_1\tp_2\tp_3\tp_4\tp_5\tp_6\tp_7\tp_8\tp_9\tp_10\n");
        out_s= fopen("out_s_10.txt","w");
        fprintf(out_s, "T\tp_1\tp_2\tp_3\tp_4\tp_5\tp_6\tp_7\tp_8\tp_9\tp_10\n");
        out_porp= fopen("out_porp_10.txt", "w");
        fprintf(out_porp, "T\tp_1\tp_2\tp_3\tp_4\tp_5\tp_6\tp_7\tp_8\tp_9\tp_10\n");
        out_SD_porp= fopen("out_SD_porp_10.txt", "w");
        fprintf(out_SD_porp, "T\tp_1\tp_2\tp_3\tp_4\tp_5\tp_6\tp_7\tp_8\tp_9\tp_10\n");
    
    // Initialization
        t=0;
        for(rep=0; rep< N_rep; ++rep){
            globalQ_avg[rep]= 0.0;
            for(i=0; i<G; ++i){
                for(j=0; j<n; ++j){
                    // Assigning initial target proportion of helpers
                    if(s3== 0) Popn[rep][i][j][id_q]=   0.5;
                    if(s3== 1) Popn[rep][i][j][id_q]=   (n*(e-1)+e*(n/N_ini_idvl-1))/(e*(n+n/N_ini_idvl-2));
                    if(s3== 2) Popn[rep][i][j][id_q]=   0.4;
                    // Degree of helping in two divisions- assuming div 1 sterile helper, div 2 pure reproductive
                    Popn[rep][i][j][id_h1]=   1.0;
                    Popn[rep][i][j][id_h2]=   0.0;
                    // Assigning phenotype to undetermined
                    Popn[rep][i][j][id_type]= 0.0;
                    //
                    if(s1== 1) Popn[rep][i][j][id_s]=   1.0;
                    else Popn[rep][i][j][id_s]=         0.0;
        }}}
        // Start of simulation
        for(t=1; t<=T; ++t){
            for(rep=0; rep<N_rep; rep++){
            // Establish links, allocate divisions of labour, calculate fitness
                for(i=0; i<G; ++i){
                // Establishing the coordination matrix (full coordination)
                    if(s1==1){
                        for(j=0; j<n; ++j){
                            for(k=0; k<n; ++k){
                                if(j!=k) link[rep][i][j][k]= 1.0;
                                else link[rep][i][j][k]= NAN;
                    }}}
                    // Establishing the coordination matrix (evoving coordination)
                    if(s1==2){
                        for(j=0; j<n; ++j){
                            for(k=0; k<n; ++k){
                                if(j!=k){
                                    if(dsfmt_genrand_open_close(&dsfmt)< Popn[rep][i][j][id_s]) \
                                        link[rep][i][j][k]= 1.0;
                                    else link[rep][i][j][k]= 0.0;
                                }
                                else link[rep][i][j][k]= NAN;
                            }
                    }}
                // Assigning the division of labour
                    // No coordination
                    if(s1==0){
                        for(j=0; j<n; ++j){
                            // Become division 1
                            if(dsfmt_genrand_open_close(&dsfmt)< Popn[rep][i][j][id_q]) Popn[rep][i][j][id_type]= 1.0;
                            // Become division 2
                            else Popn[rep][i][j][id_type]= 2.0;
                        }
                    }
                    // Full and intermediate coordination
                    if(s1==1 || s1==2){
                        for(j=0; j<N_alloc; ++j){
                            n_i= k_i= 0;
                            // Pick a random individual from the group
                            n_samp= floor(dsfmt_genrand_open_close(&dsfmt)*n);
                            // Calculate the provability of being a helper from all connected neighbours
                            for(k=0; k<n; ++k){
                                if(n_samp!=k && link[rep][i][n_samp][k]== 1.0){
                                    n_i+= 1;
                                    if(Popn[rep][i][k][id_type]==1.0) k_i+= 1;
                                }
                            }
                            // Conditional allocation
                            if(n_i> 0){
                                // Calculate observed proportion of helpers
                                p_i= (double)k_i/ n_i;
                                // If observed is lower than its target
                                if(Popn[rep][i][n_samp][id_q]> p_i) Popn[rep][i][n_samp][id_type]= 1.0;
                                if(Popn[rep][i][n_samp][id_q]< p_i) Popn[rep][i][n_samp][id_type]= 2.0;
                                if(Popn[rep][i][n_samp][id_q]== p_i){
                                    if(dsfmt_genrand_open_close(&dsfmt)< Popn[rep][i][n_samp][id_q]) Popn[rep][i][n_samp][id_type]= 1.0;
                                    else Popn[rep][i][n_samp][id_type]= 2.0;
                                }
                            }
                            else{
                                if(dsfmt_genrand_open_close(&dsfmt)< Popn[rep][i][n_samp][id_q]) \
                                    Popn[rep][i][n_samp][id_type]= 1.0;
                                else Popn[rep][i][n_samp][id_type]= 2.0;
                            }
                        }
                    }
                // Fitness calculation
                        k_num[rep]=0;
                        social_temp[rep]= 0.0;
                        // Count the number of division 1
                        for(j=0; j<n; ++j){
                            if(Popn[rep][i][j][id_type]==1.0){
                                k_num[rep]+= 1;
                                social_temp[rep]+= Popn[rep][i][j][id_h1];
                            }
                            else social_temp[rep]+= Popn[rep][i][j][id_h2];
                        }
                        // Calculate the actual proportion of division 1
                        PorpDiv1[rep][i]= (double)k_num[rep]/n;
                        // Calcualte the fitness of each individual
                        for(j=0; j<n; ++j){
                            // Calculate the amount of public good excluding the focal individual (i.e. social_wofocal)
                            if(Popn[rep][i][j][id_type]== 1.0) social_wofocal[rep]= social_temp[rep]- Popn[rep][i][j][id_h1];
                            else social_wofocal[rep]= social_temp[rep]- Popn[rep][i][j][id_h2];
                            // Assigning the cost of coordination
                            if(s2== 0) cost= theta*Popn[rep][i][j][id_s];
                            if(s2== 1) cost= theta*(1-exp(-5.0*Popn[rep][i][j][id_s]));
                            if(cost> 1.0) cost= 1.0;
                            // Calcualting fitness
                            // Division type 1
                            if(Popn[rep][i][j][id_type]==1.0){
                                ProbIdvl[rep][i*n+j]= (1- cost)*(1-Popn[rep][i][j][id_h1])*(1-e+    \
                                    e*( (1-lambda)*Popn[rep][i][j][id_h1]+ lambda/(n-1)*social_wofocal[rep] ));
                            }
                            // Division type 2
                            else{
                                ProbIdvl[rep][i*n+j]= (1- cost)*(1-Popn[rep][i][j][id_h2])*(1-e+    \
                                    e*( (1-lambda)*Popn[rep][i][j][id_h2]+ lambda/(n-1)*social_wofocal[rep] ));
                            }
                        }
                }

            // Mutation, forming the next population
                for(i=0; i<G; ++i){
                    for(j=0; j<N_ini_idvl; ++j){
                        samp= RandFromProb(ProbIdvl[rep], G*n, dsfmt_genrand_open_close(&dsfmt));
                        g_samp= samp/n;
                        n_samp= samp%n;
                        // Mutation of q
                            if((s3==0) && dsfmt_genrand_open_close(&dsfmt)< mut_rate){
                                mut_size[rep]= MutStep*normal_dist_BM(0,1, \
                                    dsfmt_genrand_open_close(&dsfmt),dsfmt_genrand_open_close(&dsfmt));
                            }
                            else mut_size[rep]= 0.0;
                            q_temp[rep]= Popn[rep][g_samp][n_samp][id_q]+ mut_size[rep];
                        // Assuming no mutation in h1 and h2
                            h1_temp[rep]= Popn[rep][g_samp][n_samp][id_h1];
                            h2_temp[rep]= Popn[rep][g_samp][n_samp][id_h2];
                        // Mutation of s
                            if(s1==2){
                                if(dsfmt_genrand_open_close(&dsfmt)< mut_rate){
                                    mut_size[rep]= MutStep*normal_dist_BM(0,1, \
                                        dsfmt_genrand_open_close(&dsfmt),dsfmt_genrand_open_close(&dsfmt));
                                }
                                else mut_size[rep]= 0.0;
                                s_temp[rep]= Popn[rep][g_samp][n_samp][id_s]+ mut_size[rep];
                            }
                        // Dealing with the boundary conditions
                            if(q_temp[rep]> 1) q_temp[rep]= 1.0;
                            if(q_temp[rep]< 0) q_temp[rep]= 0.0;
                            if(h1_temp[rep]> 1) h1_temp[rep]= 1.0;
                            if(h1_temp[rep]< 0) h1_temp[rep]= 0.0;
                            if(h2_temp[rep]> 1) h2_temp[rep]= 1.0;
                            if(h2_temp[rep]< 0) h2_temp[rep]= 0.0;
                            if(s_temp[rep]> 1) s_temp[rep]= 1.0;
                            if(s_temp[rep]< 0) s_temp[rep]= 0.0;
                        // logging
                            log_ini_idvl[rep][j][id_q]= q_temp[rep];
                            log_ini_idvl[rep][j][id_h1]= h1_temp[rep];
                            log_ini_idvl[rep][j][id_h2]= h2_temp[rep];
                            log_ini_idvl[rep][j][id_s]= s_temp[rep];
                    }
                    // Assigning the next generation
                        // Clonal group, filling all group member for the next generation with the same values
                        if(N_ini_idvl== 1){
                            for(j=0; j<n; j++){
                                PopnNext[rep][i][j][id_q]= q_temp[rep];
                                PopnNext[rep][i][j][id_h1]= h1_temp[rep];
                                PopnNext[rep][i][j][id_h2]= h2_temp[rep];
                                PopnNext[rep][i][j][id_type]= 0.0;
                                if(s1==2) PopnNext[rep][i][j][id_s]= s_temp[rep];
                                else if(s1==0) PopnNext[rep][i][j][id_s]= 0.0;
                                else if(s1==1) PopnNext[rep][i][j][id_s]= 1.0;
                        }}
                        // Non-clonal group
                        else{
                            num_div_idvl= N_ini_idvl*(int)(n/N_ini_idvl);
                            for(j=0; j<n; j++){
                                if(j<= num_div_idvl){
                                    for(k=0; k<N_ini_idvl; k++){
                                        if(j%N_ini_idvl== k){
                                            PopnNext[rep][i][j][id_q]= log_ini_idvl[rep][k][id_q];
                                            PopnNext[rep][i][j][id_h1]= log_ini_idvl[rep][k][id_h1];
                                            PopnNext[rep][i][j][id_h2]= log_ini_idvl[rep][k][id_h2];
                                            PopnNext[rep][i][j][id_type]= 0.0;
                                            if(s1==2) PopnNext[rep][i][j][id_s]= log_ini_idvl[rep][k][id_s];
                                            else if(s1==0) PopnNext[rep][i][j][id_s]= 0.0;
                                            else if(s1==1) PopnNext[rep][i][j][id_s]= 1.0;
                                }}}
                                else{
                                    k= floor(dsfmt_genrand_open_close(&dsfmt)*N_ini_idvl);
                                    PopnNext[rep][i][j][id_q]= log_ini_idvl[rep][k][id_q];
                                    PopnNext[rep][i][j][id_h1]= log_ini_idvl[rep][k][id_h1];
                                    PopnNext[rep][i][j][id_h2]= log_ini_idvl[rep][k][id_h2];
                                    PopnNext[rep][i][j][id_type]= 0.0;
                                    if(s1==2) PopnNext[rep][i][j][id_s]= log_ini_idvl[rep][k][id_s];
                                    else if(s1==0) PopnNext[rep][i][j][id_s]= 0.0;
                                    else if(s1==1) PopnNext[rep][i][j][id_s]= 1.0;
                                }
                }}}
            }
            // Calculate population values
            if(t%10==0){
                for(rep=0; rep<N_rep; rep++){
                    for(i=0; i<G; i++){
                        for(j=0; j<n; j++){
                            for(k=0; k<5; k++) Popn[rep][i][j][k]= PopnNext[rep][i][j][k];
                    }}
                    if(N_ini_idvl==1){
                        for(i=0; i<G; i++) tmp[i]= Popn[rep][i][0][id_q];
                        globalQ[rep]= Mean_array(tmp, G);
                        globalQ_sd[rep]= SD_array(tmp, G);
                        for(i=0; i<G; i++) tmp[i]= Popn[rep][i][0][id_h1];
                        globalH1[rep]= Mean_array(tmp, G);
                        globalH1_sd[rep]= SD_array(tmp, G);
                        for(i=0; i<G; i++) tmp[i]= Popn[rep][i][0][id_h2];
                        globalH2[rep]= Mean_array(tmp, G);
                        globalH2_sd[rep]= SD_array(tmp, G);
                        for(i=0; i<G; i++) tmp[i]= Popn[rep][i][0][id_s];
                        globalS[rep]= Mean_array(tmp, G);
                        globalS_sd[rep]= SD_array(tmp, G);
                    }
                    else{
                        // Q
                        for(i=0; i<G; ++i){
                            for(j=0; j<n; ++j) Poptmp[rep][i*n+j]= Popn[rep][i][j][id_q];
                        }
                        globalQ[rep]= Mean_array(Poptmp[rep], G*n);
                        globalQ_sd[rep]= SD_array(Poptmp[rep], G*n);
                        // H1
                        for(i=0; i<G; ++i){
                            for(j=0; j<n; ++j) Poptmp[rep][i*n+j]= Popn[rep][i][j][id_h1];
                        }
                        globalH1[rep]= Mean_array(Poptmp[rep], G*n);
                        globalH1_sd[rep]= SD_array(Poptmp[rep], G*n);
                        // H2
                        for(i=0; i<G; ++i){
                            for(j=0; j<n; ++j) Poptmp[rep][i*n+j]= Popn[rep][i][j][id_h2];
                        }
                        globalH2[rep]= Mean_array(Poptmp[rep], G*n);
                        globalH2_sd[rep]= SD_array(Poptmp[rep], G*n);
                        // S
                        for(i=0; i<G; ++i){
                            for(j=0; j<n; ++j) Poptmp[rep][i*n+j]= Popn[rep][i][j][id_s];
                        }
                        globalS[rep]= Mean_array(Poptmp[rep], G*n);
                        globalS_sd[rep]= SD_array(Poptmp[rep], G*n);
                    }
                    // Also calculate the actual divisions
                        for(i=0; i<G; i++) tmp[i]= PorpDiv1[rep][i];
                        globalDiv1[rep]= Mean_array(tmp, G);
                        globalDiv1_sd[rep]= SD_array(tmp, G);
                    // Swapping the trait ids, making trait 1 always does more helping
                    if((globalH2[rep]-globalH1[rep])> 0.1){
                        swap(&globalH1[rep], &globalH2[rep]);
                        swap(&globalH1_sd[rep], &globalH2_sd[rep]);
                        globalQ[rep]= 1-globalQ[rep];
                        globalDiv1[rep]= 1-globalDiv1[rep];
                    }
                }
            // Print
                fprintf(out_q, "%d\t", t);
                fprintf(out_s, "%d\t", t);
                fprintf(out_porp, "%d\t", t);
                fprintf(out_SD_porp, "%d\t", t);
                for(rep=0; rep<N_rep; rep++){
                    if(rep<(N_rep-1)){
                        fprintf(out_q, "%lf\t", globalQ[rep]);
                        fprintf(out_s, "%lf\t", globalS[rep]);
                        fprintf(out_porp, "%lf\t", globalDiv1[rep]);
                        fprintf(out_SD_porp, "%lf\t", globalDiv1_sd[rep]);
                    }
                    if(rep==(N_rep-1)) {
                        fprintf(out_q, "%lf\n", globalQ[rep]);
                        fprintf(out_s, "%lf\n", globalS[rep]);
                        fprintf(out_porp, "%lf\n", globalDiv1[rep]);
                        fprintf(out_SD_porp, "%lf\n", globalDiv1_sd[rep]);
                    }
        }}}
    // Release memory and termination
    free_d_array4d(Popn);
    free_d_array4d(PopnNext);
    free_d_array2d(Poptmp);
    free_d_array2d(ProbIdvl);
    free_d_array2d(PorpDiv1);
    free_d_array4d(link);
    fclose(out_q);
    fclose(out_s);
    fclose(out_porp);
    fclose(out_SD_porp);
    return 0;
}
////////////////////// Functions are defined in below //////////////////////
double Mean_array(double p[], int length){
    int i;
    double temp;
    temp= 0.0;
    if(length> 0){
        for (i=0; i<length; ++i) temp+= p[i]/length;
        return temp;
    }
    else return NAN;
}

double SD_array(double p[], int length){
    if(length> 1){
        int i;
        double avg, temp;
        avg= temp= 0.0;
        for (i=0; i<length; ++i) avg+= p[i]/length;
        for (i=0; i<length; ++i) temp+= pow(p[i]-avg, 2);
        return sqrt(temp/(length-1));
    }
    else return NAN;
}

double sum(double p[], int length){
    int i;
    double sum=0.0;
    for(i=0; i< length; ++i) sum+= p[i];
    return sum;
}

void swap(double *p1, double *p2){
    double temp= *p1;
    *p1= *p2;
    *p2= temp;
}

int RandFromProb(double p[], int length, double RandNum){
    int i,temp,check;
    double sum, pp;
    sum= pp= 0.0;
    for(i=0; i< length; ++i) sum+= p[i];
    if(sum> 0.0){
        for(i=0; i< length; ++i){
            if(RandNum> pp/sum) check=1;
            pp+= p[i];
            if(RandNum< pp/sum&& check==1){
                temp= i;
                i= length;
            }
            check=0;
        }
        return temp;
    }
    else{
        return floor(length*RandNum);
    }
}

double normal_dist_BM (double mean, double sd, double u1, double u2){
    // Using Box-Muller method to generate pseudo-normal distributed numbers in [0,1]
    double z1;
	z1= sqrt(-2* log(u1))* cos(2* M_PI* u2);
	return z1*sd+ mean;
}

double **d_array2d(long size_1, long size_2){
    double **x;
    long i;
    x= (double **) malloc((size_t)(size_1*sizeof(double)));
    for(i=0; i< size_1; ++i) x[i]= (double *) malloc((size_t)(size_2*sizeof(double)));
    return x;
}
void free_d_array2d(double **x){
    free(x[0]);
    free(x);
}

double ***d_array3d(long size_1, long size_2, long size_3){
    double ***x;
    long i,j;
    x= (double ***) malloc((size_t)(size_1*sizeof(double **)));
    for(i=0; i< size_1; ++i) {
        x[i]= (double **) malloc((size_t)(size_2*sizeof(double *)));
        for(j=0; j< size_2; ++j) x[i][j]= (double *) malloc((size_t)(size_3*sizeof(double)));
    }
    return x;
}

void free_d_array3d(double ***x){
    free(x[0][0]);
    free(x[0]);
    free(x);
}


double ****d_array4d(long size_1, long size_2, long size_3, long size_4){
    double ****x;
    long i,j,k;
    x= (double ****) malloc((size_t)(size_1*sizeof(double ***)));
    for(i=0; i< size_1; ++i) {
        x[i]= (double ***) malloc((size_t)(size_2*sizeof(double **)));
        for(j=0; j< size_2; ++j){
            x[i][j]= (double **) malloc((size_t)(size_3*sizeof(double *)));
            for(k=0; k< size_3; ++k) x[i][j][k]= (double *) malloc((size_t)(size_4*sizeof(double)));
        }
    }
    return x;
}

void free_d_array4d(double ****x){
    free(x[0][0][0]);
    free(x[0][0]);
    free(x[0]);
    free(x);
}
