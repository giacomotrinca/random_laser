///////////////////////////////////////////////////////////
////////////////////GTC - 2022/////////////////////////////
///////////////////////////////////////////////////////////

//////////// STANDARD LIBRARIES////////////////////////////
#include<stdlib.h>
#include<stdio.h>
#include<ctime>
#include<string.h>
#include<iostream>

//////////// MATH & STATISTICS LIBRARIES //////////////////
#include<math.h>

//////////// CUSTOM LIBRARIES & STRUCTURES/////////////////
#include "gtcStructures.h"
#include "SMrandomTetradsRUNCHECK.h"
#include "SMrandomTetradsSettings.h"
#include "SMrandomTetrads_structures.h"

#include "gtcTools.h"

using namespace std;


//////////// MAIN ////////////////////////////////////////

int main(int argc, char** argv) {


	welcomeScreen();
	
	//INIT VARIABLES;
	
	int N;				//SIZE
	double T_min;			//T_MIN
	double T_max;			//T_MAX
	int nJUMPS;			//JUMPS OF TEMPERATURES
	int dITER; 			//INPUT FILE ITERATION DISTANCE
	int IT_MIN;			//FIRST ITERATION BEFORE TERMALIZATION
	int IT_MAX;			//LAST ITERATION PRINTED
	int nRep;			//REPLICAS NUMBER (2,4)
	double* temp;			//temperature array
	double Gamma;
	int pRINTALLCONFIG;
	
	clock_t start = clock();
	char* input_file_name = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(input_file_name);
	sprintf(input_file_name,"input_analysis.dat");
	FILE* file_input = fopen(input_file_name, "r");
	checkFile(file_input, input_file_name);
	free(input_file_name);
	fscanf(file_input, "%d\t%lf\t%lf\t%d\t%d\t%d\t%d\t%d\t%le\t%d", &N, &T_min, &T_max, &nJUMPS, &dITER, &IT_MIN, &IT_MAX, &nRep, &Gamma, &pRINTALLCONFIG);
	fclose (file_input);
	
	//FURTHER VARIABLES
	int nFILE = (int)((IT_MAX-IT_MIN)/dITER) + 1; //number of configuration files
	double dT = (double)(T_max - T_min)/NPT;   //TEMPERATURE DISTANCE
	temp = (double*)calloc(NPT,sizeof(double));
	checkDouble(temp);

	for(int i = 0; i<NPT; i++) {
		temp[NPT-i-1] = T_max - (double)i*dT;
	}

	
	double* freq = getFrequencies(N);
	
	int* iter = (int*)calloc(nFILE, sizeof(int));
	checkInt(iter);
	for(int i = 0; i < nFILE; i++) {
		iter[i] = IT_MIN + i*dITER;
	}
	
	printf("Reading Parallel_Tempering Files...\n");
	double**** wholeDynamics = (double****)calloc(nRep, sizeof(double***));
    	checkHyp2Mat(wholeDynamics);
    	
    	
    	for(int r = 0; r < NREPLICAS; r++) {
    		char* temp_pipe = (char*)calloc(STRING_SIZE, sizeof(char));
		checkString(temp_pipe);
		sprintf(temp_pipe, "gnuplot");
    		FILE* p = popen(temp_pipe, "w");
        	checkFile(p, temp_pipe);
		free(temp_pipe);
		wholeDynamics[r] = readDynamics(r, NPT, IT_MAX, dITER, N, 0);
		wholeDynamics[r][0] = TransposeMat(wholeDynamics[r][0], IT_MAX + dITER, NPT-1);
       		wholeDynamics[r][1] = TransposeMat(wholeDynamics[r][1], IT_MAX + dITER, NPT-1);
       		
       		for(int k = 0; k < NPT - 1; k++) {
            		meanLogaritmicWindow(wholeDynamics[r][1][k], IT_MAX+dITER, r, k, 'E');
            		meanLogaritmicWindow(wholeDynamics[r][0][k], IT_MAX+dITER, r, k, 'R');
        	}//END TEMPERATURE FOR
        	
        	TemperingCheckGnuPipe(p, r, NPT-1);
        	pclose(p);
    	}//END REPLICA FOR
	
	
	printf("//////////////////////////////////////////////////////\n");
	printf("Computing Specific Heats...\n");
   	printf("//////////////////////////////////////////////////////\n");
	compute_specific_heat(wholeDynamics,nRep, NPT, N, dITER, IT_MAX, temp);
	printf("Free Dynamics memory..\n");
	
	for(int r = 0; r < NREPLICAS; r++) {
		for(int l = 0; l < 2; l++) {
			for(int k = 0; k < NPT - 1; k++) {
				free(wholeDynamics[r][l][k]);
			}
			free(wholeDynamics[r][l]);
		}
		free(wholeDynamics[r]);
	}//wholeDynamics[r][1:erg/0:acc][k][t]
	
	free(wholeDynamics);
	
	printf("//////////////////////////////////////////////////////\n");
	printf("Reading Configurations...\n");
	printf("//////////////////////////////////////////////////////\n");
	
	complex**** conf = (complex****)calloc(NREPLICAS, sizeof(complex***));
	checkHyp2ComplexMat(conf);
	
	for(int r = 0; r < NREPLICAS; r++) {
		conf[r] = (complex***)calloc(nFILE, sizeof(complex**));
		checkHypComplexMat(conf[r]);
		for(int i = 0; i < nFILE; i++) {
			conf[r][i] = readConfigFILE_v2(r, iter[i], N);	
		}
	}
	
	printf("Configuration files read...\n");
	//conf[r][i][k][j]
	printf("Allowing configuration indices...\n");
	complex**** new_conf = allow_conf(conf, nFILE, N);
	printf("Free old conf memory\n");
	for(int r = 0; r < NREPLICAS; r++) {
		for(int i = 0; i < nFILE; i++) {
			for(int k = 0; k < NPT; k++) {
				free(conf[r][i][k]);
			}
			free(conf[r][i]);
		}
		free(conf[r]);
	}
	free(conf);
	//new_conf[r][k][j][i]
	printf("Computing intensities...\n");
	double**** intensity = (double****)calloc(NREPLICAS, sizeof(double***));
	checkHyp2Mat(intensity);
	for(int r = 0; r < NREPLICAS; r++) {
		intensity[r] = (double***)calloc(NPT, sizeof(double**));
		checkHypMat(intensity[r]);
		for(int k = 0; k < NPT; k++) {
			intensity[r][k] = (double**)calloc(N, sizeof(double*));
			checkMat(intensity[r][k]);
			for(int j = 0; j < N; j++) {
				intensity[r][k][j] = (double*)calloc(nFILE, sizeof(double));
				checkDouble(intensity[r][k][j]);
				for(int i = 0; i < nFILE; i++) {
					intensity[r][k][j][i] = new_conf[r][k][j][i].re * new_conf[r][k][j][i].re + new_conf[r][k][j][i].im * new_conf[r][k][j][i].im;
				}
			}
		}
	}
	printf("Intensity computed.\n");
	printf("Building spectrum...\n");
	printf("Computing mean intensities...\n");
	
	int i_min;
	if(NITER_MIN_PRINT == 0) {
		i_min = (int)(nFILE/2) + 1;	
	}else{
		i_min = 0;
	}
	double**** eq_intensity = (double****)calloc(NREPLICAS, sizeof(double***));
	checkHyp2Mat(eq_intensity);
	double*** eq_mean_intensity = (double***)calloc(NREPLICAS, sizeof(double**));
	checkHypMat(eq_mean_intensity);
	double*** eq_sigma_intensity = (double***)calloc(NREPLICAS, sizeof(double**));
	checkHypMat(eq_sigma_intensity);

	for(int r = 0; r < NREPLICAS; r++) {
		eq_mean_intensity[r] = (double**)calloc(NPT, sizeof(double*));
		checkMat(eq_mean_intensity[r]);
		eq_intensity[r] = (double***)calloc(NPT, sizeof(double**));
		checkHypMat(eq_intensity[r]);
		eq_sigma_intensity[r] = (double**)calloc(NPT, sizeof(double*));
		checkMat(eq_sigma_intensity[r]);

		for(int k = 0; k < NPT; k++) {
			eq_mean_intensity[r][k] = (double*)calloc(N, sizeof(double));
			checkDouble(eq_mean_intensity[r][k]);
			eq_intensity[r][k] = (double**)calloc(N, sizeof(double*));
			checkMat(eq_intensity[r][k]);
			eq_sigma_intensity[r][k] = (double*)calloc(N, sizeof(double));
			checkDouble(eq_sigma_intensity[r][k]);
			for(int j = 0; j < N; j++) {
				eq_intensity[r][k][j] = (double*)calloc(nFILE - i_min, sizeof(double));
				for(int i = 0; i < nFILE - i_min; i++) {
					eq_intensity[r][k][j][i] = intensity[r][k][j][i+i_min];
					eq_mean_intensity[r][k][j] += eq_intensity[r][k][j][i] / ((nFILE - i_min)*sqrt(temp[k]));
				}
				eq_sigma_intensity[r][k][j] = StandardDeviation(eq_intensity[r][k][j], eq_mean_intensity[r][k][j], nFILE - i_min)/sqrt(temp[k]);
			}
			build_spectrum(eq_mean_intensity[r][k], eq_sigma_intensity[r][k], N, freq, r, k);
		}
		PlotSpectrumScript(r, NPT, N);
	}
	
	for(int r = 0; r < NREPLICAS; r++) {
		for(int k = 0; k < NPT; k++) {
			for(int j = 0; j < N; j++) {
				eq_mean_intensity[r][k][j] *= sqrt(temp[k]);
			}
		}
	}
	
	//new_conf[r][k][j][i]
	
	printf("//////////////////////////////////////////////////////\n");
	printf("Computing equilibrium theoretical Ifos...\n");
	printf("//////////////////////////////////////////////////////\n");
	
	//intensity[r][k][j][i]
	double** t_ifo = (double**)calloc(NPT, sizeof(double*));
	checkMat(t_ifo);
	for(int r = 0; r < NREPLICAS; r++) {
		for(int k = 0; k < NPT; k++) {
			eq_intensity[r][k] = TransposeMat(eq_intensity[r][k], N, nFILE-i_min);
		}
	}
	
	for(int k = 0; k < NPT; k++) {
		t_ifo[k] = (double*)calloc((nFILE - i_min)*(int)(NREPLICAS*(NREPLICAS - 1)/2), sizeof(double));
		checkDouble(t_ifo[k]);
		int l = 0;
		char* name_ifo = (char*)calloc(STRING_SIZE, sizeof(char));
    		checkString(name_ifo);
    		sprintf(name_ifo, "ifo_overlaps_temp%d.dat", k);

    		FILE* file_ifo = fopen(name_ifo, "w");
    		checkFile(file_ifo, name_ifo);
    		free(name_ifo);

    		fprintf(file_ifo, "#n_iter\tr1\tr2\tC\n");
		for(int i = 0; i < nFILE - i_min; i++) {
			for(int r1 = 0; r1 < NREPLICAS; r1++) {
				for(int r2 = 0; r2 < r1; r2++) {
					t_ifo[k][l] = OverlapIFO_Element(eq_intensity[r1][k][i], eq_intensity[r2][k][i], eq_mean_intensity[r1][k], eq_mean_intensity[r2][k], N);
					fprintf(file_ifo, "%d\t%d\t%d\t%le\n", i, r1, r2, t_ifo[k][l]);
					l++;
				}
			}
		}
		fclose(file_ifo);
	}
	
	
	for(int k = 0; k < NPT; k++) {
		HistogramIFO(t_ifo[k],(nFILE - i_min)*(int)(NREPLICAS*(NREPLICAS - 1)/2), k, temp);
		free(t_ifo[k]);
	}
	free(t_ifo);
	
	
	printf("Experimental Equilibrium Ifos...\n");
	double** e_ifo = (double**)calloc(NPT, sizeof(double*));
	double** eq_mean_rep_intensity = (double**)calloc(NPT, sizeof(double*));
	
	for(int k = 0; k < NPT; k++) {
		eq_mean_rep_intensity[k] = (double*)calloc(N, sizeof(double));
		checkDouble(eq_mean_rep_intensity[k]);
		for(int j = 0; j < N; j++) {
			for(int r = 0; r < NREPLICAS; r++) {
				eq_mean_rep_intensity[k][j] += eq_mean_intensity[r][k][j]/NREPLICAS;
			}
		}
	}
	for(int k = 0; k < NPT; k++) {
		int l = 0;
  		char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
  		checkString(filename);
  		sprintf(filename, "experimental_ifo_eq_temp%d.dat", k);
  		FILE* fw = fopen(filename, "w");
  		checkFile(fw, filename);
  		free(filename);
  		e_ifo[k] = (double*)calloc((int)(NREPLICAS*(NREPLICAS-1)/2), sizeof(double));
  		checkDouble(e_ifo[k]);
		for(int r1 = 0; r1 < NREPLICAS; r1++) {
			for(int r2 = 0; r2 < r1; r2++) {
				e_ifo[k][l] = OverlapIFO_Element(eq_mean_intensity[r1][k], eq_mean_intensity[r2][k], eq_mean_rep_intensity[k], eq_mean_rep_intensity[k], N);
        			fprintf(fw, "%d\t%d\t%le\n", r1, r2, e_ifo[k][l]);
        			l++;
			}
		}
		fclose(fw);
	}
	for(int k = 0; k < NPT; k++) {
		HistogramExpIFO(e_ifo[k], (int)(nRep*(nRep-1)/2), k, temp, -1, 0);
		free(e_ifo[k]);
	}
	free(e_ifo);
	for(int r = 0; r < NREPLICAS; r++) {
		for(int k = 0; k < NPT; k++) {
			eq_intensity[r][k] = TransposeMat(eq_intensity[r][k], nFILE-i_min, N);
		}
	}
	printf("Free equilibrium intensity memory...\n");
	for(int k = 0; k < NPT; k++) {
		free(eq_mean_rep_intensity[k]);
	}
	free(eq_mean_rep_intensity);
	for(int r = 0; r < NREPLICAS; r++) {
		for(int k = 0; k < NPT; k++) {
			free(eq_sigma_intensity[r][k]);
			free(eq_mean_intensity[r][k]);
			for(int j = 0; j < N; j++) {
				free(eq_intensity[r][k][j]);
			}
			free(eq_intensity[r][k]);
		}
		free(eq_sigma_intensity[r]);
		free(eq_mean_intensity[r]);
		free(eq_intensity[r]);
	}
	free(eq_sigma_intensity);
	free(eq_mean_intensity);
	free(eq_intensity);
	
	//intensity[r][k][j][i]
	
	if(NITER_MIN_PRINT == 0) {
		printf("Computing Experimental Dynamical Ifos...\n");
		int it_min_temp = IT_MIN;
  		int it_max_temp = 0;
  		int nfile_temp = 0;
		int bmax = (int)log2((double)nFILE);
		
		for(int b = 0; b < bmax; b++) {
			
			if(b == 0) {
      				it_max_temp = dITER;
    			}else{
     				it_max_temp = 2*it_min_temp - dITER;
   			}
   			//INDEPENDENT BLOCK MEAN
   			nfile_temp = (int)((it_max_temp - it_min_temp)/dITER) + 1;
   			double*** mean_ind_block_int = (double***)calloc(NREPLICAS, sizeof(double**));
   			checkHypMat(mean_ind_block_int);

   			for(int r = 0; r < NREPLICAS; r++) {
   				mean_ind_block_int[r] = (double**)calloc(NPT, sizeof(double*));
   				checkMat(mean_ind_block_int[r]);
   				for(int k = 0; k < NPT; k++) {
   					mean_ind_block_int[r][k] = (double*)calloc(N, sizeof(double));
   					checkDouble(mean_ind_block_int[r][k]);
   					for(int j = 0; j < N; j++) {
   						
   						for(int i = (int)(it_min_temp/dITER); i < (int)(it_max_temp/dITER)+1; i++) {
   							mean_ind_block_int[r][k][j] += intensity[r][k][j][i]/nfile_temp;
   						}
   					}
   				}
   			}
   			//CUMULATIVE BLOCK MEAN
   			nfile_temp = (int)(it_max_temp/dITER) + 1;
   			
   			double*** mean_cum_block_int = (double***)calloc(NREPLICAS, sizeof(double**));
   			checkHypMat(mean_cum_block_int);
   			for(int r = 0; r < NREPLICAS; r++) {
   				mean_cum_block_int[r] = (double**)calloc(NPT, sizeof(double*));
   				checkMat(mean_cum_block_int[r]);
   				for(int k = 0; k < NPT; k++) {
   					mean_cum_block_int[r][k] = (double*)calloc(N, sizeof(double));
   					checkDouble(mean_cum_block_int[r][k]);
   					for(int j = 0; j < N; j++) {
   						for(int i = 0; i < (int)(it_max_temp/dITER)+1; i++) {
   							mean_cum_block_int[r][k][j] += intensity[r][k][j][i]/nfile_temp;
   						}
   						
   					}
   					
   				}
   				
   			}
   			
   			
   			
   			//REPLICAS MEANS
   			double** mean_rep_ind_int = (double**)calloc(NPT, sizeof(double*));
   			checkMat(mean_rep_ind_int);
   			double** mean_rep_cum_int = (double**)calloc(NPT, sizeof(double*));
   			checkMat(mean_rep_cum_int);
   			
   			for(int k = 0; k < NPT; k++) {
   				mean_rep_ind_int[k] = (double*)calloc(N, sizeof(double));
   				checkDouble(mean_rep_ind_int[k]);
   				mean_rep_cum_int[k] = (double*)calloc(N, sizeof(double));
   				checkDouble(mean_rep_cum_int[k]);
   				for(int j = 0; j < N; j++) {
   					for(int r = 0; r < NREPLICAS; r++) {
   						mean_rep_ind_int[k][j] += mean_ind_block_int[r][k][j]/NREPLICAS;
   						mean_rep_cum_int[k][j] += mean_cum_block_int[r][k][j]/NREPLICAS;
   					}
   				}
   			}
   			
   			
   			
   			double** ind_ifo = (double**)calloc(NPT, sizeof(double*));
   			checkMat(ind_ifo);
   			double** cum_ifo = (double**)calloc(NPT, sizeof(double*));
   			checkMat(cum_ifo);
   			for(int k = 0; k < NPT; k++) {	
   				FILE* f_ind;
   				FILE* f_cum;
   				char* filename_i = (char*)calloc(STRING_SIZE, sizeof(char));
   				checkString(filename_i);
   				char* filename_c = (char*)calloc(STRING_SIZE, sizeof(char));
   				checkString(filename_c);
   				sprintf(filename_i, "ifo_overlaps_block_%d_temp%d.dat", b, k);
   				sprintf(filename_c, "ifo_overlaps_cumulative_block_%d_temp%d.dat", b, k);
   				f_ind = fopen(filename_i, "w");
   				f_cum = fopen(filename_c, "w");
   				checkFile(f_ind, filename_i);
   				checkFile(f_cum, filename_c);
   				free(filename_i);
   				free(filename_c);
   				fprintf(f_ind, "#r1\tr2\tC\n");
   				fprintf(f_cum, "#r1\tr2\tC\n");
   				
   				ind_ifo[k] = (double*)calloc((int)(NREPLICAS*(NREPLICAS-1)/2), sizeof(double));
   				checkDouble(ind_ifo[k]);
   				cum_ifo[k] = (double*)calloc((int)(NREPLICAS*(NREPLICAS-1)/2), sizeof(double));
   				checkDouble(cum_ifo[k]);
   				int l = 0;
   				for(int r1 = 0; r1 < NREPLICAS; r1++) {
   					for(int r2 = 0; r2 < r1; r2++) {
   						ind_ifo[k][l] =  OverlapIFO_Element(mean_ind_block_int[r1][k], mean_ind_block_int[r2][k], mean_rep_ind_int[k], mean_rep_ind_int[k], N);
   						cum_ifo[k][l] =  OverlapIFO_Element(mean_cum_block_int[r1][k], mean_cum_block_int[r2][k], mean_rep_cum_int[k], mean_rep_cum_int[k], N);
   						fprintf(f_ind, "%d\t%d\t%le\n", r1, r2, ind_ifo[k][l]);
   						fprintf(f_cum, "%d\t%d\t%le\n", r1, r2, cum_ifo[k][l]);
   						l++;
   					}
   				}
   				fclose(f_ind);
   				fclose(f_cum);
   				HistogramExpIFO(ind_ifo[k], (int)(NREPLICAS*(NREPLICAS-1)/2), k, temp, b, 0);
   				HistogramExpIFO(cum_ifo[k], (int)(NREPLICAS*(NREPLICAS-1)/2), k, temp, b, 1);
   				
   				
   			}

   			for(int k = 0; k < NPT; k++) {
   				free(ind_ifo[k]);
   				free(cum_ifo[k]);
   				free(mean_rep_ind_int[k]);
   				free(mean_rep_cum_int[k]);
   			}
   			free(ind_ifo);
   			free(cum_ifo);
   			free(mean_rep_ind_int);
   			free(mean_rep_cum_int);
   			
   			for(int r = 0; r < NREPLICAS; r++) {
   				for(int k = 0; k < NPT; k++) {
   					free(mean_ind_block_int[r][k]);
   					free(mean_cum_block_int[r][k]);
   				}
   				free(mean_ind_block_int[r]);
   				free(mean_cum_block_int[r]);
   			}
   			free(mean_ind_block_int);
   			free(mean_cum_block_int);
   			
   			it_min_temp = it_max_temp + dITER;
		}
	
	}
	
	for(int r = 0; r < NREPLICAS; r++) {
		for(int k = 0; k < NPT; k++) {
			for(int j = 0; j < N; j++) {
				free(intensity[r][k][j]);
			}
			free(intensity[r][k]);
		}
		free(intensity[r]);
	}
	free(intensity);
	
	printf("//////////////////////////////////////////////////////\n");
	printf("Computing Parisi Overlaps...\n");
	printf("//////////////////////////////////////////////////////\n");
	
	
	
	OrganizeFiles(pRINTALLCONFIG, nRep);
	clock_t end = clock();
	
	double duration = (double)(end-start)/CLOCKS_PER_SEC;
	
	if(duration > 3600) {
		printf("time:\t%lfh\n", duration/3600);
	}else if(duration > 60) {
		printf("time:\t%lfm\n", duration/60);
	}else{
		printf("time:\t%.0lfs\n", duration);
	}
	
	
	
	return 0;
}
