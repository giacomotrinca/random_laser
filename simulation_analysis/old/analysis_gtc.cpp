///////////////////////////////////////////////////////////
////////////////////GTC - 2021/////////////////////////////
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

//WELCOME////////////////////////////////////////////////
	//system("clear");
	printf("//////////////////////////////////////////////////////\n");
	printf("##SMrandomTetrads analysis: IN SPIN-GLASSES WE TRUST!#\n");
	printf("//////////////////////////////////////////////////////\n");
//VARIABLES/////////////////////////////////////////////

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
//printf("TUTTO OK\n");
	//spin s;
	int choose_mod; //choose if re-analyze the sample
	// 0: analyze the sample after simulation
	// 1 re-analyze the sample
//SYNTAX CONTROL//////////////////////////////////////

	checkSyntax(argc, 3);
	//printf("TUTTO OK\n");
	choose_mod = atoi(argv[1]);
//////////////////////////////////////////////////////

int freq_mode = atoi(argv[2]);
FILE* file_input;
if(choose_mod == 0) {
	file_input= fopen("input_analysis.dat", "r");
}else if(choose_mod == 1){
	file_input = fopen("stuffs/input_analysis.dat", "r");
}else{
	printf("WRONG MOD! EXIT.\n");
	exit(BAD_SIMULATION);
}
char* input_file_name = (char*)calloc(STRING_SIZE, sizeof(char));
checkString(input_file_name);
sprintf(input_file_name,"input_analysis.dat");
checkFile(file_input, input_file_name);
fscanf(file_input, "%d\t%lf\t%lf\t%d\t%d\t%d\t%d\t%d\t%le\t%d", &N, &T_min, &T_max, &nJUMPS, &dITER, &IT_MIN, &IT_MAX, &nRep, &Gamma, &pRINTALLCONFIG);
fclose (file_input);
/////////////////////////////////////////////////////

	int nPT = nJUMPS+1;
	int nFILE = (int)((IT_MAX-IT_MIN)/dITER) + 1; //number of configuration files
	double dT = (double)(T_max - T_min)/nPT;   //TEMPERATURE DISTANCE
	temp = (double*)calloc(nPT,sizeof(double));

	checkDouble(temp);

	for(int i = 0; i<nPT; i++) {

		temp[nPT-i-1] = T_max - (double)i*dT;
	}

	double* freq = (double*)calloc(N, sizeof(double));
	checkDouble(freq);
	char* temp_freq_filename = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(temp_freq_filename);

	FILE* freq_file;
	if(choose_mod == 0) {
		sprintf(temp_freq_filename, "frequencies.dat");
	}else if(choose_mod == 1) {
		sprintf(temp_freq_filename, "stuffs/frequencies.dat");
	}
	freq_file = fopen(temp_freq_filename, "r");

	checkFile(freq_file, temp_freq_filename);
	free(temp_freq_filename);
	int index_temp;
	double gain_temp;
	for(int i = 0; i < N; i++) {
		if(freq_mode == 0) {
			fscanf(freq_file, "%d\t%le\t%le\n", &index_temp, &freq[i], &gain_temp);
		}else if(freq_mode == 1) {
			fscanf(freq_file, "%le\n",&freq[i]);
			index_temp = 0;
			gain_temp = 0;
		}
	}

	fclose(freq_file);

	int* iter = (int*)calloc(nFILE, sizeof(int));
	checkInt(iter);

	for(int i = 0; i < nFILE; i++) {
		iter[i] = IT_MIN + i*dITER;
		//printf("%d\n", iter[i]);
	}

    clock_t start, end;
    double duration;

    start = clock();

    ///////////////////////////////////////////////////////////////////////////////
    ////////// PRELIMINARY ANALYSIS ///////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    system("mkdir graphs");
    printf("PRELIMINARY ANALYSIS BEGINS...\n");
    FILE** p = (FILE**)calloc(nRep, sizeof(FILE*));

    char dec;
    printf("I'M ACQUIRING THE DATA OF ALL THE DYNAMICS FOR EACH REPLICA...\n");
    double**** wholeDynamics = (double****)calloc(nRep, sizeof(double***));
    checkHyp2Mat(wholeDynamics);


    printf("I'M CHECKING WHETHER THE SYSTEM HAS REACHED EQUILIBRIUM..\n");
    for(int r = 0; r < nRep; r++) {
				char* temp_pipe = (char*)calloc(STRING_SIZE, sizeof(char));
				checkString(temp_pipe);
				sprintf(temp_pipe, "gnuplot");

        p[r] = popen(temp_pipe, "w");
        checkFile(p[r], temp_pipe);
				free(temp_pipe);
        //reading from file
        wholeDynamics[r] = readDynamics(r, nPT, IT_MAX, dITER, N, choose_mod);

        //transposing matrices to have temporal vectors
        wholeDynamics[r][0] = TransposeMat(wholeDynamics[r][0], IT_MAX + dITER, nPT-1);
        wholeDynamics[r][1] = TransposeMat(wholeDynamics[r][1], IT_MAX + dITER, nPT-1);

        //wholeDynamics[r][1:erg/0:acc][k][t]

        for(int k = 0; k < nPT - 1; k++) {
            meanLogaritmicWindow(wholeDynamics[r][1][k], IT_MAX+dITER, r, k, 'E');
            meanLogaritmicWindow(wholeDynamics[r][0][k], IT_MAX+dITER, r, k, 'R');
        }
        //printing energy in log scale

        TemperingCheckGnuPipe(p[r], r, nPT-1);

        //printf("Has the replica n°%d reached equilibrium?[y/n] ", r+1);

        //dec = getDecision();
        dec = 'y';
        if(dec == 'n' || dec == 'N') {
            printf("SO SAD...\n");
            exit(NOT_TERMALIZED);
        }

        pclose(p[r]);
        //printf("\n");
    }

    printf("\n");
    printf("//////////////////////////////////////////////////////\n");
    printf("PRELIMINARY ANALYSIS DONE.\n");
    printf("//////////////////////////////////////////////////////\n");
    ///////////////////////////////////////////////////////////////////////////////
    ////////// SPECIFIC HEAT //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    printf("NOW I COMPUTE SPECIFIC HEATS vs TEMPERATURES...\n");
    //ERG
    compute_specific_heat(wholeDynamics,nRep, nPT, N, dITER, IT_MAX, temp);

	free(wholeDynamics); //libero la memoria per l'intera dinamica
    printf("SPECIFIC HEATS COMPUTED. I PRINTED ALL GRAPHS \nIN THIS DIRECTORY.\n");
    printf("//////////////////////////////////////////////////////\n");

    //printf("Do you want to continue with spectrum building?[y/n] ");
    //dec = getDecision();
    if(dec == 'n' || dec == 'N') {
            printf("SO SAD...\n");
            exit(0);
        }
    printf("\n");
    ///////////////////////////////////////////////////////////////////////////////
    ////////// SPECTRUM BUILDING //////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    complex*** EqDynamics = (complex***)calloc(nRep, sizeof(complex**));
    checkHypComplexMat(EqDynamics);

    double*** A2 = (double***)calloc(nRep, sizeof(double**));
    checkHypMat(A2);
    double*** phiDyn = (double***)calloc(nRep, sizeof(double**));
    checkHypMat(phiDyn);
    printf("EQUILIBRIUM CONFIGURATION FILES READING....\n");
    //long long int tot_cycle = nRep*nFILE;
    //long long int counter = 0;
    read_equilibrium_configurations(EqDynamics, A2, phiDyn,nRep, nFILE, nPT, N,iter, choose_mod);
    complex*** AllDynamics;
    if(NITER_MIN_PRINT == 0) {
    	printf("NON EQUILIBRIUM CONFIGURATION FILES READING....\n");

	AllDynamics = (complex***)calloc(nRep, sizeof(complex**));
    	checkHypComplexMat(AllDynamics);
	read_all_configurations(AllDynamics, nRep, nFILE, nPT, N, iter, choose_mod);
    }


    printf("\nALL FILE READ, NOW I'M COMPUTING AVERAGES...\n");
		printf("//////////////////////////////////////////////////////\n");
		//PrintComplexMat(EqDynamics[0], 1, N);

		//A2[r][nFILE][N*nPT]
  		printf("Check configurations...\n");
		double**** checkA = (double****)calloc(nRep, sizeof(double***));
		checkHyp2Mat(checkA);
		double**** checkPhi = (double****)calloc(nRep, sizeof(double***));
		checkHyp2Mat(checkPhi);

		check_spherical_constrait(checkA, checkPhi, A2, phiDyn, nRep, nFILE, nPT, N);

	//checkA[replica][iter][npt][mode]

	if(PRINTALLCONFIG == 1) {
		printAllSpectra(checkA, checkPhi, freq, iter, nFILE, nPT, nRep, N);
	}
	free(checkA);
	free(checkPhi);
	free(phiDyn);
	printf("\n");


	for(int r = 0; r < nRep; r++) {
		if(NITER_MIN_PRINT == 0) {
			A2[r] = TransposeMat(A2[r], (int)(nFILE/2) +1, N*nPT);
		}else{
			A2[r] = TransposeMat(A2[r], nFILE, N*nPT);
		}
	}

	//A2[r][N*nPT][nFILE]

	double** A2M = (double**)calloc(nRep, sizeof(double*));
	checkMat(A2M);

	//compute averages
	for(int r = 0; r < nRep; r++) {
		A2M[r] = (double*)calloc(N*nPT, sizeof(double));
		checkDouble(A2M[r]);
		for(int i = 0; i < N*nPT; i++) {
			if(NITER_MIN_PRINT == 0) {
				A2M[r][i] = MeanVec(A2[r][i], (int)(nFILE/2) +1);
			}else{
				A2M[r][i] = MeanVec(A2[r][i], nFILE);
			}
		}
	}

	double** sigma = (double**)calloc(nRep, sizeof(double*));
	checkMat(sigma);

	//compute standard deviations
	for(int r = 0; r < nRep; r++) {
		sigma[r] = (double*)calloc(N*nPT, sizeof(double));
		checkDouble(sigma[r]);
		for(int i = 0; i < N*nPT; i++) {
			if(NITER_MIN_PRINT == 0) {
				sigma[r][i] = StandardDeviation(A2[r][i], A2M[r][i], (int)(nFILE/2) +1);
			}else{
				sigma[r][i] = StandardDeviation(A2[r][i], A2M[r][i], nFILE);
			}
		}
	}

	//free(A2);
	printf("AVERAGES COMPUTED, NOW I'M PRINTING INTENSITIES...\n");

	double*** intensity = (double***)calloc(nRep, sizeof(double**));
	double*** sigma_block = (double***)calloc(nRep, sizeof(double**));
	checkHypMat(intensity);
	checkHypMat(sigma_block);
	//printf("TUTTO OK\n");

	for(int r = 0; r < nRep; r++) {
		intensity[r] = (double**)calloc(nPT, sizeof(double*));
		checkMat(intensity[r]);
		sigma_block[r] = (double**)calloc(nPT, sizeof(double*));
		checkMat(sigma_block[r]);
		for(int k = 0; k < nPT; k++) {
			intensity[r][k] = DeBlock(A2M[r], N, nPT, k);
			sigma_block[r][k] = DeBlock(sigma[r], N, nPT, k);
			Deploy(intensity[r][k], sigma_block[r][k], N, temp[k], r,  k, choose_mod, freq_mode);

		}
		PlotSpectrumScript(r, nPT, N);
	}

	free(intensity);
	free(sigma_block);
	printf("INTENSITIES COMPUTED. I PRINTED ALL GRAPHS \nIN THIS DIRECTORY.\n");
	printf("//////////////////////////////////////////////////////\n");
	printf("//////////////////////////////////////////////////////\n");



	if(nRep > 1) {
		//printf("Do you continue with overlaps?[y/n] ");
    //dec = getDecision();
    		if(dec == 'n' || dec == 'N') {
            		printf("SO SAD...\n");
            		exit(0);
      		}
   		 printf("\n");


		printf("NOW I'M COMPUTING OVERLAPS\n");
	 	printf("//////////////////////////////////////////////////////\n");
		printf("NOW IT'S TIME FOR IFO...\n");
		printf("//////////////////////////////////////////////////////\n");


   		printf("Deblocking Spectrum....\n");
    		for(int r = 0; r < nRep; r++) {
        		A2[r] = TransposeMat(A2[r], N*nPT, nFILE);
    		}
    
    		complex**** FixedDynamics;
		complex**** FixedAllDynamics;
    
    		if(NITER_MIN_PRINT == 0) {
			FixedDynamics = fix_dynamics(EqDynamics, nRep, nPT,  N,(int)(nFILE/2) +1);
			printf("Equilibrium...");
			FixedAllDynamics = fix_dynamics(AllDynamics, nRep, nPT,  N, nFILE);
			printf("Non Equilibrium...");
			free(AllDynamics);
		}else{
			FixedDynamics = fix_dynamics(EqDynamics, nRep, nPT,  N, nFILE);
		}
		free(EqDynamics);
		//DEBLOCKING INTENSITIES
		//counter = 0;
    		//tot_cycle = nPT*nFILE*nRep;
		//A2[r][nFILE][N*nPT]
		//A2M[r][N*nPT]
		double** cVec = compute_intensity_fluctuations_overlaps(A2,A2M, nRep, nFILE, nPT, N);
		free(A2M);
		free(A2);
		for(int k = 0; k < nPT; k++) {
       			if(NITER_MIN_PRINT == 0 ) {
				HistogramIFO(cVec[k],(int)(nRep*(nRep-1)/2)*((int)(nFILE/2) +1), k, temp);
			}else{
				HistogramIFO(cVec[k],(int)(nRep*(nRep-1)/2)*nFILE, k, temp);
			}
        	}
        	free(cVec);
		
		printf("Experimental IFOs...\n");
		double** cVec_experimental = experimental_equilibrium_ifo(FixedDynamics, nRep, (int)(nFILE/2) +1, nPT, N);
		for(int k = 0; k < nPT; k++) {
			HistogramExpIFO(cVec_experimental[k], (int)(nRep*(nRep-1)/2), k, temp, -1, 0);
		}
		free(cVec_experimental);

		if(NITER_MIN_PRINT == 0) {
			dynamics_experimental_ifo(FixedAllDynamics, nRep, nFILE, nPT, N,IT_MIN, dITER, temp, 0);
			dynamics_experimental_ifo(FixedAllDynamics, nRep, nFILE, nPT, N,IT_MIN, dITER, temp, 1);

		}

	
    printf("HISTOGRAMS PRINTED.\n");
	
	complex**** q;

	if(NITER_MIN_PRINT == 0) {
		q = compute_parisi_overlaps(FixedAllDynamics,nRep,nFILE,nPT,N);
	    free(FixedAllDynamics);
	}else{
		q = compute_parisi_overlaps(FixedDynamics,nRep,nFILE,nPT,N);
	}
	
	//free(FixedDynamics);
	printf("\n");
	//PrintComplexMat(q[0][0], nRep, nRep);
	printf("//////////////////////////////////////////////////////\n");
	printf("NOW I COMPUTE P(q) FOR THIS SAMPLE...\n");
	printf("//////////////////////////////////////////////////////\n");

	if(NITER_MIN_PRINT == 0) {
		//andiamo a campionare le P(q) per finestre lunghe 2^l indipendenti.

		int b_max = (int)log2((double)nFILE);
		//printf("%d\t%d\n", iter[(int)(nFILE/2)], iter[(int)(nFILE/2)+(int)(nFILE/2)-1]);
		//printf("FACCIO CONTI....\n");
		int it_min_temp = IT_MIN;
		int it_max_temp = 0;
		int nfile_temp = 0;
		for(int b = 0; b < b_max; b++) {
			if(b == 0) {
				it_max_temp = dITER;
			}else{
				it_max_temp = 2*it_min_temp - dITER;
			}
			//leggere i file da it_min_temp a it_max_temp
			//readConfigFILE(r, iter[i], N*nPT);
			nfile_temp = (int)((it_max_temp - it_min_temp)/dITER) + 1;

			//q[k][i][r1][r2]

			complex**** q_temp = (complex****)calloc(nPT, sizeof(complex***));
			checkHyp2ComplexMat(q_temp);

			q_temp = overlaps_dynamics(q, nPT, nfile_temp,nRep,b);
			
			double** qMatRe_temp = compute_equilibrium_pq_real_part(q_temp, nRep, nPT, nfile_temp);
			double** qMatIm_temp = (double**)calloc(nPT, sizeof(double*));
			checkMat(qMatIm_temp);
			if(IM_PART_PQ == 1) {
				qMatIm_temp = compute_equilibrium_pq_imaginary_part(q_temp, nRep, nPT, nfile_temp);
			}
			free(q_temp);

			for(int k = 0; k < nPT; k++) {
	        	HistogramQ(qMatRe_temp[k],(int)(nRep*(nRep-1)/2)*nfile_temp, k, 0, temp, b);
				if(IM_PART_PQ == 1) {
					HistogramQ(qMatIm_temp[k],nRep*(nRep-1)*nfile_temp, k, 1, temp, b);
				}
	      //counter++;
	      //printf("%.1lf%%\r", ((double)counter/tot_cycle)*100);
	    	}

			printf("%d\t%d\t%d\n", it_min_temp, it_max_temp, nfile_temp);
			free(qMatIm_temp);
			free(qMatRe_temp);
			it_min_temp = it_max_temp + dITER;

		}
		free(q);

	}else{

  	printf("Preparing matrix for histograms...\n");
		double** qMatRe = compute_equilibrium_pq_real_part(q, nRep, nPT, nFILE);
		double** qMatImFixed;
		if(IM_PART_PQ == 1) {
			qMatImFixed = compute_equilibrium_pq_imaginary_part(q, nRep, nPT, nFILE);
		}
		free(q);

    for(int k = 0; k < nPT; k++) {
      HistogramQ(qMatRe[k],(int)(nRep*(nRep-1)/2)*nFILE, k, 0, temp, -1);
			if(IM_PART_PQ == 1) {
				HistogramQ(qMatImFixed[k],nRep*(nRep-1)*nFILE, k, 1, temp, -1);
			}
      //counter++;
      //printf("%.1lf%%\r", ((double)counter/tot_cycle)*100);
    }
    	printf("Done histograms for P(q).\n");

		
    printf("HISTOGRAMS PRINTED.\n");
    printf("//////////////////////////////////////////////////////\n");


	} //FINE ELSE IF


    printf("\n");
    //printf("Do you want to continue with IFOs?[y/n] ");
    //dec = getDecision();
    if(dec == 'n' || dec == 'N') {
            printf("SO SAD...\n");
            exit(0);
    }
    printf("\n");

    ///////////////////////////////////////////////////////////////////////////////////////////////
 
    printf("//////////////////////////////////////////////////////\n");

     //printf("Do you want to continue with plaquettes overlaps?[y/n] ");
    //dec = getDecision();
    if(dec == 'n' || dec == 'N') {
            printf("SO SAD...\n");
            exit(0);
    }
    printf("\n");

    printf("//////////////////////////////////////////////////////\n");
    printf("PLAQUETTE OVERLAPS\n");
    printf("//////////////////////////////////////////////////////\n");


    //CALCOLO IL NUMERO DI PLACCHETTE DOPO LA DILUZIONE DOVUTA ALLA FMC
    printf("Computing number of tetrads...\n");


    int Nplaq = NumTetrads(freq, N, Gamma);
		Nplaq = NearPowTwo(Nplaq);
    printf("Reduced number of tetrads: %d\n", Nplaq);
    printf("Reading interactions_file.dat\n");

    double** intData = ReadInteractionFile(Nplaq, choose_mod);


    printf("Computing epsilons...\n");

    //counter = 0;
    //tot_cycle = nPT*nRep*(nRep-1)*nFILE/2;

    double** eMat = compute_plaquette_overlaps(FixedDynamics, intData,nRep,nPT,Nplaq,nFILE);
	free(FixedDynamics);
	free(intData);
    printf("\n");
    printf("Printing histograms....\n");

    for(int k = 0; k < nPT; k++) {
      if(NITER_MIN_PRINT == 0) {
				HistogramEPS(eMat[k],(int)(nRep*(nRep-1)/2)*(int)(nFILE/2) +1, k, temp);
			}else{
				HistogramEPS(eMat[k],(int)(nRep*(nRep-1)/2)*nFILE, k, temp);
			}
    }
	} else { //FINE IF SUL NUMERO DI REPLICHE
		free(A2M);
		free(A2);
	}
	end = clock();
	duration = ((double)(end-start))/CLOCKS_PER_SEC;
	printf("\nTotal time:\t%.1lf s\n", duration);

	OrganizeFiles(PRINTALLCONFIG, nRep);

	free(temp);
	free(iter);

	return 0;
}
