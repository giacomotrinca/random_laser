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

	
	//Stavolta le temperature verranno immesse dalla riga di comando
	
	checkSyntax(3, argc);
	
	welcomeScreen(3.1);
	
	//Cercherò di fare a meno del file di input. Quindi, utilizzerò i parametri del file di configurazione delle simulazione
	
	int dITER = NITER_PRINT_CONF;
	int N = Size;
	int IT_MIN = NITER_MIN_PRINT;
	int IT_MAX = NITERATIONS-NITER_PRINT_CONF;
	int nFILE = (int)((IT_MAX-IT_MIN)/dITER) + 1; //number of configuration files
	double t_min = atof(argv[1]);path = data.choose_path(N, eq)
	double t_max = atof(argv[2]);
	double t_spacing = (double)(t_max-t_min)/NPT;
	
	
	double* temp = (double*)calloc(NPT, sizeof(double));
	checkDouble(temp);
	
	for(int i = 0; i < NPT; i++) {
	//Le temperature sono adesso in ordine crescente, più comodo no?
		temp[NPT - i- 1] = t_max - (double)i * t_spacing;
	}
	
	//Adesso acquisiamo le frequenze.
	
	double* freq = getFrequencies(N);
	
	//Mi serve un vettore che tenga in memoria le iterazioni, così poi da leggere il giusto file di configurazione
	int* iter = (int*)calloc(nFILE, sizeof(int));
	checkInt(iter);
	for(int i = 0; i < nFILE; i++) {
		iter[i] = IT_MIN + i*dITER;
	}
	
	//Siamo pronti ad acquisire i file parallel_tempering
	printf("Reading Parallel_Tempering Files...\n");
	double*** erg = read_parallel_tempering(IT_MAX + dITER);
	
	//Possiamo calcolare ora il calore specifico sulla metà dei dati che abbiamo acquisito 
	printf("Computing specific heat...\n");
	compute_specific_heat_V3(erg, IT_MAX+dITER, temp);
	
	//erg non serve più: lo libero.
	freeMemory3(erg, NREPLICAS, NPT-1);
	
	
	//acquisisco le configurazioni
	printf("Configuration files reading...\n");
	
	//da qui il programma prende due strade, a seconda se studiamo solo le proprietà di equilibrio (PT) o le proprietà dinamiche (NOPT)
	
	if(PT_EXCHANGE == 1) {
		//NEL CASO IN CUI VI SIA IL PT i dati che abbiamo stampato sono tutti di equilibrio	
		//spin[r][k][j][i]
		complex**** spin = load_configurations_PT(iter, nFILE);
		printf("Configurations loaded!\n");
		
		//spettro per ogni passo all'equilibrio
		printf("Computing spectrum....\n");
		double**** spectrum = instant_spectrum_PT(spin, nFILE);	
		printf("Instant spectrum computed!\n");
		printf("Computing mean spectrum....\n");
		double*** intensity = mean_spectrum_PT(spectrum, nFILE, freq, temp); //media dell'intensità sui passi all'equilibrio
		printf("Computing theoretical IFOs...\n");
		compute_theo_ifo_PT(spectrum, intensity, nFILE, temp);
		//theo_ifo[k][r1][r2]
		printf("Computing experimental IFOs...\n");
		compute_exp_ifo_PT(intensity, nFILE, temp);
		freeMemory4(spectrum, NREPLICAS, NPT, Size);
		freeMemory3(intensity, NREPLICAS, NPT);
		printf("Computing Parisi overlaps...\n");
		compute_parisi_overlap_PT(spin, nFILE, temp);
		freeMemory4C(spin, NREPLICAS, NPT, Size);
		
		
	}else{
		//SENZA PT stampiamo tutti i dati, quindi li carichiamo in blocchi di lunghezza esponenziale (il blocco b è lungo 2^b)
		//spin[r][k][j][b][i]
		int b_max = (int)log2(nFILE);
		complex***** spin = load_configurations_NOPT(iter, nFILE);
		printf("Configurations loaded!\n");

		//dinamica degli spettri di emissione
		printf("Computing spectrum....\n");
		double***** spectrum = instant_spectrum_NOPT(spin, nFILE);
		printf("Instant spectrum computed!\n");		
		printf("Computing mean spectrum...\n");
		double**** intensity = mean_spectrum_NOPT(spectrum, nFILE, freq, temp); //media sui blocchi temporali
		printf("Computing theoretical IFOs...\n");
		compute_theo_ifo_NOPT(spectrum, intensity, nFILE, temp);
		printf("Computing experimental IFOs...\n");
		compute_exp_ifo_NOPT(intensity, nFILE, temp);
		freeMemory5(spectrum, NREPLICAS, NPT, Size, b_max);
		freeMemory4(intensity, NREPLICAS, NPT, Size);
		printf("Computing Parisi overlaps...\n");
		compute_parisi_overlap_NOPT(spin, nFILE, temp);
		freeMemory5C(spin, NREPLICAS, NPT, Size, b_max);
		
	}

	return 0;
}
