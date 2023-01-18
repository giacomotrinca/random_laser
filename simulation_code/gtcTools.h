//////////////////////FUNCTIONS/////////////////////////////
void welcomeScreen(double v) {
	
	printf("//////////////////////////////////////////////////////\n");
	printf("##SMrandomTetrads analysis: IN SPIN-GLASSES WE TRUST!#\n");
	printf("////////////////// V - %g ///////////////////////////\n", v);
	printf("//////////////////////////////////////////////////////\n");


}
/// SYNTAX CONTROL////////////
void checkSyntax(int a, int n) {
    if(a != n){

		printf("SYNTAX ERROR! Check Readme.txt file for usage.\n");
		exit(SYNTAX_ERROR);
	}

}
//TIME TOOLS
void printTime(clock_t run_time, clock_t start_time) {
  if((double)(run_time - start_time)/CLOCKS_PER_SEC > 3600.) {
    printf("time: %.3lf h---------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n", (double)(run_time - start_time)/(3600*CLOCKS_PER_SEC));
  }else if((double)(run_time - start_time)/CLOCKS_PER_SEC > 60.) {
    printf("time: %.3lf m---------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n", (double)(run_time - start_time)/(60*CLOCKS_PER_SEC));
  }else{
    printf("time: %.3lf s---------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n", (double)(run_time - start_time)/(CLOCKS_PER_SEC));
  }
}


//CHECKING ALLOCATION FUNCTIONS
void checkDouble(double* v) {
	if(v == NULL) {
		printf("Memory allocation problem! Exit.\n");
		exit(MEMORY_ALLOC_ERROR);
	}
}
void checkFile(FILE* f, char* name) {
	if(f == NULL) {
		printf("File %s not found! Exit.\n", name);
		exit(FILE_NOT_FOUND);
	}
}
void checkInt(int* v) {
	if(v == NULL) {
		printf("Memory allocation problem! Exit.\n");
		exit(MEMORY_ALLOC_ERROR);
	}
}
void checkString(char* v) {
	if(v == NULL) {
		printf("Memory allocation problem! Exit.\n");
		exit(MEMORY_ALLOC_ERROR);
	}
}
void checkSentence(char** v) {
	if(v == NULL) {
		printf("Memory allocation problem! Exit.\n");
		exit(MEMORY_ALLOC_ERROR);
	}
}
void checkMat(double** m) {
	if(m == NULL) {
		printf("Memory allocation problem! Exit.\n");
		exit(MEMORY_ALLOC_ERROR);
	}
}
void checkHypMat(double*** m) {
	if(m == NULL) {
		printf("Memory allocation problem! Exit.\n");
		exit(MEMORY_ALLOC_ERROR);
	}
}
void checkHyp2Mat(double**** m) {
	if(m == NULL) {
		printf("Memory allocation problem! Exit.\n");
		exit(MEMORY_ALLOC_ERROR);
	}
}
void checkComplex(complex* s) {
	if(s == NULL) {
		printf("Memory allocation problem! Exit.\n");
		exit(MEMORY_ALLOC_ERROR);
	}
}
void checkComplexMat(complex** m) {
	if(m == NULL) {
		printf("Memory allocation problem! Exit.\n");
		exit(MEMORY_ALLOC_ERROR);
	}
}
void checkHypComplexMat(complex*** m) {
	if(m == NULL) {
		printf("Memory allocation problem! Exit.\n");
		exit(MEMORY_ALLOC_ERROR);
	}
}
void checkHyp2ComplexMat(complex**** m) {
	if(m == NULL) {
		printf("Memory allocation problem! Exit.\n");
		exit(MEMORY_ALLOC_ERROR);
	}
}
////////////////////////////////////////////////////////////
// COMPLEX NUMBERS ALGEBRA /////////////////////////////////
////////////////////////////////////////////////////////////
double Mod2Value(complex z) {  //THIS FUNCTION COMPUTES |Z|^2
	return z.re * z.re + z.im*z.im;
}

double ArgValue(complex z) { //this function computes phi t.c z = |z|exp(i phi)
  return atan2(z.im, z.re);
}
double ModValue(complex z) {  //THIS FUNCTION COMPUTES |Z|
	return sqrt(z.re * z.re + z.im*z.im);
}
complex Conjugate(complex z) {
	complex zc;
	zc.re = z.re;
	zc.im = -z.im;

	return zc;
}

int intPow(int a, int b) {
	
	int p = 1;
	for(int k = 0; k < b; k++) {
		p *= a;
	}
	
	return p;

}
complex ComplexProd(complex a, complex b) {
	complex p;
	p.re = a.re*b.re - a.im*b.im;
	p.im = a.im*b.re + a.re*b.im;

	return p;
}
double* ComplexVec2(complex* vec, int size) { //THIS FUNCTION COMPUTES |z|^2 AS AN Array
	double* v2 = (double*)calloc(size, sizeof(double));

	for(int i = 0; i < size; i++) {
		v2[i] = Mod2Value(vec[i]);
	}

	return v2;
}

double* ComplexVecPhi(complex* vec, int size) {
  double* phi = (double*)calloc(size, sizeof(double));

  for(int i = 0; i < size; i++) {
    phi[i] = ArgValue(vec[i]);
  }

  return phi;
}
/// INPUT/OUTPUT TOOLS ////////////////////////////////////////////////

char getDecision() {
	char dec;
	do{
			dec = getchar();
			if (dec == '\n') dec = getchar();
			while(dec != 'n' && dec != 'N' && dec != 'y' && dec != 'Y')
			{
					//printf("invalid input, enter the choice(y/Y/n/N) again : ");
					dec = getchar();
					if (dec == '\n') dec = getchar();
			}
	}while(dec != 'Y' && dec != 'y' && dec != 'n' && dec != 'N');

	return dec;
}

double* getFrequencies(int size) {

	double* freq = (double*)calloc(size, sizeof(double));
	checkDouble(freq);
	char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(filename);
	sprintf(filename, "frequencies.dat");
	
	FILE* fr = fopen(filename, "r");
	checkFile(fr, filename);
	free(filename);
	
	for(int i = 0; i < size; i++) {
		int tempIndex = 0;
		double tempGain = 0.;
		fscanf(fr, "%d\t%le\t%le\n", &tempIndex, &freq[i], &tempGain);
	}
	fclose(fr);
	
	return freq;

}
//FILE TOOLS /////////////////////////////////////////////////////////
void skipLine(FILE* fr) {

	fscanf(fr, "%*[^\n]\n");

}
void OrganizeFiles(int p, int nRep) {
  //system("mkdir configs");
	//system("mv config_*.dat configs");

	system("mkdir check_tempering");
	system("mv check_tempering_*.dat check_tempering");
if(nRep > 1) {
	system("mkdir PQ_plaqs");
  system("mkdir PQ_plaqs/histograms");

	system("mv EPS_*.dat PQ_plaqs/histograms");
  system("mkdir PQ_plaqs/plaq_overlaps");
  system("mv plaq_overlaps_*.dat PQ_plaqs/plaq_overlaps");
  if(MIXED == 1) {
	   system("mkdir PG");
	   system("mv G_*.dat PG");
   }
	system("mkdir IFOs");
  system("mkdir IFOs/histograms");
	system("mv IFO_temp*.dat IFOs/histograms");
  system("mkdir IFOs/histograms/experimental");
  system("mkdir IFOs/histograms/experimental/equilibrium");
  system("mv IFO_exp_temp*.dat IFOs/histograms/experimental/equilibrium");
  system("mkdir IFOs/histograms/experimental/blocks_indep");
  system("mv IFO_exp_block*.dat IFOs/histograms/experimental/blocks_indep");
  system("mkdir IFOs/histograms/experimental/cumulative");
  system("mv IFO_cumulative_exp_block*.dat IFOs/histograms/experimental/cumulative");
  system("mkdir IFOs/overlaps");
  system("mv ifo_overlaps_temp*.dat IFOs/overlaps");
  system("mkdir IFOs/overlaps/experimental");
  system("mkdir IFOs/overlaps/experimental/equilibrium");
  system("mv experimental_ifo_eq_*.dat IFOs/overlaps/experimental/equilibrium");
  system("mkdir IFOs/overlaps/experimental/blocks_indep");
  system("mv ifo_overlaps_block*.dat IFOs/overlaps/experimental/blocks_indep");
  system("mkdir IFOs/overlaps/experimental/cumulative");
  system("mv ifo_overlaps_cumulative_*.dat IFOs/overlaps/experimental/cumulative");


  system("mkdir PqRE");
  system("mkdir PqRE/histograms");
	system("mv p_q_re_*.dat PqRE/histograms");
  system("mkdir PqRE/overlaps");
  system("mv parisi_overlaps_*dat PqRE/overlaps");
  if(IM_PART_PQ == 1) {
	   system("mkdir PqIM");
     system("mkdir PqIM/histograms");
	   system("mv p_q_im_*.dat PqIM/histograms");
    }
  system("mkdir graphs/PQ_plaqs");
  system("mv graphs/EPS_dist_*.png graphs/PQ_plaqs");
  if(MIXED == 1) {
    system("mkdir graphs/PG");
    system("mv graphs/G_dist_*.png graphs/PG");
  }
  system("mkdir graphs/IFOs");
  system("mv graphs/IFO_dist_temp*.png graphs/IFOs");
  system("mkdir graphs/IFOs/experimental");
  system("mkdir graphs/IFOs/experimental/equilibrium");
  system("mv IFO_dist_exp_temp*.png graphs/IFOs/experimental/equilibrium");
  system("mkdir graphs/IFOs/experimental/blocks_indep");
  system("mv IFO_dist_exp_block*.png graphs/IFOs/experimental/blocks_indep");
  system("mkdir graphs/IFOs/experimental/cumulative");
  system("mv IFO_dist_cumulative_exp_block*.png graphs/IFOs/experimental/cumulative");
  if(IM_PART_PQ == 1){
    system("mkdir graphs/PqIM");
    system("mv graphs/p_q_im_dist_*.png graphs/PqIM");
  }
  system("mkdir graphs/PqRE");
  system("mv graphs/p_q_re_dist_*.png graphs/PqRE");
}
	system("mkdir mean_spectrum");
	system("mv intensity_*.dat mean_spectrum");
	if(p == 1) {
		system("mkdir instant_spectrum");
		system("mv instant_spectra_*.dat instant_spectrum");
    system("mkdir graphs/instant_spectra");
		system("mv instant_spectra_*.png graphs/instant_spectra");
		system("mkdir graphs/instant_phases");
		system("mv instant_phase_*.png graphs/instant_phases");

	}
	system("rm *.p");

	system("mkdir specific_heat");
	system("mv specific_heat_*.dat specific_heat");

	system("mkdir sources");
	system("mv *.cpp sources");
	system("mv *.cu sources");
	system("mv *.h sources");

	//system("mkdir stuffs");
	//system("mv *.dat stuffs");
	//system("mv *.txt stuffs");

	system("mkdir graphs/mean_spectrum");
	system("mv graphs/spectrum_*.png graphs/mean_spectrum");
	system("mkdir graphs/specific_heat");
	system("mv graphs/specific_heat_*.png graphs/specific_heat");

	system("mkdir graphs/check_tempering");
	system("mv check_tempering_*.png graphs/check_tempering");

}

complex* readConfigFILE(int r, int iter, int nRows, int choose_mod) { //THIS FUNCTION READ A SINGLE CONFIGURATION FILE

	complex* conf = (complex*)calloc(nRows, sizeof(complex));
	checkComplex(conf);

	double T;
	int I;

	char* fileName = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(fileName);
	if(choose_mod == 0) {
    sprintf(fileName, "config_nrep%d_iter_%d.dat", r, iter);
  }else if(choose_mod == 1) {
    sprintf(fileName, "configs/config_nrep%d_iter_%d.dat", r, iter);
  }

	FILE* fp = fopen(fileName, "r");
	checkFile(fp, fileName);


	for(int i = 0; i < nRows; i++) {
			fscanf(fp, "%lf %d %le %le", &T, &I, &conf[i].re, &conf[i].im);
	}

	fclose(fp);
	free(fileName);

	return conf;

}



complex**** allow_conf(complex**** old, int i_max, int j_max) {

	complex**** new_conf = (complex****)calloc(NREPLICAS, sizeof(complex***));
	checkHyp2ComplexMat(new_conf);
	//conf[r][i][k][j]
	for(int r = 0; r < NREPLICAS; r++) {
		new_conf[r] = (complex***)calloc(NPT, sizeof(complex**));
		checkHypComplexMat(new_conf[r]);
		
		for(int k = 0; k < NPT; k++) {
			new_conf[r][k] = (complex**)calloc(j_max, sizeof(complex*));
			checkComplexMat(new_conf[r][k]);
			
			for(int j = 0; j < j_max; j++) {
				new_conf[r][k][j] = (complex*)calloc(i_max, sizeof(complex));
				checkComplex(new_conf[r][k][j]);
				
				for(int i = 0; i < i_max; i++) {
				
					new_conf[r][k][j][i] = old[r][i][k][j];
				}
			}
		}
	}
	
	return(new_conf);

}

double*** readDynamics(int r, int nPT, int IT_MAX, int dITER, int size, int choose_mod) {//THIS FUNCTION READ PARALLEL_TEMPERING%d.dat

	double*** data = (double***)calloc(2, sizeof(double**));
	checkHypMat(data);
	data[0] = (double**)calloc(IT_MAX+dITER, sizeof(double*));  //Acceptation rates
	checkMat(data[0]);
	data[1] = (double**)calloc(IT_MAX+dITER, sizeof(double*)); //Energies
	checkMat(data[1]);


	char* fileName = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(fileName);
	if(choose_mod == 0) {
    sprintf(fileName, "parallel_tempering%d.dat", r);
  }else if(choose_mod == 1) {
    sprintf(fileName, "stuffs/parallel_tempering%d.dat", r);
  }

	FILE* fr = fopen(fileName, "r");
	checkFile(fr, fileName);
	free(fileName);
	//TEMP VARIABLES
	int mcs;
	double T, a, e;
	int nCol = nPT - 1;


	for(int j = 0; j < IT_MAX + dITER; j++) { // CICLO SULLE RIGHE
		fscanf(fr, " %d \t", &mcs);
		//printf("\n%d\n\n", mcs);

		data[0][j] = (double*)calloc(nCol, sizeof(double));
		checkDouble(data[0][j]);
		data[1][j] = (double*)calloc(nCol, sizeof(double));
		checkDouble(data[1][j]);

		for(int i = 0; i < nCol; i++) { //CICLO SULLE COLONNE
			fscanf(fr, " %lg %le %le \t", &T, &a, &e);
			//printf("%lf %le %le \n\n", T,a,e);
			//LE TEMPERATURE SONO INDICIZZATE AL CONTRARIO: LE RIMETTO A POSTO....
			data[0][j][nCol-1-i] = a;
			data[1][j][nCol-1-i] = size * e;

			//printf("%le %le\t", data[0][j][i], data[1][j][i]);

		}//FINE CICLO SULLE COLONNE
    fscanf(fr,"\n");
		//printf("\n");
	} //FINE CICLO SULLE RIGHE

	fclose(fr);


	return data;


}



void PlotSpectrumScript(int r, int nPT, int N){
	char* fileName = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(fileName);
	char* command = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(command);
	sprintf(fileName, "plot_spectrum_script_N%d_nrep%d.p", N, r);
	sprintf(command, "gnuplot plot_spectrum_script_N%d_nrep%d.p", N, r);
	FILE* script = fopen(fileName, "w");
	checkFile(script, fileName);
	fprintf(script, "set terminal png\n");
	fprintf(script, "set xlabel 'w'\n");
	fprintf(script, "set ylabel 'I(w)'\n");
	fprintf(script, "unset key\n");
	fprintf(script, "set output 'spectrum_N%d_nrep%d.png' \n", N, r);
	fprintf(script, "plot 'intensity_nrep%d_temp%d.dat' w lp pt 7, ", r, 0);

	for(int i = 1; i < nPT-1; i++) {
		fprintf(script, "'intensity_nrep%d_temp%d.dat' w lp pt 7, ", r, i);
	}

	fprintf(script, "'intensity_nrep%d_temp%d.dat' w lp pt 7", r, nPT-1);
	free(fileName);
	fclose(script);
	system(command);
	char* move = (char*)calloc(STRING_SIZE,sizeof(char));
	checkString(move);
	sprintf(move, "mv spectrum_N%d_nrep%d.png graphs", N, r);
	system(move);
	free(command);
	free(move);

}
void PlotSpecificHeatScript(int r, int N) {
	char* fileName = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(fileName);
	char* command = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(command);
	sprintf(fileName, "plot_specific_heat_script_N%d_nrep%d.p", N, r);
	sprintf(command, "gnuplot plot_specific_heat_script_N%d_nrep%d.p", N, r);
	FILE* script = fopen(fileName, "w");
	checkFile(script, fileName);
	fprintf(script, "set terminal png\n");
	fprintf(script, "set xlabel 'T'\n");
	fprintf(script, "set ylabel 'C_v(T)'\n");
	fprintf(script, "unset key\n");
	fprintf(script, "set output 'specific_heat_N%d_nrep%d.png' \n", N, r);
	fprintf(script, "plot 'specific_heat_nrep%d_size%d.dat' w lp pt 7,'specific_heat_nrep%d_size%d.dat' w errorbars\n", r, N, r, N);

	free(fileName);
	fclose(script);
	system(command);
	char* move = (char*)calloc(STRING_SIZE,sizeof(char));
	checkString(move);
	sprintf(move, "mv specific_heat_N%d_nrep%d.png graphs", N, r);
	system(move);
	free(command);
	free(move);

}
void PrintSpecificHeat(double* sh, double* errSh, int size,int r,double* temp, int N) {
	FILE* fw;
	char* fileName = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(fileName);


	sprintf(fileName, "specific_heat_nrep%d_size%d.dat", r, N);

	fw = fopen(fileName, "w");
	checkFile(fw, fileName);

	for(int i = 0; i < size; i++) {
		fprintf(fw, "%lf\t%le\t%le\n", temp[i], sh[i], errSh[i]);
	}

	free(fileName);
	fclose(fw);


}


void PrintPQ(double* q, double* p_q, int size, int k, int p, double mean,double devStd, double skew,double errS, double kurt, double errK, int block) {
	char* fileName = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(fileName);

	if(p == 0) {

    if(block == -1) {
      sprintf(fileName, "p_q_re_temp%d.dat", k);
    }else{
      sprintf(fileName, "p_q_re_temp%d_block%d.dat", k, block);
    }
  }else if(p == 1) {
    if(block == -1) {
      sprintf(fileName, "p_q_im_temp%d.dat", k);
    }else{
      sprintf(fileName, "p_q_im_temp%d_block%d.dat", k, block);
    }
  }

	FILE* fw = fopen(fileName, "w");

  fprintf(fw, "#mean: %g +/- %g\tskewness: %g +/- %g kurtosis: %g +/- %g\n", mean,devStd, skew, errS, kurt, errK);
	for(int i = 0; i < size; i++) {
		fprintf(fw, "%le\t%le\n", q[i], p_q[i]);
	}

	fclose(fw);
	free(fileName);


}
void PrintIFO(double* q, double* p_q, int size, int k) {
	char* fileName = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(fileName);

	sprintf(fileName, "IFO_temp%d.dat", k);
	FILE* fw = fopen(fileName, "w");

	for(int i = 0; i < size; i++) {
		fprintf(fw, "%le\t%le\n", q[i], p_q[i]);
	}

	fclose(fw);
	free(fileName);


}

void PrintExpIFO(double* q, double* p_q, int size, int k, int block, int c) {
	char* fileName = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(fileName);

	if(block == -1) {
    sprintf(fileName, "IFO_exp_temp%d.dat", k);
  }else{
    if(c == 0) {
      sprintf(fileName, "IFO_exp_block%d_temp%d.dat", block, k);
    }else if(c == 1) {
      sprintf(fileName, "IFO_cumulative_exp_block%d_temp%d.dat", block, k);
    }
  }
	FILE* fw = fopen(fileName, "w");

	for(int i = 0; i < size; i++) {
		fprintf(fw, "%le\t%le\n", q[i], p_q[i]);
	}

	fclose(fw);
	free(fileName);


}



void PlotPQScript(int k, int p, double* temp, double mean, double devStd, double skew, double errS, double kurt, double errK, int block) {
	char* fileName = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(fileName);
	char* command = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(command);

	if(p == 0) {
		if(block == -1) {
      sprintf(fileName, "plot_p_q_re_script_temp%d.p", k);
		  sprintf(command, "gnuplot plot_p_q_re_script_temp%d.p", k);
    }else{
      sprintf(fileName, "plot_p_q_re_script_temp%d_block%d.p", k, block);
		  sprintf(command, "gnuplot plot_p_q_re_script_temp%d_block%d.p", k, block);
    }
	}
	if(p == 1) {
		if(block == -1) {
      sprintf(fileName, "plot_p_q_im_script_temp%d.p", k);
		  sprintf(command, "gnuplot plot_p_q_im_script_temp%d.p", k);
    }else{
      sprintf(fileName, "plot_p_q_im_script_temp%d_block%d.p", k, block);
		  sprintf(command, "gnuplot plot_p_q_im_script_temp%d_block%d.p", k, block);
    }
	}
	FILE* script = fopen(fileName, "w");
	checkFile(script, fileName);

	fprintf(script, "set terminal png\n");
	fprintf(script, "set xlabel 'q'\n");
	fprintf(script, "set ylabel 'P(q)'\n");
  char* title = (char*)calloc(10*STRING_SIZE, sizeof(char));
  checkString(title);
  sprintf(title,"mean: %g+/-%g, skewness: %g+/-%g kurtosis %g+/-%g", mean,devStd, skew, errS, kurt, errK );
  fprintf(script, "set title '%s' font ',8'\n", title);
  free(title);
	//fprintf(script, "unset key\n");
	if(p == 0) {
    if(block == -1) {
      fprintf(script, "set output 'p_q_re_dist_temp%d.png'\n", k);
    }else{
      fprintf(script, "set output 'p_q_re_dist_temp%d_block%d.png'\n", k, block);
    }
  }
	if(p == 1) {
    if(block == -1) {
      fprintf(script, "set output 'p_q_im_dist_temp%d.png'\n", k);
    }else{
      fprintf(script, "set output 'p_q_im_dist_temp%d_block%d.png'\n", k, block);
    }
  }

	//fprintf(script, "set xrange [-1:1]\n");
	//fprintf(script, "set nonlinear y via log10(y) inverse 10**y\n");
	fprintf(script, "set logscale y 10\n");
	fprintf(script, "set yrange [0.0001:]\n");
	//fprintf(script, "set offset graph 0.05,0.05,0.05,0.0\n");

	//fprintf(script, "set boxwidth width*0.9\n");
	fprintf(script, "set style fill solid 0.5\n");
	//fprintf(script, "set tics out nomirror\n");
	//fprintf(script, "stats 'p_q_temp%d.dat'\n", k);
	//fprintf(script, "sum=STATS_records\n");
	if(p == 0) {
    if(block == -1) {
      fprintf(script, "plot 'p_q_re_temp%d.dat' w boxes lc rgb'red' title 'T=%.2lf'",k, temp[k]);
    }else{
      fprintf(script, "plot 'p_q_re_temp%d_block%d.dat' w boxes lc rgb'red' title 'T=%.2lf'",k,block, temp[k]);
    }
  }
	if(p == 1) {
    if(block == -1) {
      fprintf(script, "plot 'p_q_im_temp%d.dat' w boxes lc rgb'red' title 'T=%.2lf'",k, temp[k]);
    }else{
      fprintf(script, "plot 'p_q_im_temp%d_block%d.dat' w boxes lc rgb'red' title 'T=%.2lf'",k,block, temp[k]);
    }
  }

	free(fileName);
	fclose(script);
	system(command);
	char* move = (char*)calloc(STRING_SIZE,sizeof(char));
	checkString(move);
	if(p == 0) {
    if(block == -1) {
      sprintf(move, "mv p_q_re_dist_temp%d.png graphs", k);
    }else{
      sprintf(move, "mv p_q_re_dist_temp%d_block%d.png graphs", k, block);
    }
  }
	if(p == 1) {
    if(block == -1) {
      sprintf(move, "mv p_q_im_dist_temp%d.png graphs", k);
    }else{
      sprintf(move, "mv p_q_im_dist_temp%d_block%d.png graphs", k, block);
    }
  }
	system(move);
	free(command);
	free(move);

}

void PlotIFOScript(int k, double* temp) {
	char* fileName = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(fileName);
	char* command = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(command);
	sprintf(fileName, "plot_IFO_script_temp%d.p", k);
	sprintf(command, "gnuplot plot_IFO_script_temp%d.p", k);
	FILE* script = fopen(fileName, "w");
	checkFile(script, fileName);

	fprintf(script, "set terminal png\n");
	fprintf(script, "set xlabel 'C'\n");
	fprintf(script, "set ylabel 'P(C)'\n");
	//fprintf(script, "unset key\n");
	fprintf(script, "set output 'IFO_dist_temp%d.png'\n", k);
	//fprintf(script, "set xrange [-1:1]\n");
	//fprintf(script, "set nonlinear y via log10(y) inverse 10**y\n");
	fprintf(script, "set logscale y 10\n");
	fprintf(script, "set yrange [0.0001:]\n");
	//fprintf(script, "set offset graph 0.05,0.05,0.05,0.0\n");

	//fprintf(script, "set boxwidth width*0.9\n");
	fprintf(script, "set style fill solid 0.5\n");
	//fprintf(script, "set tics out nomirror\n");
	//fprintf(script, "stats 'p_q_temp%d.dat'\n", k);
	//fprintf(script, "sum=STATS_records\n");
	fprintf(script, "plot 'IFO_temp%d.dat' w boxes lc rgb'red' title 'T = %.2lf'",k, temp[k]);

	free(fileName);
	fclose(script);
	system(command);
	char* move = (char*)calloc(STRING_SIZE,sizeof(char));

	checkString(move);
	sprintf(move, "mv IFO_dist_temp%d.png graphs", k);
	system(move);
	free(command);
	free(move);

}

void PlotExpIFOScript(int k, double* temp, int block, int c) {
	char* fileName = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(fileName);
	char* command = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(command);
	sprintf(fileName, "plot_IFO_script_temp%d.p", k);
	sprintf(command, "gnuplot plot_IFO_script_temp%d.p", k);
	FILE* script = fopen(fileName, "w");
	checkFile(script, fileName);

	fprintf(script, "set terminal png\n");
	fprintf(script, "set xlabel 'C'\n");
	fprintf(script, "set ylabel 'P(C)'\n");
	//fprintf(script, "unset key\n");
	if(block == -1) {
    fprintf(script, "set output 'IFO_dist_exp_temp%d.png'\n", k);
  }else{
    if(c == 0) {
      fprintf(script, "set output 'IFO_dist_exp_block%d_temp%d.png'\n", block, k);
    }else if(c == 1) {
      fprintf(script, "set output 'IFO_dist_cumulative_exp_block%d_temp%d.png'\n", block, k);

    }
  }
	//fprintf(script, "set xrange [-1:1]\n");
	//fprintf(script, "set nonlinear y via log10(y) inverse 10**y\n");
	fprintf(script, "set logscale y 10\n");
	fprintf(script, "set yrange [0.0001:]\n");
	//fprintf(script, "set offset graph 0.05,0.05,0.05,0.0\n");

	//fprintf(script, "set boxwidth width*0.9\n");
	fprintf(script, "set style fill solid 0.5\n");
	//fprintf(script, "set tics out nomirror\n");
	//fprintf(script, "stats 'p_q_temp%d.dat'\n", k);
	//fprintf(script, "sum=STATS_records\n");
	if(block == -1) {
    fprintf(script, "plot 'IFO_exp_temp%d.dat' w boxes lc rgb'red' title 'T = %.2lf'",k, temp[k]);
  }else{
    if(c == 0) {
      fprintf(script, "plot 'IFO_exp_block%d_temp%d.dat' w boxes lc rgb'red' title 'T = %.2lf'",block, k, temp[k]);
    }else if(c == 1) {
      fprintf(script, "plot 'IFO_cumulative_exp_block%d_temp%d.dat' w boxes lc rgb'red' title 'T = %.2lf'",block, k, temp[k]);
    }
  }
	free(fileName);
	fclose(script);
	system(command);

	free(command);
	//free(move);

}

void PlotGScript(int k, double* temp) {
	char* fileName = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(fileName);
	char* command = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(command);
	sprintf(fileName, "plot_G_script_temp%d.p", k);
	sprintf(command, "gnuplot plot_G_script_temp%d.p", k);
	FILE* script = fopen(fileName, "w");
	checkFile(script, fileName);

	fprintf(script, "set terminal png size 1024,1024\n");
	fprintf(script, "set xlabel 'G'\n");
	fprintf(script, "set ylabel 'P(G)'\n");
	//fprintf(script, "unset key\n");
	fprintf(script, "set output 'G_dist_temp%d.png'\n", k);
	//fprintf(script, "set xrange [-1:1]\n");
	//fprintf(script, "set nonlinear y via log10(y) inverse 10**y\n");
	fprintf(script, "set logscale y 10\n");
	fprintf(script, "set yrange [0.0001:]\n");
	//fprintf(script, "set offset graph 0.05,0.05,0.05,0.0\n");

	//fprintf(script, "set boxwidth width*0.9\n");
	fprintf(script, "set style fill solid 0.5\n");
	//fprintf(script, "set tics out nomirror\n");
	//fprintf(script, "stats 'p_q_temp%d.dat'\n", k);
	//fprintf(script, "sum=STATS_records\n");
	fprintf(script, "plot 'G_temp%d.dat' w boxes lc rgb'red' title 'T = %.2lf'",k, temp[k]);

	free(fileName);
	fclose(script);
	system(command);
	char* move = (char*)calloc(STRING_SIZE,sizeof(char));

	checkString(move);
	sprintf(move, "mv G_dist_temp%d.png graphs", k);
	system(move);
	free(command);
	free(move);

}
void plotInstantPhase(int r, int k, int iter) {
  char* fileName = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(fileName);

  sprintf(fileName, "instant_spectra_nrep%d_temp%d_iter%d.dat", r, k, iter);
  char* plotFilename = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(plotFilename);

  char* imgFilename = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(imgFilename);
  sprintf(imgFilename, "instant_phase_nrep%d_temp%d_iter%d.png", r, k, iter);

  char* tempfile_name = (char*)calloc(STRING_SIZE, sizeof(char*));
  checkString(tempfile_name);
  sprintf(tempfile_name, "instant_phase_plot.p");

  FILE* plot = fopen(tempfile_name, "w");
  checkFile(plot, tempfile_name);
  free(tempfile_name);

  fprintf(plot, "set term png\n");
  fprintf(plot, "set xlab 'w'\n");
  fprintf(plot, "unset key\n");
  fprintf(plot, "set yrange[-pi:pi]\n");
  fprintf(plot, "set ytics pi\n");
  fprintf(plot, "set format y '%%.0P{/Symbol p}'\n");
  fprintf(plot, "set ylab 'phi(w)'\n");
  fprintf(plot, "set xrange[0:1]\n");
  //fprintf(plot, "set yrange[0:0.4]\n");
  fprintf(plot, "set output '%s'\n", imgFilename);
  fprintf(plot, "set title 't = %d'\n", iter);
  fprintf(plot, "plot '%s' u 1:3 w p pt 7 lt rgb 'black'\n", fileName);

  fclose(plot);

  system("gnuplot instant_phase_plot.p");

  system("rm instant_phase_plot.p");

  free(fileName);
  free(imgFilename);

}

void plotInstantSpectra(int r, int k, int iter) {
  char* fileName = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(fileName);

  sprintf(fileName, "instant_spectra_nrep%d_temp%d_iter%d.dat", r, k, iter);
  char* plotFilename = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(plotFilename);

  char* imgFilename = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(imgFilename);
  sprintf(imgFilename, "instant_spectra_nrep%d_temp%d_iter%d.png", r, k, iter);
  char* tempfile_name = (char*)calloc(STRING_SIZE, sizeof(char*));
  checkString(tempfile_name);
  sprintf(tempfile_name, "instant_spectra_plot.p");
  FILE* plot = fopen(tempfile_name, "w");
  checkFile(plot, tempfile_name);
  free(tempfile_name);
  fprintf(plot, "set term png\n");
  fprintf(plot, "set xlab 'w'\n");
  fprintf(plot, "unset key\n");
  fprintf(plot, "set ylab 'A^2(w)'\n");
  fprintf(plot, "set xrange[0:1]\n");
  fprintf(plot, "set yrange[0:0.4]\n");
  fprintf(plot, "set output '%s'\n", imgFilename);
  fprintf(plot, "set title 't = %d'\n", iter);
  fprintf(plot, "plot '%s' w lp pt 12 lt rgb 'blue'\n", fileName);

  fclose(plot);

  system("gnuplot instant_spectra_plot.p");

  system("rm instant_spectra_plot.p");

  free(fileName);
  free(imgFilename);

}

void PrintG(double* q, double* p_q, int size, int k) {
	char* fileName = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(fileName);

	sprintf(fileName, "G_temp%d.dat", k);
	FILE* fw = fopen(fileName, "w");

	for(int i = 0; i < size; i++) {
		fprintf(fw, "%le\t%le\n", q[i], p_q[i]);
	}

	fclose(fw);
	free(fileName);


}

void PrintEPS(double* q, double* p_q, int size, int k) {
	char* fileName = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(fileName);

	sprintf(fileName, "EPS_temp%d.dat", k);
	FILE* fw = fopen(fileName, "w");

	for(int i = 0; i < size; i++) {
		fprintf(fw, "%le\t%le\n", q[i], p_q[i]);
	}

	fclose(fw);
	free(fileName);


}

void PlotEPSScript(int k, double* temp) {
	char* fileName = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(fileName);
	char* command = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(command);
	sprintf(fileName, "plot_EPS_script_temp%d.p", k);
	sprintf(command, "gnuplot plot_EPS_script_temp%d.p", k);
	FILE* script = fopen(fileName, "w");
	checkFile(script, fileName);

	fprintf(script, "set terminal png\n");
	fprintf(script, "set xlabel 'Q'\n");
	fprintf(script, "set ylabel 'P(Q)'\n");
	//fprintf(script, "unset key\n");
	fprintf(script, "set output 'EPS_dist_temp%d.png'\n", k);
	fprintf(script, "set xrange [-0.2:]\n");
	//fprintf(script, "set nonlinear y via log10(y) inverse 10**y\n");
	fprintf(script, "set logscale y 10\n");
	fprintf(script, "set yrange [0.0001:]\n");
	//fprintf(script, "set offset graph 0.05,0.05,0.05,0.0\n");

	//fprintf(script, "set boxwidth width*0.9\n");
	fprintf(script, "set style fill solid 0.5\n");
	//fprintf(script, "set tics out nomirror\n");
	//fprintf(script, "stats 'p_q_temp%d.dat'\n", k);
	//fprintf(script, "sum=STATS_records\n");
	fprintf(script, "plot 'EPS_temp%d.dat' w boxes lc rgb'red' title 'T = %.2lf'",k, temp[k]);

	free(fileName);
	fclose(script);
	system(command);
	char* move = (char*)calloc(STRING_SIZE,sizeof(char));

	checkString(move);
	sprintf(move, "mv EPS_dist_temp%d.png graphs", k);
	system(move);
	free(command);
	free(move);

}

//PIPES FOR DATA PLOTTING RUNTIME/////////////////////////////////////

void TemperingCheckGnuPipe(FILE* p, int r, int nPT) {
	char** fileName = (char**)calloc(nPT, sizeof(char*));
	checkSentence(fileName);



	//fprintf(p, "set terminal x11\n");
	fprintf(p, "unset key\n");
  fprintf(p, "set terminal png\n");
  char* img_name = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(img_name);
  sprintf(img_name, "check_tempering_nrep%d.png", r);
  fprintf(p, "set output '%s'\n", img_name);
	fprintf(p, "set logscale x 2\n");
  //fprintf(p, "unset ytics\n");
  fprintf(p, "set y2tics 0.1\n");
	fprintf(p, "set xlabel'Monte Carlo Steps'\n");
	fprintf(p, "set ylabel 'E(t)'\n");
	fileName[0] = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(fileName[0]);
	sprintf(fileName[0], "check_tempering_E_nrep_%d_temp%d.dat", r, 0);

	fprintf(p, "plot '%s' w lp pt 7, '%s' w errorbars, ", fileName[0], fileName[0]);

	for(int i = 1; i < nPT-1; i++) {
		fileName[i] = (char*)calloc(STRING_SIZE, sizeof(char));
		checkString(fileName[i]);
		sprintf(fileName[i], "check_tempering_E_nrep_%d_temp%d.dat", r, i);
		fprintf(p, "'%s' w lp pt 7, '%s' w errorbars, ", fileName[i], fileName[i]);
	}

	fileName[nPT-1] = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(fileName[nPT-1]);
	sprintf(fileName[nPT-1], "check_tempering_E_nrep_%d_temp%d.dat", r, nPT-1);
	fprintf(p, "'%s' w lp pt 7, '%s' w errorbars\n", fileName[nPT-1], fileName[nPT-1]);


	fflush(p);

	free(fileName);

}

//ARRAY TOOLS/////////////////////////////////////////////////
static int compare (const void * a, const void * b){
  if (*(double*)a > *(double*)b) return 1;
  else if (*(double*)a < *(double*)b) return -1;
  else return 0;
}
double ArraySum(double* vec, int size) { //THIS FUNCTION SUMS THE ELEMENTS OF AN ARRAY
	double m = .0;

	for(int i = 0; i < size; i++) {
		m += vec[i];
	}

	return m;
}

double MaxVec(double* vec, int size) {
	
	double max = vec[0];
	
	for(int i = 1; i < size;i++) {
		
		if(max < vec[i]) {
			max = vec[i];
		}
	
	}
	
	return max;


}

double* DoubleAddElement(double* vec, double e, int size) {
  double* r = (double*)calloc(size+1, sizeof(double));
  checkDouble(r);
  for(int i = 0; i < size; i++) {
    r[i] = vec[i];
  }
  r[size] = e;
  return r;
}

int* IntAddElement(int* vec, int e, int size) {
  int* r = (int*)calloc(size+1, sizeof(int));
  checkInt(r);
  for(int i = 0; i < size; i++) {
    r[i] = vec[i];
  }
  r[size] = e;
  return r;
}

double* sqVec(double* vec, int size) {
	double* v = (double*)calloc(size, sizeof(double));
	checkDouble(v);
	for(int i = 0; i < size; i++) {
		v[i] = vec[i]*vec[i];
	}
	return v;
}
double* TranslateVec(double* vec, double value, int size){
  double* tVec = (double*)calloc(size, sizeof(double));
  checkDouble(tVec);

  for(int i = 0; i < size; i++) {
        tVec[i] = vec[i] - value;
  }

  return tVec;
}
double* sumTwoArrays(double* v, double* w, int size) { //THIS FUNCTION SUM TWO ARRAYS

	double* s = (double*)calloc(size, sizeof(double));
	checkDouble(s);

	for(int i = 0; i < size; i++) {
		s[i] = v[i] + w[i];
	}

	return s;
}
double MeanVec(double* vec, int size) { //THIS FUNCTION COMPUTER THE MEAN OF AN ARRAY

	double m = .0;

	for(int i = 0; i < size; i++) {
		m += vec[i];
	}

	return m/size;
}
double StandardDeviation(double* vec, double mean, int size) { //THIS FUNCTION COMPUTER THE STANDARD DEVIATION OF AN ARRAY

  /*
  double sigma = .0;

	for(int i = 0; i < size; i++) {
		sigma += (vec[i]-mean)*(vec[i]-mean);
	}
	sigma = sigma/(size-1);

	sigma = sqrt(sigma);
  */

  double* prod = (double*)calloc(size, sizeof(double));
  checkDouble(prod);

  for(int i = 0; i < size; i++) {
    prod[i] = vec[i]*vec[i];
  }

  double m2 = ArraySum(prod, size) / size;
  free(prod);

  double s = m2 - mean*mean;

	return sqrt(s/size);
}

double mean_of_distribution (double* x, double* p, double bin_size, int size) {

  double m = 0.;
  for(int i = 0; i < size; i++) {
    m += x[i]*p[i]*bin_size;
  }
  return m;
}
double variance (double* x, double* p, double bin_size, int size) {

  double m = 0.;
  double m2 = 0.;
  for(int i = 0; i < size; i++) {
    m += x[i]*p[i]*bin_size;
  }
  for(int i = 0; i < size; i++) {

   m2 += (x[i]-m)*(x[i]-m)*p[i]*bin_size;
  }
  return m2;
}
double Skewness(double* x, double* p, double bin_size, int size) {


  double m3 = 0.;
  double m = 0.;
  double m2 = 0.;
  for(int i = 0; i < size; i++) {
    m += x[i]*p[i]*bin_size;
  }
  for(int i = 0; i < size; i++) {

   m2 += (x[i]-m)*(x[i]-m)*p[i]*bin_size;
  }
  for(int i = 0; i < size; i++) {
    m3 += (x[i]-m)*(x[i]-m)*(x[i]-m)*p[i]*bin_size;
  }

  double s = m3/(sqrt(m2)*sqrt(m2)*sqrt(m2));

  return s;
}

double errSkewness(double* x, double* p, double bin_size, int size, double meanSkew) {

  double err = 0.;
  double m = 0.;


  for(int i = 0; i < size; i++) {
    m += x[i]*p[i]*bin_size;
  }

  for(int i = 0; i < size; i++) {
    err += ((x[i] - m) - meanSkew)*((x[i] - m) - meanSkew)*((x[i] - m) - meanSkew)*((x[i] - m) - meanSkew)*((x[i] - m) - meanSkew)*((x[i] - m) - meanSkew)*p[i]*bin_size;
  }


  return sqrt(err/size);
}

double errKurtosis(double* x, double* p, double bin_size, int size, double meanKurt) {
  double err = 0.;
  double m = 0.;

  for(int i = 0; i < size; i++) {
    m += x[i]*p[i]*bin_size;
  }

  for(int i = 0; i < size; i++) {
    err += ((x[i] - m) - meanKurt)*((x[i] - m) - meanKurt)*((x[i] - m) - meanKurt)*((x[i] - m) - meanKurt)*((x[i] - m) - meanKurt)*((x[i] - m) - meanKurt)*((x[i] - m) - meanKurt)*((x[i] - m) - meanKurt)*p[i]*bin_size;
  }



  return sqrt(err/size);
}


double Kurtosis(double*x, double* p, double bin_size, int size) {


  double m4 = 0.;
  double m2 = 0.;
  double m = 0.;
  for(int i = 0; i < size; i++) {
    m += x[i]*p[i]*bin_size;
  }
  for(int i = 0; i < size; i++) {
   m2 += (x[i]-m)*(x[i]-m)*p[i]*bin_size;

  }
  for(int i = 0; i < size; i++) {
    m4 += (x[i]-m)*(x[i]-m)*(x[i]-m)*(x[i]-m)*p[i]*bin_size;
  }

  double k = m4/(m2*m2);

  return k;

}

double est_variance(double* vec, int size) {



	double m2 = 0.;
	
	double m = MeanVec(vec, size);
	
	for(int i = 0; i < size; i++) {
		m2 += (vec[i] - m)*(vec[i] - m) / size;
	
	}
	
	return m2;

}

double est_skewness(double* vec, int size) {

	double m3 = 0.;
	
	double m = MeanVec(vec, size);
	double m2 =est_variance(vec, size);
	
	for(int i = 0; i < size; i++) {
	
		m3 += (vec[i] - m)*(vec[i] - m)*(vec[i] - m) / size;
	}
	
	return m3/(sqrt(m2)*sqrt(m2)*sqrt(m2));


}

double est_kurtosis(double* vec, int size) {

	double m4 = 0.;
	
	double m = MeanVec(vec, size);
	double m2 = est_variance(vec, size);
	
	
	for(int i = 0; i < size; i++) {
	
		m4 += (vec[i] - m)*(vec[i] - m)*(vec[i] - m)*(vec[i] - m) / size;
	
	}
	
	return m4 / (m2 * m2);

}






double* ScalarMoltVec(double* vec, int size, double b) {
  double* r = (double*)calloc(size, sizeof(double));
  checkDouble(r);
  for(int i = 0; i < size; i++) {
    r[i] = b * vec[i];
  }
  return r;
}

double* ScaleVector(double* vec, int size, double scale) { //THIS FUNCTION RESCALE AN ARRAY OF A FACTOR
	double* Resized = (double*)calloc(size, sizeof(double));
	checkDouble(Resized);

	for(int i = 0; i < size; i++) {
		Resized[i] = vec[i]/scale;
	}

	return Resized;
}
double** TransposeMat(double**mat, int R, int C) { //THIS FUNCTION TRANSPOSE A MATRIX

	double** TMat = (double**)calloc(C, sizeof(double*));
	checkMat(TMat);

	for(int i = 0; i < C; i++) {
		TMat[i] = (double*)calloc(R, sizeof(double));
	}

	for(int j = 0; j < C; j++) {
		for(int i = 0; i < R; i++) {
			TMat[j][i] = mat[i][j];
		}
	}

	return TMat;
}
double* DeBlock(double* vec, int BlockSize, int nBlock, int k) { //THIS FUNCTION CHOOSE A BLOCK FROM AN ARRAY
	double* v = (double*)calloc(BlockSize, sizeof(double));
	checkDouble(v);

	int l = 0;

	for(int i = 0; i < nBlock*BlockSize; i++) {
		if(i >= k*BlockSize && l< BlockSize) {
			v[l] = vec[i];
			l++;
		}
	}
	return v;
}
complex* ComplexDeBlock(complex* vec, int BlockSize, int nBlock, int k) { //THIS FUNCTION CHOOSE A BLOCK FROM AN ARRAY
	complex* v = (complex*)calloc(BlockSize, sizeof(complex));
	checkComplex(v);

	int l = 0;

	for(int i = 0; i < nBlock*BlockSize; i++) {
		if(i >= k*BlockSize && l< BlockSize) {
			v[l].re = vec[i].re;
			v[l].im = vec[i].im;
			l++;
		}
	}
	//printf("OK COMPLEXDEBLOCK\n");
	return v;
}
void ReBlock(double* vec, double* result, int size, int nBlock){
	for(int j = 0; j < size; j++) {
			result[nBlock*size + j] = vec[j];
	}
}

void PrintMat(double** mat, int R, int C) {
	for(int i = 0; i < R; i ++) {
		for(int j = 0; j < C; j ++) {
			printf("%lf\t", mat[i][j]);
		}
		printf("\n");
	}
}
void PrintComplexMat(complex** mat, int R, int C) {
	printf("Real Part\n");
	for(int i = 0; i < R; i ++) {
		for(int j = 0; j < C; j ++) {
			printf("%lf\t", mat[i][j].re);
		}
		printf("\n");
	}
	printf("Imaginary part\n");
	for(int i = 0; i < R; i ++) {
		for(int j = 0; j < C; j ++) {
			printf("%lf\t", mat[i][j].im);
		}
		printf("\n");
	}
}

//////////////////// HISTOGRAMS //////////////////////////////////////
void HistogramQ(double* vec, int size, int temp, int p, double* T, int block) {

    double* counter = (double*)calloc(HISTBARS, sizeof(double));
    checkDouble(counter);

    double* qMean = (double*)calloc(HISTBARS, sizeof(double));
    checkDouble(qMean);

    double min = -1.01;
    double max = 1.01;

    for(int i = 0; i < size; i++) {
      if(vec[i] > 1 || vec[i] < -1) {
        printf("OUT OF RANGE!");
      }
    }

    /*
    for(int i = 0; i < size; i++) {

        if(min>vec[i]){
            min = vec[i];
        }
        if(max<vec[i]) {
            max = vec[i];
        }

    }
*/
    double bin_size = (max - min)/HISTBARS;

    for(int i = 0; i < HISTBARS; i++) {

        for(int j = 0; j < size; j++) {

            if(vec[j] > min + i*bin_size && vec[j] <= min + (i+1)*bin_size) {
                    counter[i] = counter[i] + 1.;
            }

        }

        qMean[i] = min + bin_size/2 + i*bin_size;
    }

    double norm = ArraySum(counter, HISTBARS)*bin_size;

    for(int i = 0; i < HISTBARS; i++) {

        counter[i] = counter[i]/norm;
    }
    double skew = Skewness(qMean, counter, bin_size, HISTBARS);
    double errS = errSkewness(qMean, counter, bin_size, HISTBARS, skew);
    double kurt = Kurtosis(qMean, counter, bin_size, HISTBARS);
    double errK = errKurtosis(qMean, counter, bin_size, HISTBARS, kurt);
    double mean = mean_of_distribution(qMean, counter, bin_size, HISTBARS);
    double devStd = StandardDeviation(qMean, mean, HISTBARS);

    PrintPQ(qMean, counter, HISTBARS, temp, p, mean,devStd, skew, errS, kurt, errK, block);
    PlotPQScript(temp, p, T, mean,devStd, skew,errS, kurt, errK, block);

}
void HistogramIFO(double* vec, int size, int temp, double* T) {

    double* counter = (double*)calloc(HISTBARS, sizeof(double));
    checkDouble(counter);

    double* qMean = (double*)calloc(HISTBARS, sizeof(double));
    checkDouble(qMean);

    double min = -1.01;
    double max = 1.01;
    for(int i = 0; i < size; i++) {
      if(vec[i] > 1 || vec[i] < -1) {
        printf("OUT OF RANGE!");
      }
    }

/*
    for(int i = 0; i < size; i++) {

        if(min>vec[i]){
            min = vec[i];
        }
        if(max<vec[i]) {
            max = vec[i];
        }

    }
*/

    double bin_size = (max - min)/HISTBARS;

    for(int i = 0; i < HISTBARS; i++) {

        for(int j = 0; j < size; j++) {

            if(vec[j] > min + i*bin_size && vec[j] <= min + (i+1)*bin_size) {
                    counter[i] = counter[i] + 1.;
            }

        }

        qMean[i] = min + bin_size/2 + i*bin_size;
    }

    double norm = ArraySum(counter, HISTBARS)*bin_size;

    for(int i = 0; i < HISTBARS; i++) {

        counter[i] = counter[i]/norm;
    }

    PrintIFO(qMean, counter, HISTBARS, temp);
    PlotIFOScript(temp, T);

}




void HistogramExpIFO(double* vec, int size, int temp, double* T, int block, int c) {

    double* counter = (double*)calloc(HISTBARS, sizeof(double));
    checkDouble(counter);

    double* qMean = (double*)calloc(HISTBARS, sizeof(double));
    checkDouble(qMean);

    double min = -1.01;
    double max = 1.01;
    for(int i = 0; i < size; i++) {
      if(vec[i] > 1 || vec[i] < -1) {
        printf("OUT OF RANGE!");
      }
    }

/*
    for(int i = 0; i < size; i++) {

        if(min>vec[i]){
            min = vec[i];
        }
        if(max<vec[i]) {
            max = vec[i];
        }

    }
*/

    double bin_size = (max - min)/HISTBARS;

    for(int i = 0; i < HISTBARS; i++) {

        for(int j = 0; j < size; j++) {

            if(vec[j] > min + i*bin_size && vec[j] <= min + (i+1)*bin_size) {
                    counter[i] = counter[i] + 1.;
            }

        }

        qMean[i] = min + bin_size/2 + i*bin_size;
    }

    double norm = ArraySum(counter, HISTBARS)*bin_size;

    for(int i = 0; i < HISTBARS; i++) {

        counter[i] = counter[i]/norm;
    }

    PrintExpIFO(qMean, counter, HISTBARS, temp, block, c);
    PlotExpIFOScript(temp, T, block, c);

}


void HistogramG(double* vec, int size, int temp, double* T) {

    double* counter = (double*)calloc(HISTBARS, sizeof(double));
    checkDouble(counter);

    double* qMean = (double*)calloc(HISTBARS, sizeof(double));
    checkDouble(qMean);

    double min = -1.01;
    double max = 1.01;

    for(int i = 0; i < size; i++) {
      if(vec[i] > 1 || vec[i] < -1) {
        printf("OUT OF RANGE!");
      }
    }
/*
    for(int i = 0; i < size; i++) {

        if(min>vec[i]){
            min = vec[i];
        }
        if(max<vec[i]) {
            max = vec[i];
        }

    }
*/

    double bin_size = (max - min)/HISTBARS;

    for(int i = 0; i < HISTBARS; i++) {

        for(int j = 0; j < size; j++) {

            if(vec[j] > min + i*bin_size && vec[j] <= min + (i+1)*bin_size) {
                    counter[i] = counter[i] + 1.;
            }

        }

        qMean[i] = min + bin_size/2 + i*bin_size;
    }

    double norm = ArraySum(counter, HISTBARS)*bin_size;

    for(int i = 0; i < HISTBARS; i++) {

        counter[i] = counter[i]/norm;
    }

    PrintG(qMean, counter, HISTBARS, temp);
    PlotGScript(temp, T);

}

void HistogramEPS(double* vec, int size, int temp, double* T) {
	
	double min = -0.025;
    double max = 3;
    double* counter = (double*)calloc(max*HISTBARS, sizeof(double));
    checkDouble(counter);

    double* qMean = (double*)calloc(max*HISTBARS, sizeof(double));
    checkDouble(qMean);

/*
    double min = 0.;
    double max = 0.;
*/


    


    
/*
    for(int i = 0; i < size; i++) {

        if(min>vec[i]){
            min = vec[i];
        }
        if(max<vec[i]) {
            max = vec[i];
        }

    }
*/

    double bin_size = (max - min)/(max*HISTBARS);

    for(int i = 0; i < max*HISTBARS; i++) {

        for(int j = 0; j < size; j++) {

            if(vec[j] > min + i*bin_size && vec[j] <= min + (i+1)*bin_size) {
                    counter[i] = counter[i] + 1.;
            }

        }

        qMean[i] = (min + bin_size/2 + i*bin_size);
    }

    double norm = ArraySum(counter, max*HISTBARS) * bin_size;

    for(int i = 0; i < max*HISTBARS; i++) {

        counter[i] = counter[i]/norm;
    }

    PrintEPS(qMean, counter, max*HISTBARS, temp);
    PlotEPSScript(temp, T);

}
////////////////////SPECTRUM TOOLS////////////////////////////////////

 
void PrintIntensity(double* vec, double* sigma, double* freq, int r, int temp, int N) {
	char* fileName = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(fileName);
	sprintf(fileName, "intensity_nrep%d_temp%d.dat", r, temp);
  	
  	
	FILE* fp = fopen(fileName, "w");
  	
	for(int i = 0; i < N; i++) {
		fprintf(fp, "%le\t%le\t%le\n", freq[i], vec[i], sigma[i]);
	}

	fclose(fp);
	free(fileName);
}

double sumKnotJ(double* vec, int j, int size) {
  double sum = 0.;
  for(int k = 0; k < size; k++) {
    if(k != j) {
      sum+=vec[k];
    }
  }
  return sum;
}
double* sigmaIntensity(double* Intensity, double* A2M, double* sigma, int size) {

	double* sigmaIk = (double*)calloc(size, sizeof(double));
	checkDouble(sigmaIk);
  double** sigmaIj = (double**)calloc(size, sizeof(double*));
  double A2MSum = ArraySum(A2M, size);

  checkMat(sigmaIj);
  for(int k = 0; k < size; k++) {
    sigmaIj[k] = (double*)calloc(size, sizeof(double));
    checkDouble(sigmaIj[k]);
    for(int j = 0; j < size; j++) {
      sigmaIj[k][j] = (sigma[j]*Intensity[k]/A2MSum)*(sigma[j]*Intensity[k]/A2MSum);
    }
  }


  for(int k = 0; k < size; k++) {
    sigmaIk[k] = (sigma[k]*(1-Intensity[k])/A2MSum)*(sigma[k]*(1-Intensity[k])/A2MSum)+sumKnotJ(sigmaIj[k], k, size);
    sigmaIk[k] = sqrt(sigmaIk[k]);
  }

  for(int k = 0; k < size; k ++) {
  	free(sigmaIj[k]);
  }
  free(sigmaIj);
  return sigmaIk;

}



void build_spectrum(double* vec, double* sigma, int size, double* freq, int r, int k) {
	double sum = ArraySum(vec, size);
	double* I = ScaleVector(vec, size, sum);
	double* sigmaI = sigmaIntensity(I, vec, sigma, size);
	PrintIntensity(I, sigmaI, freq, r, k, size);
	free(sigmaI);
}

void checkConfiguration(double* vec, int size, int lim) {
  double s = ArraySum(vec, size);
  //printf("SUM: %.20lf\n", s);
  if(s < 1.*lim - 0.2 && s  > 1.*lim + 0.2) {
    printf("Constraint not respected! Exit.\n");
    exit(BAD_SIMULATION);
  }

}
///////////////// TERMALIZATION CHECK TOOLS/////////////////////////////////

void meanLogaritmicWindow(double* vec, int size, int r , int k, char type) {


	int temp = size;
	int n = 0;

	do{
			//printf("temp: %d\n", temp);
			temp = temp/2;
			n++;
	}while(temp > 1);


	//printf("size: %d\tn:%d\n", size, n);

	double** meanArray = (double**)calloc(2,sizeof(double*));
	checkMat(meanArray);
	meanArray[0] = (double*)calloc(n, sizeof(double));
	checkDouble(meanArray[0]);
	meanArray[1] = (double*)calloc(n, sizeof(double));
	checkDouble(meanArray[1]);

	int w = 1;

	char* fileName = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(fileName);

	sprintf(fileName, "check_tempering_%c_nrep_%d_temp%d.dat",type,  r, k);
	FILE* fw = fopen(fileName, "w");

  int min = 0;
  int max = 0;


	for(int i = 0; i < n; i++) {
    if(i == 0) {
      max = min + 1;
    }else{
      max = 2*min -1;
    }
    //printf("%d\t%d\t%d\n", min, max, w);

    w = max - min +1;
    double* temp_vec = (double*)calloc(w, sizeof(double));
    checkDouble(temp_vec);

    for(int j = 0; j < w; j++) {
      temp_vec[j] = vec[j+2*i];
    }
    meanArray[0][i] = MeanVec(temp_vec, w);
    meanArray[1][i] = StandardDeviation(temp_vec, meanArray[0][i], w);
    fprintf(fw, "%d\t%lf\t%lf\n", max*NSTEP, meanArray[0][i], meanArray[1][i]);
    min = max + 1;
    free(temp_vec);
	}

  free(meanArray);
	fclose(fw);


  /*
  meanArray[0][i] = MeanVec(temp_vec, w);
  meanArray[1][i] = StandardDeviation(temp_vec, meanArray[0][i], w);
  fprintf(fw, "%d\t%lf\t%lf\n", w, meanArray[0][i], meanArray[1][i]);
  */

}
////// SPECIFIC HEAT TOOLS ///////////////////////////////////////////////

double* EquilibriumEnergy(double* erg, int t_max) {
	//This function extrapolates the data of the energy at equilibrium
	int size = (int)(t_max/2);
	double* eq_erg = (double*)calloc(size, sizeof(double));
	checkDouble(eq_erg);

	for(int i = 0; i < size; i++) {
		eq_erg[i] = erg[i+size];
	}

	return eq_erg;

}

////// FREQUENCY MANIPULATIONS ////////////////////////////////////////////
double* GenCombLikeFrequencies(double min, double max, int size) {

    double* w = (double*)calloc(size, sizeof(double));
    checkDouble(w);

    double step = (max - min)/size;

    for(int i = 0; i < size; i ++) {
      w[i] = step + min + i*step;
      //printf("%g\n", w[i]);
    }

    return w;

}

long long int NumTetrads(double* w, int size, double Gamma) {
    long long int tetrads = 0;
    double a,b,c,d;
    int count = 0;
    //printf("FREQUENZE\n");
    char* interaction_file_name = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(interaction_file_name);
    sprintf(interaction_file_name, "interactions_file_size%d_gamma%g.dat", size, Gamma);
    FILE* graph_file = fopen(interaction_file_name, "w");
    checkFile(graph_file, interaction_file_name);
    free(interaction_file_name);
    //printf("\n");
    //printf("GAMMA:%g\n", Gamma);
    if(Gamma > 0) {
      for(int l = 0; l < size; l++) { // FOR L
        for(int k = 0; k < l; k++) { // FOR k
          for(int j = 0; j < k; j++) { //FOR J
            for(int i = 0; i < j; i++) { //FOR I
              a = w[i];
              b = w[l];
              c = w[k];
              d = w[j];
              //printf("|%g + %g - %g - %g| = %g\n", a,b,c,d, fabs(a+b-c-d));
              count++;
              if(fabs(a+b-c-d) <= Gamma) {
                //printf("N\tGamma\tFMC\n");
                //printf("%g+%g-%g-%g\n", w[i], w[l], w[k], w[j]);
                //printf("\n\nFMC:%d\t%lf\t%lf\n\n\n", size, Gamma, abs(w[i]+w[l]-w[k]-w[j]));
                fprintf(graph_file, "%d\t%d\t%d\t%d\n", l, k, j, i);
                tetrads++;
              }else if(fabs(a+c-b-d) <= Gamma) {
                tetrads++;
                fprintf(graph_file, "%d\t%d\t%d\t%d\n", l, k, j, i);

              }else if(fabs(a+d-b-c) <= Gamma) {
                tetrads++;
                fprintf(graph_file, "%d\t%d\t%d\t%d\n", l, k, j, i);

              }//FINE IF CHECK FMC
            } // FINE FOR I
          }//FINE FOR J
        }//FINE FOR k
      }//FINE FOR L
    }else{
      for(int l = 0; l < size; l++) { // FOR L
        for(int k = 0; k < l; k++) { // FOR k
          for(int j = 0; j < k; j++) { //FOR J
            for(int i = 0; i < j; i++) { //FOR I
              a = w[i];
              b = w[l];
              c = w[k];
              d = w[j];
              //printf("|%g + %g - %g - %g| = %g\n", a,b,c,d, fabs(a+b-c-d));
              count++;
              if(fabs(a+b-c-d) <= ZERO) {
                //printf("N\tGamma\tFMC\n");
                //printf("%g+%g-%g-%g\n", w[i], w[l], w[k], w[j]);
                //printf("\n\nFMC:%d\t%lf\t%lf\n\n\n", size, Gamma, abs(w[i]+w[l]-w[k]-w[j]));
                tetrads++;
                fprintf(graph_file, "%d\t%d\t%d\t%d\n", l, k, j, i);

              }else if(fabs(a+c-b-d) <= ZERO) {
                tetrads++;
                fprintf(graph_file, "%d\t%d\t%d\t%d\n", l, k, j, i);

              }else if(fabs(a+d-b-c) <= ZERO) {
                tetrads++;//FINE IF CHECK FMC
                fprintf(graph_file, "%d\t%d\t%d\t%d\n", l, k, j, i);

              }
            }// FINE FOR I
          }//FINE FOR J
        }//FINE FOR k
      }//FINE FOR L
    }

    fclose(graph_file);
    //system("mkdir stuffs");
    //system("mv interactions_file_gtc.dat stuffs");
    //printf("COUNT:%d\n\n", count);
    //printf("N4:%lld\n\n", tetrads);
    //printf("|%g+%g-%g-%g|=%g\n", a,b,c,d, fabs(a+b-c-d));
    return tetrads;
}
long long int* AllNumTetrads(double* w, int size, double Gamma) {
    long long int* tetrads = (long long int*)calloc(3, sizeof(long long int));

    double a,b,c,d;
    int count = 0;
    //printf("FREQUENZE\n");

    for(int i = 0; i < size; i++) {
      //printf("%g\n", w[i]);
    }
    //printf("\n");
    //printf("GAMMA:%g\n", Gamma);
    if(Gamma > 0) {
      for(int l = 0; l < size; l++) { // FOR L
        for(int k = 0; k < l; k++) { // FOR k
          for(int j = 0; j < k; j++) { //FOR J
            for(int i = 0; i < j; i++) { //FOR I
              a = w[i];
              b = w[l];
              c = w[k];
              d = w[j];
              //printf("|%g + %g - %g - %g| = %g\n", a,b,c,d, fabs(a+b-c-d));
              count++;
              if(fabs(a+b-c-d) <= Gamma) { //primo ordinamento
                //printf("N\tGamma\tFMC\n");
                //printf("%g+%g-%g-%g\n", w[i], w[l], w[k], w[j]);
                //printf("\n\nFMC:%d\t%lf\t%lf\n\n\n", size, Gamma, abs(w[i]+w[l]-w[k]-w[j]));
                tetrads[0]++;
              }else if(fabs(a+c-b-d) <= Gamma) {
                tetrads[1]++;
              }else if(fabs(a+d-b-c) <= Gamma) {
                tetrads[2]++;
              }//FINE IF CHECK FMC
            } // FINE FOR I
          }//FINE FOR J
        }//FINE FOR k
      }//FINE FOR L
    }else{
      for(int l = 0; l < size; l++) { // FOR L
        for(int k = 0; k < l; k++) { // FOR k
          for(int j = 0; j < k; j++) { //FOR J
            for(int i = 0; i < j; i++) { //FOR I
              a = w[i];
              b = w[l];
              c = w[k];
              d = w[j];
              //printf("|%g + %g - %g - %g| = %g\n", a,b,c,d, fabs(a+b-c-d));
              count++;
              if(fabs(a+b-c-d) <= ZERO) { //ZERO is defined in gtcStructures.h
                //printf("N\tGamma\tFMC\n");
                //printf("%g+%g-%g-%g\n", w[i], w[l], w[k], w[j]);
                //printf("\n\nFMC:%d\t%lf\t%lf\n\n\n", size, Gamma, abs(w[i]+w[l]-w[k]-w[j]));
                tetrads[0]++;
              }else if(fabs(a+c-b-d) <= ZERO) {
                tetrads[1]++;
              }else if(fabs(a+d-b-c) <= ZERO) {
                tetrads[2]++;
              }//FINE IF CHECK FMC
            }// FINE FOR I
          }//FINE FOR J
        }//FINE FOR k
      }//FINE FOR L
    }
    //printf("COUNT:%d\n\n", count);
    //printf("N4:%lld\n\n", tetrads);
    //printf("|%g+%g-%g-%g|=%g\n", a,b,c,d, fabs(a+b-c-d));
    return tetrads;
}


long long int NearPowTwo(int n) {
  long long int p = 1;

  while(p <= n) {
          p = p*2;
  }

   return (int)p/2;
}

int NumTetradsFromFILE(int size, int Gamma) {

    int tetrads = 0;
    FILE* fr = fopen("frequencies.dat", "r");

    double* w1 = (double*)calloc(size, sizeof(double));
    checkDouble(w1);
    double* w2 = (double*)calloc(size, sizeof(double));
    checkDouble(w2);
    double* w3 = (double*)calloc(size, sizeof(double));
    checkDouble(w3);
    double* w4 = (double*)calloc(size, sizeof(double));
    checkDouble(w4);


    for(int j = 0; j < size; j++) {
        fscanf(fr, "%le\n", &w1[j]);
        //printf("%lf\n", w[j]);
        w2[j] = w1[j];
        w3[j] = w1[j];
        w4[j] = w1[j];

    }

    fclose(fr);

    //COMPUTING NPLAQ

    for(int i = 0; i < size; i++) {
        for(int j = 0; j < i; j++) {
            for(int k = 0; k < j; k++) {
                for(int l = 0; l < k; l++) {

                    if(i != j && i != k && i != l && j != k && j != l && k != l) {

                        if(fabs(w1[i]-w2[j]+w3[k]-w4[l]) <= Gamma) {
                            tetrads++;
                        }
                    }
                }
            }
        }
    }

    return tetrads;

}

double** ReadInteractionFile(int Nplaq, int choose_mod) {

    FILE* fr;
    char* tempfile_name = (char*)calloc(STRING_SIZE, sizeof(char*));
    checkString(tempfile_name);
    if(choose_mod == 0) {
      sprintf(tempfile_name, "interactions_file.dat");
    }else if(choose_mod == 1) {
      sprintf(tempfile_name, "stuffs/interactions_file.dat");
    }
    fr = fopen(tempfile_name, "r");
    checkFile(fr, tempfile_name);
    free(tempfile_name);
    double** data = (double**)calloc(Nplaq, sizeof(double*));
    checkMat(data);
    //double J = 0.;
    for(int i = 0; i < Nplaq; i++) {

        data[i] = (double*)calloc(5, sizeof(double));
        checkDouble(data[i]);

        fscanf(fr, "%lf %lf %lf %lf %le\n", &data[i][0], &data[i][1], &data[i][2], &data[i][3], &data[i][4]);


    }

    fclose(fr);
    return data;
}

double* GenUniformFrequencies(double min, double max, int size) {

    double* w = (double*)calloc(size, sizeof(double));
    checkDouble(w);

    for(int k = 0; k < size; k++) {
        w[k] = min + (max - min)*drand48();
        //printf("%lf\n", w[k]);
    }
    qsort(w, size, sizeof(double), compare);
    return w;
}


void printFreq(double* w, int size) {

    char* fileName = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(fileName);

    sprintf(fileName, "frequencies_N%d.dat", size);

    FILE* fw = fopen(fileName, "w");
    checkFile(fw, fileName);

    free(fileName);
    for(int k = 0; k < size; k++) {

        fprintf(fw, "%le\n", w[k]);
    }

    fclose(fw);

}

//double* freqOnGraph()


////// OVERLAP TOOLS //////////////////////////////////////////////////////
complex** OverlapMatrix(complex** mat, int nRep, int N) {

	// mat[replicas][modes]
	complex** OVMatrix = (complex**)calloc(nRep, sizeof(complex*));
	checkComplexMat(OVMatrix);
	//allocation
	for(int r = 0; r < nRep; r++) {
		OVMatrix[r] = (complex*)calloc(nRep, sizeof(complex));
		checkComplex(OVMatrix[r]);
	}

	for(int alpha = 0; alpha < nRep; alpha++) {
		for(int gamma = 0; gamma < nRep; gamma++) {
			for(int i = 0; i < N; i++) {
				OVMatrix[alpha][gamma].re = (OVMatrix[alpha][gamma].re + ComplexProd(Conjugate(mat[alpha][i]), mat[gamma][i]).re);
				OVMatrix[alpha][gamma].im = (OVMatrix[alpha][gamma].im + ComplexProd(Conjugate(mat[alpha][i]), mat[gamma][i]).im);
			}
		}
	}


	return OVMatrix;
}
double* ComplexSupTriang(complex** q, int size, int p) {
    //p = 0 -> real part, p = 1 - > imaginary part
	double* vec = (double*)calloc((int)(size*(size-1)/2), sizeof(double));
	checkDouble(vec);
    if(p != 0 && p != 1) {
		printf("Choose p = 0 for real part or p = 1 for imaginary part!\n");
		exit(SYNTAX_ERROR);
    }
	int l = 0;
	for(int alpha = 0; alpha < size; alpha++) {
		for(int gamma = 0; gamma < size; gamma++) {
			if(gamma>alpha) {
				if(p == 0) vec[l] = q[alpha][gamma].re;
                if(p == 1) vec[l] = q[alpha][gamma].im;
				l++;
			}
		}
	}

 return vec;
}
double* SupTriang(double** q, int size) {

	double* vec = (double*)calloc((int)(size*(size+1)/2), sizeof(double));
	checkDouble(vec);
	int l = 0;
	for(int alpha = 0; alpha < size; alpha++) {
		for(int gamma = 0; gamma < size; gamma++) {
			if(gamma>=alpha) {
				vec[l] = q[alpha][gamma];
				l++;
			}
		}
	}

 return vec;
}
double OverlapIFO_Element(double* Ir1, double* Ir2, double* mr1, double* mr2, int size) {
	
	double ifo = 0.;
	double sum1 = 0.;
	double sum2 = 0.;
	
	for(int j = 0; j < size; j++) {
		ifo += (Ir1[j] - mr1[j])*(Ir2[j] - mr2[j]);
		sum1 += (Ir1[j] - mr1[j])*(Ir1[j] - mr1[j]);
		sum2 += (Ir2[j] - mr2[j])*(Ir2[j] - mr2[j]);
	}
	
	return ifo / (sqrt(sum1)*sqrt(sum2));
   
   
}

double MixedOverlap_Element(complex* vec1, complex* vec2, int size) {

    //mat[replica][mode]
	double Element = 0.;

    double* sigma = (double*)calloc(size, sizeof(double));
    checkDouble(sigma);

    double* tau = (double*)calloc(size, sizeof(double));
    checkDouble(tau);
	for(int j = 0; j < size; j++) {

        sigma[j] = vec1[j].re;
        tau[j] = vec2[j].im;

		Element += sigma[j]*tau[j];
	}

	free(sigma);
    free(tau);

	return Element/size;
}
double Plaq(complex* a, double* index) {  ///STIAMO USANDO LA CONVENZIONE A*_I A_J A_k A*_L

    double e = 0.;

    //e = (a[(int)index[0]].re * a[(int)index[3]].re - a[(int)index[0]].im * a[(int)index[3]].im)*(a[(int)index[2]].re * a[(int)index[1]].re - a[(int)index[2]].im * a[(int)index[1]].im) + (a[(int)index[0]].im * a[(int)index[3]].re + a[(int)index[0]].re * a[(int)index[3]].im ) * (a[(int)index[2]].im * a[(int)index[1]].re + a[(int)index[2]].re * a[(int)index[1]].im );
//    e = (a[(int)index[0]].re * a[(int)index[1]].re - a[(int)index[0]].im * a[(int)index[1]].im)*(a[(int)index[2]].re * a[(int)index[3]].re - a[(int)index[2]].im * a[(int)index[3]].im) + (a[(int)index[0]].im * a[(int)index[1]].re + a[(int)index[0]].re * a[(int)index[1]].im ) * (a[(int)index[2]].im * a[(int)index[3]].re + a[(int)index[2]].re * a[(int)index[3]].im );
     //Ai1 Ai2 Ai3 Ai4 cos φi2 − φi1 + φi3 − φi4 ,
     //e = ModValue(a[(int)index[0]])*ModValue(a[(int)index[3]])*ModValue(a[(int)index[2]])*ModValue(a[(int)index[1]])*cos(ArgValue(a[(int)index[3]]) - ArgValue(a[(int)index[0]]) + ArgValue(a[(int)index[2]]) -  ArgValue(a[(int)index[1]]));

     //e = (complex a, complex b)ComplexProd(complex a, complex b)

     complex ai = a[(int)index[0]];
     complex aj = a[(int)index[3]];
     complex ak = a[(int)index[2]];
     complex al = a[(int)index[1]];

     //double J = index[4];

     double ti = ai.im;
     double tj = aj.im;
     double tk = ak.im;
     double tl = al.im;
     double si = ai.re;
     double sj = aj.re;
     double sk = ak.re;
     double sl = al.re;
     e = ((sl*si-tl*ti)*sj+(tl*si+sl*ti)*tj)*sk + ((tl*ti-sl*si)*tj+(sl*ti+tl*si)*sj)*tk;

    return e;

}

// OTHER TOOLS

void printAllSpectra(double**** checkA, double**** checkPhi,double* freq, int* iter, int nFILE, int nPT, int nRep, int N) {
  FILE* instant_spectra_file;
	char* instant_spectra_name;

	printf("PLOTTING INSTANT SPECTRUM....\n");
	int counter = 0;
	int tot_cycle = N*nFILE*nPT*nRep;
	for(int r = 0; r < nRep; r++) {
		for(int k = 0; k < nPT; k++) {
			if(NITER_MIN_PRINT == 0) {
        for(int i = 0; i < (int)(nFILE/2) +1; i++) {
				  instant_spectra_name = (char*)calloc(STRING_SIZE, sizeof(char));
				  sprintf(instant_spectra_name, "instant_spectra_nrep%d_temp%d_iter%d.dat", r, k, iter[i+(int)(nFILE/2)]);
				  checkString(instant_spectra_name);
				  instant_spectra_file = fopen(instant_spectra_name, "w");
				  checkFile(instant_spectra_file, instant_spectra_name);
				  for(int j = 0; j < N; j++) {
					  fprintf(instant_spectra_file, "%le\t%le\t%le\n", freq[j], checkA[r][i][k][j]/ArraySum(checkA[r][i][k], N), checkPhi[r][i][k][j]);
					  counter++;
					  printf("%.1lf%%\r", ((double)counter/tot_cycle)*100);
				  }

				  fclose(instant_spectra_file);
				  plotInstantSpectra(r, k, iter[i+(int)(nFILE/2)]);
				  plotInstantPhase(r, k, iter[i+(int)(nFILE/2)]);
			  }
      }else{
        for(int i = 0; i < nFILE; i++) {
				  instant_spectra_name = (char*)calloc(STRING_SIZE, sizeof(char));
				  sprintf(instant_spectra_name, "instant_spectra_nrep%d_temp%d_iter%d.dat", r, k, iter[i]);
				  checkString(instant_spectra_name);
				  instant_spectra_file = fopen(instant_spectra_name, "w");
				  checkFile(instant_spectra_file, instant_spectra_name);
				  for(int j = 0; j < N; j++) {
					  fprintf(instant_spectra_file, "%le\t%le\t%le\n", freq[j], checkA[r][i][k][j]/ArraySum(checkA[r][i][k], N), checkPhi[r][i][k][j]);
					  counter++;
					  printf("%.1lf%%\r", ((double)counter/tot_cycle)*100);
				  }

				  fclose(instant_spectra_file);
				  //plotInstantSpectra(r, k, iter[i]);
				  //plotInstantPhase(r, k, iter[i]);
			  }
      }

		}
	}
	free(instant_spectra_name);
}





complex**** fix_dynamics(complex*** EqDynamics, int nRep, int nPT, int N, int nFILE) {
  complex**** deblockedDynamics = (complex****)calloc(nRep, sizeof(complex***));
	checkHyp2ComplexMat(deblockedDynamics);
	//printf("TUTTO OK\n");
 //ComplexDeBlock(complex* vec, int BlockSize, int nBlock, int k)
    int counter = 0;
    int tot_cycle = nRep*nFILE*nPT;
    printf("Rearranging arrays of dynamics...\n");
	for(int r = 0; r < nRep; r++) {
		deblockedDynamics[r] = (complex***)calloc(nFILE, sizeof(complex**));
		checkHypComplexMat(deblockedDynamics[r]);
		for(int i = 0; i < nFILE; i++) {
			deblockedDynamics[r][i] = (complex**)calloc(nPT, sizeof(complex*));
			checkComplexMat(deblockedDynamics[r][i]);
			for(int k = 0; k < nPT; k++) {
				deblockedDynamics[r][i][k] = ComplexDeBlock(EqDynamics[r][i], N, nPT, k);
				counter++;
				printf("%.1lf%%\r", ((double)counter/tot_cycle)*100);
			}
		}
	}
	//PrintComplexMat(deblockedDynamics[0][0], 2,2);
	printf("\n");
	//printf("TUTTO OK\n");
	//SPENDO QUALCHE SECONDO DI COMPUTAZIONE MA OTTENGO UNA MATRICE SENSATA
	complex**** FixedDynamics = (complex****)calloc(nPT, sizeof(complex***));
	checkHyp2ComplexMat(FixedDynamics);
	printf("I allow dynamics indices to calculate overlaps...\n");
	counter = 0;
    tot_cycle = nRep*nFILE*nPT*N;
	for(int k = 0; k < nPT; k++) {
		FixedDynamics[k] = (complex***)calloc(nFILE, sizeof(complex**));
		checkHypComplexMat(FixedDynamics[k]);
		for(int i = 0; i < nFILE; i++) {
			FixedDynamics[k][i] = (complex**)calloc(nRep, sizeof(complex*));
			checkComplexMat(FixedDynamics[k][i]);
			for(int r = 0; r < nRep; r++) {
				FixedDynamics[k][i][r] = (complex*)calloc(N, sizeof(complex));
				checkComplex(FixedDynamics[k][i][r]);
				for(int j = 0; j < N; j++) {
					FixedDynamics[k][i][r][j] = deblockedDynamics[r][i][k][j];
					counter++;
					printf("%.1lf%%\r", ((double)counter/tot_cycle)*100);
				}
			}
		}
	}
	free(deblockedDynamics);

  return FixedDynamics;
}

void read_equilibrium_configurations(complex*** EqDynamics, double*** A2, double*** phiDyn, int nRep, int nFILE, int nPT, int N, int* iter, int choose_mod) {


  for(int r = 0; r < nRep; r++) {
    if(NITER_MIN_PRINT == 0) {
      A2[r] = (double**)calloc((int)(nFILE/2) +1, sizeof(double*));
      checkMat(A2[r]);
      phiDyn[r] = (double**)calloc((int)(nFILE/2) +1, sizeof(double*));
      checkMat(phiDyn[r]);

      EqDynamics[r] = (complex**)calloc((int)(nFILE/2) +1, sizeof(complex));
      checkComplexMat(EqDynamics[r]);
      for(int i = 0; i < (int)(nFILE/2) +1; i++) {

        //temp_spin = readConfigFILE(r, iter[i], N*nPT);
        EqDynamics[r][i] = readConfigFILE(r, iter[i+(int)(nFILE/2)], N*nPT, choose_mod);
        A2[r][i] = ComplexVec2(EqDynamics[r][i], N*nPT);
        phiDyn[r][i] = ComplexVecPhi(EqDynamics[r][i], N*nPT);
        //free(temp_spin);
      //  counter++;
        //printf("%.1lf%%\r", ((double)counter/tot_cycle/2)*100);
      }
    }else{
      A2[r] = (double**)calloc(nFILE, sizeof(double*));
      checkMat(A2[r]);
      phiDyn[r] = (double**)calloc(nFILE, sizeof(double*));
      checkMat(phiDyn[r]);

      EqDynamics[r] = (complex**)calloc(nFILE, sizeof(complex));
      checkComplexMat(EqDynamics[r]);
      for(int i = 0; i < nFILE; i++) {

      //temp_spin = readConfigFILE(r, iter[i], N*nPT);
        EqDynamics[r][i] = readConfigFILE(r, iter[i], N*nPT, choose_mod);
        A2[r][i] = ComplexVec2(EqDynamics[r][i], N*nPT);
        phiDyn[r][i] = ComplexVecPhi(EqDynamics[r][i], N*nPT);
        //free(temp_spin);
      //  counter++;
        //printf("%.1lf%%\r", ((double)counter/tot_cycle)*100);
      }
    }
  }


}

void read_all_configurations(complex*** Dynamics, int nRep, int nFILE, int nPT, int N, int* iter, int choose_mod) {

  for(int r = 0; r < nRep; r++) {

    Dynamics[r] = (complex**)calloc(nFILE, sizeof(complex));
    checkComplexMat(Dynamics[r]);
    for(int i = 0; i < nFILE; i++) {

    //temp_spin = readConfigFILE(r, iter[i], N*nPT);
      Dynamics[r][i] = readConfigFILE(r, iter[i], N*nPT, choose_mod);
    }
  }
}

void compute_specific_heat(double**** wholeDynamics,int nRep, int nPT, int N, int dITER, int IT_MAX, double* temp) {



  double*** eqErg = (double***)calloc(nRep, sizeof(double**));
  checkHypMat(eqErg);

  double* meanErg = (double*)calloc(nPT-1, sizeof(double));
  checkDouble(meanErg);


  //ERG2
  /*
  double*** eqErg2 = (double***)calloc(nRep, sizeof(double**));
  checkHypMat(eqErg2);

  double** meanErg2 = (double**)calloc(2, sizeof(double*));
  checkMat(meanErg2);

  meanErg2[0] = (double*)calloc(nPT-1, sizeof(double));
  checkDouble(meanErg2[0]);

  meanErg2[1] = (double*)calloc(nPT-1, sizeof(double));
  checkDouble(meanErg2[1]);
*/
  //specific HEAT

  //double** sHeat = (double**)calloc(2, sizeof(double*));

/////////////////////////////////////////////////

  for(int r = 0; r < nRep; r++) { //CICLO REPLICHE

      eqErg[r] = (double**)calloc(nPT-1, sizeof(double*));
      checkMat(eqErg[r]);
      /*
      eqErg2[r] = (double**)calloc(nPT-1, sizeof(double*));
      checkMat(eqErg2[r]);
      */
      double* sh = (double*)calloc(nPT-1, sizeof(double));
      checkDouble(sh);
      double* errSh = (double*)calloc(nPT-1, sizeof(double));

      for(int k = 0; k < nPT-1; k++) { //CICLO TEMPERATURE
          //calcolo media e errori SOLO sui dati all'equilibrio, sennò
          //non vale il teorema di fluttuazione-dissipazione :)


          eqErg[r][k] = EquilibriumEnergy(wholeDynamics[r][1][k], IT_MAX+dITER);
          //eqErg2[r][k] = sqVec(eqErg[r][k], (int)((IT_MAX+dITER)/2));

          meanErg[k] = MeanVec(eqErg[r][k], (int)((IT_MAX+dITER)/2));

          for(int i = 0; i < (int)((IT_MAX+dITER)/2); i++) {
            sh[k] += (eqErg[r][k][i] - meanErg[k])*(eqErg[r][k][i] - meanErg[k]) / (int)((IT_MAX+dITER)/2);
          }

          for(int i = 0; i < (int)((IT_MAX+dITER)/2); i++) {
            errSh[k] += ((eqErg[r][k][i] - meanErg[k])*(eqErg[r][k][i] - meanErg[k]) - sh[k])*((eqErg[r][k][i] - meanErg[k])*(eqErg[r][k][i] - meanErg[k]) - sh[k])/ (int)((IT_MAX+dITER)/2);
          }

          sh[k] = sh[k] / (N * temp[k] * temp[k]);
          errSh[k] = sqrt((errSh[k] / (N * temp[k] * temp[k]))/(int)((IT_MAX+dITER)/2));
          //meanErg[k] = StandardDeviation(eqErg[r][k], meanErg[0][k],(int)((IT_MAX+dITER)/2));
          /*
          meanErg2[0][k] = MeanVec(eqErg2[r][k], (int)((IT_MAX+dITER)/2));
          meanErg2[1][k] = StandardDeviation(eqErg2[r][k], meanErg2[0][k],(int)((IT_MAX+dITER)/2));
          */

      }//FINE CICLO TEMPERATURE

      //qui calcolo il calore specifico

      //sHeat = SpecificHeat(meanErg, meanErg2, nPT-1, N, temp);

      //void PrintSpecificHeat(double** mat, int size,int r,double* temp, int N)
      PrintSpecificHeat(sh, errSh, nPT-1, r, temp, N);
      PlotSpecificHeatScript(r, N);

      free(sh);
      free(errSh);

  }//FINE CICLO REPLICHE

  free(eqErg);
  free(meanErg);

}



void check_spherical_constrait(double**** checkA, double**** checkPhi, double*** A2, double*** phiDyn, int nRep, int nFILE, int nPT, int N) {
  if(NITER_MIN_PRINT == 0) {

		for(int r = 0; r < nRep; r++) {
			checkA[r] = (double***)calloc((int)(nFILE/2) +1, sizeof(double**));
			checkHypMat(checkA[r]);
			checkPhi[r] = (double***)calloc((int)(nFILE/2) +1, sizeof(double**));
			checkHypMat(checkPhi[r]);
			for(int i = 0; i < (int)(nFILE/2) +1; i++) {

				checkA[r][i] = (double**)calloc(nPT, sizeof(double*));
				checkMat(checkA[r][i]);
				checkPhi[r][i] = (double**)calloc(nPT, sizeof(double*));
				checkMat(checkPhi[r][i]);
				for(int k = 0; k < nPT; k++) {
					checkA[r][i][k] = DeBlock(A2[r][i], N, nPT, k);
					checkPhi[r][i][k] = DeBlock(phiDyn[r][i], N, nPT, k);
					checkConfiguration(checkA[r][i][k], N, N);
				}
			}
		}
	}else{
		for(int r = 0; r < nRep; r++) {
			checkA[r] = (double***)calloc(nFILE, sizeof(double**));
			checkHypMat(checkA[r]);
			checkPhi[r] = (double***)calloc(nFILE, sizeof(double**));
			checkHypMat(checkPhi[r]);
			for(int i = 0; i < nFILE; i++) {

				checkA[r][i] = (double**)calloc(nPT, sizeof(double*));
				checkMat(checkA[r][i]);
				checkPhi[r][i] = (double**)calloc(nPT, sizeof(double*));
				checkMat(checkPhi[r][i]);
				for(int k = 0; k < nPT; k++) {
					checkA[r][i][k] = DeBlock(A2[r][i], N, nPT, k);
					checkPhi[r][i][k] = DeBlock(phiDyn[r][i], N, nPT, k);
					checkConfiguration(checkA[r][i][k], N, N);
				}
			}
		}
	}
}


complex**** compute_parisi_overlaps(complex**** FixedDynamics, int nRep, int nFILE, int nPT, int N) {

	printf("Computing overlaps matrix...\n");
	long long int counter = 0;
    long long int tot_cycle = nRep*nRep*nFILE*nPT;
	complex**** q = (complex****)calloc(nPT, sizeof(complex***));
	checkHyp2ComplexMat(q);
	for(int k = 0; k < nPT; k++) {
		q[k] = (complex***)calloc(nFILE, sizeof(complex**));
		checkHypComplexMat(q[k]);
   	    char* name_overlaps = (char*)calloc(STRING_SIZE, sizeof(char));
        checkString(name_overlaps);
        sprintf(name_overlaps, "parisi_overlaps_temp%d.dat", k);
        FILE* overlaps_file = fopen(name_overlaps, "w");
        checkFile(overlaps_file, name_overlaps);
        free(name_overlaps);
        fprintf(overlaps_file, "#n_iter\tr1\tr2\tRe(q)\tIm(q)\n");
		for(int i = 0; i < nFILE; i++) {
			q[k][i] = OverlapMatrix(FixedDynamics[k][i], nRep, N);
			for(int r1 = 0; r1<nRep; r1++) {
				for(int r2 = r1+1; r2<nRep; r2++) {
					q[k][i][r1][r2].re = q[k][i][r1][r2].re/N;
					q[k][i][r1][r2].im = q[k][i][r1][r2].im/N;
         		    if(r2>r1) {
                    	fprintf(overlaps_file, "%d\t%d\t%d\t%le\t%le\n", i, r1, r2, q[k][i][r1][r2].re, q[k][i][r1][r2].im);
                    }
					counter++;
					printf("%.1lf%%\r", ((double)counter/tot_cycle)*100);
				}
			}
            fprintf(overlaps_file, "\n");
		}
        fclose(overlaps_file);
	}
    return q;
}


double** compute_equilibrium_pq_real_part(complex**** q, int nRep, int nPT, int nFILE) {

 
  double* qTempRe = (double*)calloc((int)(nRep*(nRep-1)/2), sizeof(double));
  checkDouble(qTempRe);
  double** qMatRe = (double**)calloc(nPT, sizeof(double*));
  checkMat(qMatRe);
  //int l = 0;
  for(int k = 0; k < nPT; k++) {
    qMatRe[k] = (double*)calloc((int)(nRep*(nRep-1)/2)*nFILE, sizeof(double));
    checkDouble(qMatRe[k]);
    for(int i = 0; i < nFILE; i++) {
      qTempRe = ComplexSupTriang(q[k][i], nRep, 0); //0 sta per "parte reale"
      ReBlock(qTempRe, qMatRe[k],(int)(nRep*(nRep-1)/2), i);
    }
  }

  free(qTempRe);

  return qMatRe;


}



double** compute_equilibrium_pq_imaginary_part(complex**** q, int nRep, int nPT, int nFILE) {

  double* qTempIm = (double*)calloc((int)(nRep*(nRep-1)/2), sizeof(double));
  checkDouble(qTempIm);
  double** qMatIm = (double**)calloc(nPT, sizeof(double*));
  checkMat(qMatIm);
  //int l = 0;
  for(int k = 0; k < nPT; k++) {
    qMatIm[k] = (double*)calloc((int)(nRep*(nRep-1)/2)*nFILE, sizeof(double));
    checkDouble(qMatIm[k]);
    for(int i = 0; i < nFILE; i++) {
      qTempIm = ComplexSupTriang(q[k][i], nRep, 1); //1 sta per "parte immaginaria"
      ReBlock(qTempIm, qMatIm[k],(int)(nRep*(nRep-1)/2), i);
    }
  }

  free(qTempIm);


  double** qMatImFixed = (double**)calloc(nPT, sizeof(double*));
  checkMat(qMatImFixed);
  for(int k = 0; k < nPT; k++) {
    qMatImFixed[k] = (double*)calloc(nRep*(nRep-1)*nFILE, sizeof(double));
    checkDouble(qMatImFixed[k]);
    for(int l = 0; l < nRep*(nRep-1)*nFILE; l++) {
      if(l < (int)(nFILE*nRep*(nRep-1)/2)) {
        qMatImFixed[k][l] = qMatIm[k][l];
      }else {
        qMatImFixed[k][l] = -qMatIm[k][(int)(nRep*(nRep-1)/2)*nFILE - l];
      }
    }
  }
  free(qMatIm);

  return qMatImFixed;

}

complex**** overlaps_dynamics(complex**** q, int nPT, int nfile_temp, int nRep,int b) {

  complex**** q_temp = (complex****)calloc(nPT, sizeof(complex***));
  checkHyp2ComplexMat(q_temp);
  for(int k = 0; k < nPT; k++) {
    q_temp[k] = (complex***)calloc(nfile_temp, sizeof(complex**));
    checkHypComplexMat(q_temp[k]);
    for(int i = 0; i < nfile_temp; i++) {
      q_temp[k][i] = (complex**)calloc(nRep, sizeof(complex*));
      checkComplexMat(q_temp[k][i]);
      for(int r1 = 0; r1 < nRep; r1++) {
        q_temp[k][i][r1] = (complex*)calloc(nRep, sizeof(complex));
        checkComplex(q_temp[k][i][r1]);
        for(int r2 = 0; r2 < nRep; r2++) {
          q_temp[k][i][r1][r2].re = q[k][i+2*b][r1][r2].re;
          q_temp[k][i][r1][r2].im = q[k][i+2*b][r1][r2].im;
        }
      }
    }
  }

  return q_temp;
}

double ** compute_mixed_overlaps (complex**** FixedDynamics, int nRep, int nFILE, int nPT, int N) {

  double**** G = (double****)calloc(nPT, sizeof(double***));
  checkHyp2Mat(G);

  double** GVec = (double**)calloc(nPT, sizeof(double*));
  checkMat(GVec);
  int l = 0;

  for(int k = 0; k < nPT; k++) {
      G[k] = (double***)calloc(nFILE, sizeof(double**));
      checkHypMat(G[k]);

      GVec[k] = (double*)calloc(nRep*nRep*nFILE, sizeof(double));
      checkDouble(GVec[k]);

      l = 0;

      for(int i = 0; i < nFILE; i++) {
          G[k][i] = (double**)calloc(nRep, sizeof(double*));
          checkMat(G[k][i]);

          for(int r1 = 0; r1 < nRep; r1++) {
              G[k][i][r1] = (double*)calloc(nRep, sizeof(double));
              checkDouble(G[k][i][r1]);

              for(int r2 = 0; r2 < nRep; r2++) {
                  G[k][i][r1][r2] = MixedOverlap_Element(FixedDynamics[k][i][r1], FixedDynamics[k][i][r2], N);

                  GVec[k][l] = G[k][i][r1][r2];
                  l++;

              }
          }
      }
  }

  return GVec;

}



double** compute_intensity_fluctuations_overlaps(double*** A2,double** A2M, int nRep, int nFILE, int nPT, int N) {

  double**** A2_deBlocked = (double****)calloc(nRep, sizeof(double***));
  checkHyp2Mat(A2_deBlocked);
  int l;
  for(int r = 0; r < nRep; r++) {

    if(NITER_MIN_PRINT == 0) {
      A2_deBlocked[r] = (double***)calloc((int)(nFILE/2) +1, sizeof(double**));
      checkHypMat(A2_deBlocked[r]);
      for(int i = 0; i < (int)(nFILE/2) +1; i++) {
        A2_deBlocked[r][i] = (double**)calloc(nPT, sizeof(double*));
        checkMat(A2_deBlocked[r][i]);

        for(int k = 0; k < nPT; k++) {

          A2_deBlocked[r][i][k] = DeBlock(A2[r][i], N, nPT, k);
        //counter++;
        //printf("%.1lf%%\r", ((double)counter/tot_cycle)*100);
        }
      }
    }else{
      A2_deBlocked[r] = (double***)calloc(nFILE, sizeof(double**));
      checkHypMat(A2_deBlocked[r]);
      for(int i = 0; i < nFILE; i++) {
        A2_deBlocked[r][i] = (double**)calloc(nPT, sizeof(double*));
        checkMat(A2_deBlocked[r][i]);

        for(int k = 0; k < nPT; k++) {

          A2_deBlocked[r][i][k] = DeBlock(A2[r][i], N, nPT, k);
          //counter++;
          //printf("%.1lf%%\r", ((double)counter/tot_cycle)*100);
        }
      }
    }
  }

  //DEBLOCKING MEAN INTENSITIES
  //counter = 0;
  //tot_cycle = nPT*nRep;
  printf("Deblocking means...\n");
  double*** A2M_deBlocked = (double***)calloc(nRep, sizeof(double**));
  checkHypMat(A2M_deBlocked);

  for(int r = 0; r < nRep; r++) {
    A2M_deBlocked[r] = (double**)calloc(nPT, sizeof(double*));
    checkMat(A2M_deBlocked[r]);

    for(int k = 0; k < nPT; k++) {
      A2M_deBlocked[r][k] = DeBlock(A2M[r], N, nPT, k);
    //counter++;
    //printf("%.1lf%%\r", ((double)counter/tot_cycle)*100);
    }

  }


//A2_deBlocked[replicas][iteration][temperatures][mode]
//A2M_deBlocked[replicas][temperatures][mode]
  printf("\n");

  //OverlapIFO_Element(double* Ir1, double* Ir2, double* mr1, double* mr2, int size)
    printf("Computing IFOs....\n");
      //counter = 0;
    //tot_cycle = nPT*nRep*nRep*nFILE;
    double** cVec = (double**)calloc(nPT, sizeof(double*));
  	checkMat(cVec);
  for(int k = 0; k < nPT; k++) {
    char* name_ifo = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(name_ifo);
    sprintf(name_ifo, "ifo_overlaps_temp%d.dat", k);

    FILE* file_ifo = fopen(name_ifo, "w");
    checkFile(file_ifo, name_ifo);
    free(name_ifo);

    fprintf(file_ifo, "#n_iter\tr1\tr2\tC\n");
	 
  	
  	
    if(NITER_MIN_PRINT == 0) {
     cVec[k] = (double*)calloc((int)(nRep*(nRep-1)/2)*((int)(nFILE/2) +1), sizeof(double));
     checkDouble(cVec[k]);
     
	  l = 0;
      for(int i = 0; i < (int)(nFILE/2) +1; i++) {
       
        for(int r1 = 0; r1 < nRep; r1++) {

          for(int r2 = r1+1; r2 < nRep; r2++) {

            cVec[k][l] = OverlapIFO_Element(A2_deBlocked[r1][i][k], A2_deBlocked[r2][i][k], A2M_deBlocked[r1][k], A2M_deBlocked[r2][k], N);
            
            if(r2 > r1) {
              fprintf(file_ifo, "%d\t%d\t%d\t%le\n", i, r1, r2, cVec[k][l]);
            }
            l++;
            //counter++;
            //printf("%.1lf%%\r", ((double)counter/tot_cycle)*100);
          }
        }
      }
    }else{
      cVec[k] = (double*)calloc((int)(nRep*(nRep-1)/2)*nFILE, sizeof(double));
      checkDouble(cVec[k]);
	  l = 0;
      for(int i = 0; i < nFILE; i++) {
        

        for(int r1 = 0; r1 < nRep; r1++) {

          

          for(int r2 = r1+1; r2 < nRep; r2++) {

            cVec[k][l] = OverlapIFO_Element(A2_deBlocked[r1][i][k], A2_deBlocked[r2][i][k], A2M_deBlocked[r1][k], A2M_deBlocked[r2][k], N);
            
            if(r2 < r1) {
              fprintf(file_ifo, "%d\t%d\t%d\t%le\n", i, r1, r2, cVec[k][l]);
            }
            l++;

          }
        }
        fprintf(file_ifo, "\n");
      }
    }
    fclose(file_ifo);
  }
  free(A2_deBlocked);
  free(A2M_deBlocked);
	
	
  return cVec;
}

double** experimental_equilibrium_ifo(complex**** Dynamics, int nRep, int nFILE, int nPT, int N) {
  //FixedDynamics[k][i][r][j]

  double** cmat = (double**)calloc(nPT, sizeof(double*));
  checkMat(cmat);
  double *** A2_mean = (double***)calloc(nPT, sizeof(double**));
  checkHypMat(A2_mean);
  double** A2_rep = (double**)calloc(nPT, sizeof(double*));
  checkMat(A2_rep);


  for(int k = 0; k < nPT; k++) {
    A2_mean[k] = (double**)calloc(nRep, sizeof(double*));
    checkMat(A2_mean[k]);
    for(int r = 0; r < nRep; r++) {
      A2_mean[k][r] = (double*)calloc(N, sizeof(double));
      checkDouble(A2_mean[k][r]);
      for(int j = 0; j < N; j++) {
        for(int i = 0; i < nFILE; i++) {
          A2_mean[k][r][j] += Mod2Value(Dynamics[k][i][r][j])/nFILE;
        }
      }
    }
  }

  for(int k = 0; k < nPT; k++) {
    A2_rep[k] = (double*)calloc(N, sizeof(double));
    checkDouble(A2_rep[k]);
    for(int j = 0; j < N; j++) {
      for(int r = 0; r < nRep; r++) {
        A2_rep[k][j] += A2_mean[k][r][j]/nRep;
      }
    }
  }
  int l = 0;
  //double OverlapIFO_Element(double* Ir1, double* Ir2, double* mr1, double* mr2, int size) {

  for(int k = 0; k < nPT; k++) {
    l = 0;
  char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(filename);
  sprintf(filename, "experimental_ifo_eq_temp%d.dat", k);
  FILE* fw = fopen(filename, "w");
  checkFile(fw, filename);
  free(filename);
    cmat[k] = (double*)calloc((int)(nRep*(nRep-1)/2), sizeof(double));
    checkDouble(cmat[k]);
    for(int r1 = 0; r1 < nRep; r1++) {
      for(int r2 = r1+1; r2< nRep; r2++) {
        cmat[k][l] = OverlapIFO_Element(A2_mean[k][r1], A2_mean[k][r2], A2_rep[k], A2_rep[k], N);
        fprintf(fw, "%d\t%d\t%le\n", r1, r2, cmat[k][l]);
        l++;
      }
    }
    fclose(fw);
  }

  free(A2_mean);
  free(A2_rep);

  return cmat;
}

void dynamics_experimental_ifo(complex**** Dynamics, int nRep, int nFILE, int nPT, int N, int IT_MIN, int dITER, double* temp, int cumulative) {

  int bmax = (int)log2((double)nFILE);
  int it_min_temp = IT_MIN;
  int it_max_temp = 0;
  int nfile_temp = 0;

  //FixedDynamics[k][i][r][j]

  double **** A2_mean = (double****)calloc(bmax, sizeof(double***));
  checkHyp2Mat(A2_mean);
  double*** A2_rep = (double***)calloc(bmax, sizeof(double**));
  checkHypMat(A2_rep);

  for(int b = 0; b < bmax; b++) {
    if(b == 0) {
      it_max_temp = dITER;
    }else{
      it_max_temp = 2*it_min_temp - dITER;
    }

    nfile_temp = (int)((it_max_temp - it_min_temp)/dITER) + 1;

    A2_mean[b] = (double***)calloc(nPT, sizeof(double**));
    checkHypMat(A2_mean[b]);
    int i_min = 0;
    for(int k = 0; k < nPT; k++) {
      A2_mean[b][k] = (double**)calloc(nRep, sizeof(double*));
      checkMat(A2_mean[b][k]);
      for(int r = 0; r < nRep; r++) {
        A2_mean[b][k][r] = (double*)calloc(N, sizeof(double));
        checkDouble(A2_mean[b][k][r]);
        for(int j = 0; j < N; j++) {
          if(cumulative == 0) {
            i_min = (int)(it_min_temp/dITER);
          }else if(cumulative == 1) {
            i_min = 0;
          }
          for(int i = i_min; i < (int)(it_max_temp/dITER)+1; i++) {
            A2_mean[b][k][r][j] += Mod2Value(Dynamics[k][i][r][j])/nfile_temp;
          }
        }
      }
    }

    it_min_temp = it_max_temp + dITER;
  }

  for(int b = 0; b < bmax; b++) {
    A2_rep[b] = (double**)calloc(nPT, sizeof(double*));
    checkMat(A2_rep[b]);
    for(int k = 0; k < nPT; k++) {
      A2_rep[b][k] = (double*)calloc(N, sizeof(double));
      checkDouble(A2_rep[b][k]);
      for(int j = 0; j < N; j++) {
        for(int r = 0; r < nRep; r++) {
          A2_rep[b][k][j] += A2_mean[b][k][r][j]/nRep;
        }
      }
    }
  }

  int l = 0;

  double*** cmat = (double***)calloc(bmax, sizeof(double**));
  checkHypMat(cmat);

  for(int b = 0; b < bmax; b++) {
    cmat[b] = (double**)calloc(nPT, sizeof(double*));
    checkMat(cmat[b]);
    for(int k = 0; k < nPT; k++) {
      cmat[b][k] = (double*)calloc((int)(nRep*(nRep-1)/2), sizeof(double));
      checkDouble(cmat[b][k]);
      l = 0;
      char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
      checkString(filename);
      if(cumulative == 0) {
        sprintf(filename, "ifo_overlaps_block_%d_temp%d.dat", b, k);
      }else if(cumulative == 1) {
        sprintf(filename, "ifo_overlaps_cumulative_block_%d_temp%d.dat", b, k);
      }
      FILE* fw = fopen(filename, "w");
      checkFile(fw, filename);
      free(filename);

      for(int r1 = 0; r1 < nRep; r1++) {
        for(int r2 = r1+1; r2 < nRep; r2++) {
          cmat[b][k][l] = OverlapIFO_Element(A2_mean[b][k][r1], A2_mean[b][k][r2], A2_rep[b][k], A2_rep[b][k], N);
          fprintf(fw, "%d\t%d\t%le\n", r1, r2, cmat[b][k][l]);
          l++;
        }
      }
      fclose(fw);
    }
  }

  for(int b = 0; b < bmax; b++) {
    for(int k = 0; k < nPT; k++) {
      HistogramExpIFO(cmat[b][k], (int)(nRep*(nRep-1)/2), k, temp, b,cumulative);
    }
  }

free(A2_rep);
free(A2_mean);
free(cmat);

}
double** compute_plaquette_overlaps(complex**** FixedDynamics, double** intData, int nRep, int nPT, int Nplaq, int nFILE) {


  int nfile_temp;

  if(NITER_MIN_PRINT == 0) {
    nfile_temp = ((int)(nFILE/2) +1);
  }else{
    nfile_temp = nFILE;
  }

  double e1 = 0.;
  double e2 = 0.;

 
  double** eMat = (double**)calloc(nPT, sizeof(double*));
  checkMat(eMat);
	int l;
  for(int k = 0; k < nPT; k++) { //CICLO TEMPERATURE
    
	l = 0;
	eMat[k] = (double*)calloc((int)(nRep * (nRep-1)/2)*nfile_temp, sizeof(double));

	
    for(int i = 0; i < nfile_temp; i++) {//CICLO ITERAZIONI

      for(int r1 = 0; r1 < nRep; r1++) { //CICLO REPLICA 1

        for(int r2 = 0; r2 < r1; r2++) { //CICLO REPLICA 2
		  
		  


          for(int n = 0; n < Nplaq; n++) {
            e1 = Plaq(FixedDynamics[k][i][r1], intData[n]);
            e2 = Plaq(FixedDynamics[k][i][r2], intData[n]);

            eMat[k][l] += e1*e2/(4*Nplaq);
            
          }
          	if(eMat[k][l] > 1) {
          		printf("OUT OF RANGE!\n");
          		printf("Q[%d][%d][%d][%d] = %lf\n", k, i, r1, r2, eMat[k][l]);
          	}
			l++;
            //QMat[k][i][r1][r2] /= (sqrt(sum1*sum2));
            //printf("norma: %lf\n", sqrt(sum1*sum2));
            //sum1 = 0.;
            //sum2 = 0.;
        }//FINE CICLO REPLICA 2
      }//FINE CICLO REPLICA 1

    }//FINE CICLO ITERAZIONI
  }//FINE CICLO TEMPERATURE

  //PrintMat(QMat[nPT - 1][nfile_temp - 1], nRep, nRep);
  

  

  for(int k = 0; k < nPT; k++) {

    l = 0;
    char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(filename);
    sprintf(filename, "plaq_overlaps_temp%d.dat", k);
    FILE* fw = fopen(filename, "w");
    checkFile(fw, filename);
    fprintf(fw, "#T\tITER\tr1\tr2\tQ\n");
    for(int i = 0; i < nfile_temp; i++) {
      for(int r1 = 0; r1 < nRep; r1++) {
        for(int r2 = 0; r2 < r1; r2++) {
          
          fprintf(fw, "%d\t%d\t%d\t%d\t%le\n",k,i,r1,r2, eMat[k][l]);
          l++;
        }
      }
      fprintf(fw, "\n");
    }
    fclose(fw);
    free(filename);
  }

  return eMat;
}





//DISORDER AVERAGE TOOLS
void plot_disorder_average_pq(char* filename, int k, int block, double mean, double devstd, double kurt,double errK, double skew, double errS) {

  char* tempfile_name = (char*)calloc(STRING_SIZE, sizeof(double));
  checkString(tempfile_name);
  sprintf(tempfile_name, "plot_script.p");

  FILE* script = fopen(tempfile_name, "w");
  checkFile(script, tempfile_name);
  free(tempfile_name);
  fprintf(script, "set terminal png\n");
  if(NITER_MIN_PRINT == 0) {
    fprintf(script, "set output 'dis_ave_pq_temp%d_block%d.png'\n", k, block);
  }else{
    fprintf(script, "set output 'dis_ave_pq_temp%d.png'\n", k);
  }
  fprintf(script, "set title 'Overlaps Distribution'\n");
  fprintf(script, "set log y 10\n");
  fprintf(script, "set yrange[0.001:10]\n");
  fprintf(script, "set xlabel 'q'\n");
  fprintf(script, "set ylabel 'P(q)'\n");
  char* title = (char*)calloc(10*STRING_SIZE, sizeof(char));
  checkString(title);
  sprintf(title,"mean: %g+/-%g, skewness: %g+/-%g kurtosis %g+/-%g", mean,devstd, skew, errS, kurt, errK);
  fprintf(script, "set title '%s' font ',8'\n", title);
  free(title);
  if(block == -1) {
    fprintf(script, "plot '%s' u 1:2:3 w ye lt 7 notitle, '%s' w lp lt 7 title 'T = %d'\n", filename, filename, k);
  }else{
    fprintf(script, "plot '%s' u 1:2:3 w ye lt 7 notitle, '%s' w lp lt 7 title 'mcs = %d, T %d'\n", filename, filename, (int)pow(2,block+1), k);
  }

  fclose(script);
  system("gnuplot plot_script.p");
  system("rm plot_script.p");

}
double** acquiring_freq(int nsamples, char* location) {

  double** freq = (double**)calloc(nsamples, sizeof(double*));
  checkMat(freq);

  for(int n = 0; n < nsamples; n++) { //ciclo sui samples

    freq[n] = (double*)calloc(Size, sizeof(double));
    checkDouble(freq[n]);

    location = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(location);

    sprintf(location, "sample%d/frequencies.dat", n+1);
    FILE* freq_file = fopen(location, "r");
    checkFile(freq_file, location);
    int index_temp;
    double gain_temp;
    for(int j = 0; j < Size; j++) { //ciclo sui modi

      fscanf(freq_file, "%d\t%le\t%le\n", &index_temp, &freq[n][j], &gain_temp);
      //printf("%le\n", freq[n][j]);
    } //fine ciclo sui modi

    fclose(freq_file);
  }// fine ciclo sui samples


  return freq;

}


void plot_dis_ave_spectrum(char* location, int k, int r) {

  char* script_file_name = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(script_file_name);

  sprintf(script_file_name, "script_spectrum.p");
  FILE* script_file = fopen(script_file_name, "w");
  checkFile(script_file, script_file_name);


  fprintf(script_file, "set term png\n");
  fprintf(script_file, "set encoding utf8\n");
  fprintf(script_file, "set xlabel 'ω'\n");
  fprintf(script_file, "set ylabel 'I(ω)'\n");
  fprintf(script_file, "set yrange[0:0.5]\n");

  fprintf(script_file, "set title 'Random Laser Intensity Spectrum'\n");

  fprintf(script_file, "set output 'dis_ave_spectrum_nrep%d_temp%d.png'\n", r, k);

  fprintf(script_file, "plot '%s' w lp lt 7 title 'T %d'\n", location, k);

  fclose(script_file);

  system("gnuplot script_spectrum.p");

  system("rm script_spectrum.p");


}
void plot_dis_ave_spectrum_rebinned(char* location, int k, int r) {

  char* script_file_name = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(script_file_name);

  sprintf(script_file_name, "script_spectrum.p");
  FILE* script_file = fopen(script_file_name, "w");
  checkFile(script_file, script_file_name);


  fprintf(script_file, "set term png\n");
  fprintf(script_file, "set encoding utf8\n");
  fprintf(script_file, "set xlabel 'ω'\n");
  fprintf(script_file, "set ylabel 'I(ω)'\n");
  fprintf(script_file, "set yrange[0:0.3]\n");

  fprintf(script_file, "set title 'Random Laser Intensity Spectrum'\n");

  fprintf(script_file, "set output 'dis_ave_spectrum_rebinned_nrep%d_temp%d.png'\n", r, k);

  fprintf(script_file, "plot '%s' u 1:2:3 w ye lt 7 notitle, '%s' w lp lt 7 title 'T %d'\n", location,location, k);

  fclose(script_file);

  system("gnuplot script_spectrum.p");

  system("rm script_spectrum.p");


}
void disorder_average_spectrum_histogram(double* vec, double* intensity, int size, int r, int k) {

  double min = 0.;
  double max = 1.;
  int threshold_hist = 2;
  double bin_size = (max - min)/(HISTBARS/threshold_hist);
  double* omega_mean = (double*)calloc(HISTBARS/threshold_hist, sizeof(double));
  checkDouble(omega_mean);
  double* intensity_mean = (double*)calloc(HISTBARS/threshold_hist, sizeof(double));
  checkDouble(intensity_mean);
  double* sigma_intensity = (double*)calloc(HISTBARS/threshold_hist, sizeof(double));
  checkDouble(sigma_intensity);

  int* counter = (int*)calloc(HISTBARS/threshold_hist, sizeof(int));
  checkInt(counter);
  double* temp_int;

  char* hist_file_name = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(hist_file_name);

  sprintf(hist_file_name, "dis_ave_spectrum_rebinned_nrep%d_temp%d.dat", r, k);
  FILE* hist_file = fopen(hist_file_name, "w");
  checkFile(hist_file, hist_file_name);


  fprintf(hist_file, "#omega\tintensity\tsigma\n");


  for(int i = 0; i < HISTBARS/threshold_hist; i++) {
    for(int j = 0; j < size; j++) {
      if((vec[j] >= min + i*bin_size)&&(vec[j] < min + (i+1)*bin_size)) {
        counter[i]++;
      }
    }
    omega_mean[i] = min + bin_size/2 + i*bin_size;

    //printf("counter[%d]: %d\n", i, counter[i]);
  }
  int c;
  for(int i = 0; i < HISTBARS/threshold_hist; i++) {
    if(counter[i] > 0) {
      temp_int = (double*)calloc(counter[i], sizeof(double));
      checkDouble(temp_int);
      c = 0;
      for(int j = 0; j < size; j++) {
        if((vec[j] >= min + i*bin_size)&&(vec[j] < min + (i+1)*bin_size)) {
          temp_int[c] = intensity[j];
          c++;
        }
      }
      intensity_mean[i] = MeanVec(temp_int, c);
      sigma_intensity[i] = StandardDeviation(temp_int, intensity_mean[i], c);
      free(temp_int);

      fprintf(hist_file, "%le\t%le\t%le\n", omega_mean[i], intensity_mean[i], sigma_intensity[i]);
    }
  }

  fclose(hist_file);
  plot_dis_ave_spectrum_rebinned(hist_file_name, k, r);

}
void all_spectrum_random(int nsamples) {


  char* location;
  int counter = 0;
  double*** sp_random = (double***)calloc(NREPLICAS, sizeof(double**));
  checkHypMat(sp_random);
  double*** err_random = (double***)calloc(NREPLICAS, sizeof(double**));
  checkHypMat(err_random);
  double* freq = (double*)calloc(Size*nsamples, sizeof(double));
  checkDouble(freq);


  for(int r = 0; r < NREPLICAS; r++) { //ciclo sulle repliche

    sp_random[r] = (double**)calloc(NPT, sizeof(double*));
    checkMat(sp_random[r]);
    err_random[r] = (double**)calloc(NPT, sizeof(double*));
    checkMat(err_random[r]);
    for(int k = 0; k < NPT; k++) { //ciclo sulle temprature

      counter = 0;
      sp_random[r][k] = (double*)calloc(nsamples*Size, sizeof(double));
      checkDouble(sp_random[r][k]);
      err_random[r][k] = (double*)calloc(nsamples*Size, sizeof(double));
      checkDouble(err_random[r][k]);
      for(int n = 0; n < nsamples; n++) { //ciclo sui campioni
        location = (char*)calloc(2*STRING_SIZE, sizeof(char));
        checkString(location);
		//printf("CHECK\n");
        sprintf(location, "sample%d/mean_spectrum/intensity_nrep%d_temp%d.dat", n+1, r, k);
        FILE* spectra_file = fopen(location, "r");
        for(int j = 0; j < Size; j++) {
          double temp_f = 0.;
          double temp_i = 0.;
          double temp_e = 0.;

          fscanf(spectra_file, "%le\t%le\t%le\n", &temp_f, &temp_i, &temp_e);
          sp_random[r][k][counter] = temp_i;
          err_random[r][k][counter] = temp_e;
          freq[counter] = temp_f;

          counter++;
          //printf("%le\n", spectra[r][k][n][j]);
        }

        fclose(spectra_file);
		


      }//fine ciclo sui campioni
    }//fine ciclo sulle temperature
  }//fine ciclo sulle repliche

  double a = 0.;
  double b = 0.;
  double c = 0.;

  for(int n = 0; n < nsamples*Size; n++) {
    for(int m = n+1; m < nsamples*Size; m++) {
      if (freq[n] > freq[m]){
            a = freq[n];
            freq[n] = freq[m];
            freq[m] = a;
            for(int k = 0; k < NPT; k++) {

              for(int r = 0; r < NREPLICAS; r++) {
                b = sp_random[r][k][n];
                sp_random[r][k][n] = sp_random[r][k][m];
                sp_random[r][k][m] = b;
                c = err_random[r][k][n];
                err_random[r][k][n] = err_random[r][k][m];
                err_random[r][k][m] = c;
              }
            }
         }
    }
  }
  //REBINNING DEGLI SPETTRI
  printf("\nRebinning spectrum...\n");
  for(int r = 0; r < NREPLICAS; r++) {
    for(int k = 0; k < NPT; k++) {
      disorder_average_spectrum_histogram(freq, sp_random[r][k], Size*nsamples, r, k);
    }
  }
  FILE* mean_spectra_file;

  for(int r = 0; r < NREPLICAS; r++) {
    for(int k = 0; k < NPT; k++) {
      location = (char*)calloc(2*STRING_SIZE, sizeof(char));
      checkString(location);
      sprintf(location, "dis_ave_spectrum_nrep%d_temp%d.dat", r, k);
      mean_spectra_file = fopen(location, "w");
      checkFile(mean_spectra_file, location);
      for(int n = 0; n < nsamples*Size; n++) {

        fprintf(mean_spectra_file, "%le\t%le\t%le\n", freq[n], sp_random[r][k][n], err_random[r][k][n]);
        //printf("freq[%d]: %le\n", n, freq[n]);
      }


      fclose(mean_spectra_file);
      plot_dis_ave_spectrum(location, k, r);

      char* command = (char*)calloc(STRING_SIZE, sizeof(char));
      checkString(command);
      sprintf(command, "mv %s DIS_AVE/spectrum", location);
      system(command);
      free(command);
    }
  }

  free(freq);
  free(sp_random);
  free(err_random);
  free(location);

}
void disorder_average_spectrum_comblike(int nsamples, char* location, double** freq) {

  double**** spectra = (double****)calloc(NREPLICAS, sizeof(double***));
  checkHyp2Mat(spectra);

  for(int r = 0; r < NREPLICAS; r++) { //ciclo sulle repliche
    spectra[r] = (double***)calloc(NPT, sizeof(double**));
    checkHypMat(spectra[r]);
    for(int k = 0; k < NPT; k++) { //ciclo sulle temprature

      spectra[r][k] = (double**)calloc(nsamples, sizeof(double*));
      checkMat(spectra[r][k]);
      for(int n = 0; n < nsamples; n++) { //ciclo sui campioni
        spectra[r][k][n] = (double*)calloc(Size, sizeof(double));
        location = (char*)calloc(2*STRING_SIZE, sizeof(char));
        checkString(location);
        sprintf(location, "sample%d/mean_spectrum/intensity_nrep%d_temp%d.dat", n+1, r, k);
        FILE* spectra_file = fopen(location, "r");
        for(int j = 0; j < Size; j++) {
          double temp_f = 0.;
          double temp_i = 0.;
          double temp_e = 0.;

          fscanf(spectra_file, "%le\t%le\t%le\n", &temp_f, &temp_i, &temp_e);
          spectra[r][k][n][j] = temp_i;
          //printf("%le\n", spectra[r][k][n][j]);
        }

        fclose(spectra_file);

      }//fine ciclo sui campioni
    }//fine ciclo sulle temperature
  }//fine ciclo sulle repliche

  for(int r = 0; r < NREPLICAS; r++) {
    for(int k = 0; k < NPT; k++) {
      spectra[r][k] = TransposeMat(spectra[r][k], nsamples, Size);
    }
  }
  FILE* mean_spectra_file;
  double*** mean_spectra = (double***)calloc(NREPLICAS, sizeof(double**));
  checkHypMat(mean_spectra);
  double*** err_spectra = (double***)calloc(NREPLICAS, sizeof(double**));
  checkHypMat(err_spectra);
  for(int r = 0; r < NREPLICAS; r++) { //ciclo sulle repliche
    mean_spectra[r] = (double**)calloc(NPT, sizeof(double*));
    checkMat(mean_spectra[r]);

    err_spectra[r] = (double**)calloc(NPT, sizeof(double*));
    checkMat(err_spectra[r]);

    for(int k = 0; k < NPT; k++) { //ciclo sulle temperature
      mean_spectra[r][k] = (double*)calloc(Size, sizeof(double));
      checkDouble(mean_spectra[r][k]);

      err_spectra[r][k] = (double*)calloc(Size, sizeof(double));
      checkDouble(err_spectra[r][k]);

      sprintf(location, "dis_ave_spectrum_nrep%d_temp%d.dat", r, k);
      mean_spectra_file = fopen(location, "w");
      checkFile(mean_spectra_file, location);
      for(int j = 0; j < Size; j++) { //ciclo sui modi
        for(int n = 0; n < nsamples; n++) { //ciclo sui campioni
          mean_spectra[r][k][j] += spectra[r][k][j][n];
        } //fine ciclo sui campioni

        mean_spectra[r][k][j] /= nsamples;
        err_spectra[r][k][j] = StandardDeviation(spectra[r][k][j], mean_spectra[r][k][j], nsamples);

        fprintf(mean_spectra_file, "%le\t%le\t%le\n", freq[0][j], mean_spectra[r][k][j], err_spectra[r][k][j]);

      }//fine ciclo sui modi
      fclose(mean_spectra_file);
      plot_dis_ave_spectrum(location, k, r);

      char* command = (char*)calloc(STRING_SIZE, sizeof(char));
      checkString(command);
      sprintf(command, "mv %s DIS_AVE/spectrum", location);
      system(command);
      free(command);
    }//fine ciclo sulle temperature

  }//fine ciclo sulle repliche
  location = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(location);


  free(mean_spectra);
  free(err_spectra);
  free(spectra);
}

void plot_disorder_average_sh(char* filename, int k) {

  char* tempfile_name = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(tempfile_name);
  sprintf(tempfile_name, "plot_script.p");

  FILE* script = fopen(tempfile_name, "w");
  checkFile(script,tempfile_name);
  free(tempfile_name);
  fprintf(script, "set terminal png\n");

  fprintf(script, "set output 'dis_ave_specific_heat_nrep%d.png'\n", k);


  //fprintf(script, "set log y 10\n");
  //fprintf(script, "set yrange[0.001:100]\n");
  fprintf(script, "set xlabel 'T'\n");
  fprintf(script, "set ylabel 'C_V(T)'\n");

  fprintf(script, "plot '%s' u 1:2:3 w ye lt 7 notitle, '%s' w lp lt 7 notitle\n", filename, filename);

  fclose(script);
  system("gnuplot plot_script.p");
  system("rm plot_script.p");
}

void disorder_average_specific_heat(int nsamples, char* location) {

  double*** sh = (double***)calloc(NREPLICAS, sizeof(double**));
  checkHypMat(sh);

  double* temp = (double*)calloc(NPT-1, sizeof(double));
  checkDouble(temp);

  for(int r = 0; r < NREPLICAS; r++) { //ciclo sulle repliche
    sh[r] = (double**)calloc(nsamples, sizeof(double*));
    checkMat(sh[r]);

    for(int n = 0; n < nsamples; n++) {//ciclo sui campioni
      sh[r][n] = (double*)calloc(NPT-1, sizeof(double));
      checkDouble(sh[r][n]);

      location = (char*)calloc(STRING_SIZE, sizeof(char));
      checkString(location);
      sprintf(location, "sample%d/specific_heat/specific_heat_nrep%d_size%d.dat",n+1, r, Size);

      FILE* spec_heat_file = fopen(location, "r");
      checkFile(spec_heat_file, location);

      for(int k = 0; k < NPT-1; k++) { //ciclo sulle temperature
        double t_err;
        fscanf(spec_heat_file, "%le\t%le\t%le\n", &temp[k], &sh[r][n][k], &t_err);
      } //fine ciclo sulle temperature

      fclose(spec_heat_file);

    }//fine ciclo sui campioni

  }//fine ciclo sulle repliche

  double** mean_sh = (double**)calloc(NREPLICAS, sizeof(double*));
  checkMat(mean_sh);
  double** err_sh = (double**)calloc(NREPLICAS, sizeof(double*));
  checkMat(err_sh);


  for(int r = 0; r < NREPLICAS; r++) {
    sh[r] = TransposeMat(sh[r], nsamples, NPT-1);

    mean_sh[r] = (double*)calloc(NPT-1, sizeof(double));
    checkDouble(mean_sh[r]);

    err_sh[r] = (double*)calloc(NPT-1, sizeof(double));
    checkDouble(err_sh[r]);

    location = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(location);

    sprintf(location, "dis_ave_specific_heat_nrep%d_size%d.dat", r, Size);
    FILE* spec_heat_file = fopen(location, "w");
    checkFile(spec_heat_file, location);

    for(int k = 0; k < NPT-1; k++) {


      for(int n = 0; n < nsamples; n++) {
        mean_sh[r][k] += sh[r][k][n];
      }
      mean_sh[r][k] /= nsamples;
      err_sh[r][k] = StandardDeviation(sh[r][k], mean_sh[r][k], nsamples);
      fprintf(spec_heat_file, "%le\t%le\t%le\n", temp[k], mean_sh[r][k], err_sh[r][k]);
    }
    fclose(spec_heat_file);

    plot_disorder_average_sh(location,r);
    char* command = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(command);
    sprintf(command, "mv %s DIS_AVE/specific_heat", location);
    system(command);
    free(command);
  }

  free(mean_sh);

  free(err_sh);
  free(sh);
  free(temp);



}

void overlaps_hist_disorder_average(int nsamples, int block) {

  char* location;
  double*** pq = (double***)calloc(NPT, sizeof(double**));
  checkHypMat(pq);
  double* q = (double*)calloc(HISTBARS, sizeof(double));
  checkDouble(q);
  double max = 1.01;
  double min = -1.01;

  double bin_size = (max-min)/HISTBARS;


  for(int k = 0; k < NPT; k++) {
    pq[k] = (double**)calloc(nsamples, sizeof(double*));
    checkMat(pq[k]);

    for(int n = 0; n < nsamples; n++){
      pq[k][n] = (double*)calloc(HISTBARS, sizeof(double));
      checkDouble(pq[k][n]);

      location = (char*)calloc(STRING_SIZE, sizeof(char));
      checkString(location);
      if(block == -1) {
        sprintf(location, "sample%d/PqRE/histograms/p_q_re_temp%d.dat", n+1, k);
      }else{
        sprintf(location, "sample%d/PqRE/histograms/p_q_re_temp%d_block%d.dat", n+1, k, block);
      }
      FILE* pq_file = fopen(location, "r");
      checkFile(pq_file, location);
      fscanf(pq_file, "%*[^\n]\n");
      for(int j = 0; j < HISTBARS; j++) {
        fscanf(pq_file, "%le\t%le\n", &q[j], &pq[k][n][j]);
      }
      fclose(pq_file);

    }
  }



  double** mean_pq = (double**)calloc(NPT, sizeof(double*));
  checkMat(mean_pq);
  double** err_pq = (double**)calloc(NPT, sizeof(double*));
  checkMat(err_pq);

  for(int k = 0; k < NPT; k++) {

    pq[k] = TransposeMat(pq[k], nsamples, HISTBARS);
    mean_pq[k] = (double*)calloc(HISTBARS, sizeof(double));
    checkDouble(mean_pq[k]);
    err_pq[k] = (double*)calloc(HISTBARS, sizeof(double));
    checkDouble(err_pq[k]);
    location = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(location);
    if(block == -1) {
      sprintf(location, "dis_ave_p_q_temp%d.dat", k);
    }else{
      sprintf(location, "dis_ave_p_q_temp%d_block%d.dat", k, block);
    }
    FILE* pq_file = fopen(location, "w");
    checkFile(pq_file, location);

    for(int j = 0; j < HISTBARS; j++) {

      for(int n = 0; n < nsamples; n++) {
        mean_pq[k][j] += pq[k][j][n];
      }

      mean_pq[k][j] /= nsamples;
      err_pq[k][j] = StandardDeviation(pq[k][j], mean_pq[k][j], nsamples);

    }
    double mean = mean_of_distribution(q, mean_pq[k],bin_size, HISTBARS);
    double devStd = StandardDeviation(q, mean, HISTBARS);
    double skew = Skewness(q, mean_pq[k], bin_size, HISTBARS);
    double errS = errSkewness(q, mean_pq[k], bin_size, HISTBARS, skew);
    double kurt = Kurtosis(q, mean_pq[k], bin_size, HISTBARS);
    double errK = errKurtosis(q, mean_pq[k], bin_size, HISTBARS, kurt);

    fprintf(pq_file, "#mean: %g +/- %g\tskewness: %g +/- %g kurtosis: %g +/- %g\n", mean,devStd, skew, errS, kurt, errK);
    for(int j = 0; j < HISTBARS; j++) {
      fprintf(pq_file, "%le\t%le\t%le\n", q[j], mean_pq[k][j], err_pq[k][j]);
    }
    fclose(pq_file);
    plot_disorder_average_pq(location, k, block, mean, devStd, kurt,errK, skew, errS);

    char* command = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(command);
    sprintf(command, "mv %s DIS_AVE/PqRE/", location);
    system(command);
    free(command);
  }

  free(mean_pq);
  free(pq);
  free(q);
  free(err_pq);
  free(location);
}

void plot_disorder_average_ifo(char* filename, int k) {

  char* tempfile_name = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(tempfile_name);
  sprintf(tempfile_name, "plot_script.p");

  FILE* script = fopen(tempfile_name, "w");
  checkFile(script,tempfile_name);
  free(tempfile_name);
  fprintf(script, "set terminal png\n");

  fprintf(script, "set output 'dis_ave_ifo_temp%d.png'\n", k);


  fprintf(script, "set log y 10\n");
  fprintf(script, "set yrange[0.001:10]\n");
  fprintf(script, "set xlabel 'C'\n");
  fprintf(script, "set ylabel 'P(C)'\n");
  fprintf(script, "set title 'Intensity Fluctuations Overlaps Distribution'\n");

  fprintf(script, "plot '%s' u 1:2:3 w ye lt 7 notitle, '%s' w lp lt 7 title 'T %d'\n", filename, filename, k);

  fclose(script);
  system("gnuplot plot_script.p");
  system("rm plot_script.p");
}

void plot_exp_disorder_average_ifo(char* filename, int k, int block, int cumulative) {

  char* tempfile_name = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(tempfile_name);
  sprintf(tempfile_name, "plot_script.p");

  FILE* script = fopen(tempfile_name, "w");
  checkFile(script,tempfile_name);
  free(tempfile_name);
  fprintf(script, "set terminal png\n");

  if(cumulative == 0) {
    fprintf(script, "set output 'dis_ave_exp_IFO_block%d_temp%d.png'\n", block, k);
  }else if (cumulative == 1){
    fprintf(script, "set output 'dis_ave_exp_IFO_cumulative_block%d_temp%d.png'\n", block, k);

  }else if(cumulative == 2) {
      fprintf(script, "set output 'dis_ave_exp_IFO_equilibrium_temp%d.png'\n", k);
  }

  fprintf(script, "set log y 10\n");
  fprintf(script, "set yrange[0.001:10]\n");
  fprintf(script, "set xlabel 'C'\n");
  fprintf(script, "set ylabel 'P(C)'\n");
  fprintf(script, "set title 'Intensity Fluctuations Overlaps Distribution'\n");

  fprintf(script, "plot '%s' u 1:2:3 w ye lt 7 notitle, '%s' w lp lt 7 title 'T %d'\n", filename, filename, k);

  fclose(script);
  system("gnuplot plot_script.p");
  system("rm plot_script.p");
}
void ifo_hist_disorder_average(int nsamples) {

  char* location;
  double*** pq = (double***)calloc(NPT, sizeof(double**));
  checkHypMat(pq);
  double* q = (double*)calloc(HISTBARS, sizeof(double));
  checkDouble(q);


  for(int k = 0; k < NPT; k++) {
    pq[k] = (double**)calloc(nsamples, sizeof(double*));
    checkMat(pq[k]);

    for(int n = 0; n < nsamples; n++){
      pq[k][n] = (double*)calloc(HISTBARS, sizeof(double));
      checkDouble(pq[k][n]);

      location = (char*)calloc(STRING_SIZE, sizeof(char));
      checkString(location);
      sprintf(location, "sample%d/IFOs/histograms/IFO_temp%d.dat", n+1, k);
      FILE* pq_file = fopen(location, "r");
      checkFile(pq_file, location);
      //fscanf(pq_file, "%*[^\n]\n");
      for(int j = 0; j < HISTBARS; j++) {
        fscanf(pq_file, "%le\t%le\n", &q[j], &pq[k][n][j]);
      }
      fclose(pq_file);

    }
  }



  double** mean_pq = (double**)calloc(NPT, sizeof(double*));
  checkMat(mean_pq);
  double** err_pq = (double**)calloc(NPT, sizeof(double*));
  checkMat(err_pq);

  for(int k = 0; k < NPT; k++) {

    pq[k] = TransposeMat(pq[k], nsamples, HISTBARS);
    mean_pq[k] = (double*)calloc(HISTBARS, sizeof(double));
    checkDouble(mean_pq[k]);
    err_pq[k] = (double*)calloc(HISTBARS, sizeof(double));
    checkDouble(err_pq[k]);
    location = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(location);
    sprintf(location, "dis_ave_IFO_temp%d.dat", k);
    FILE* pq_file = fopen(location, "w");
    checkFile(pq_file, location);

    for(int j = 0; j < HISTBARS; j++) {

      for(int n = 0; n < nsamples; n++) {
        mean_pq[k][j] += pq[k][j][n];
      }

      mean_pq[k][j] /= nsamples;
      err_pq[k][j] = StandardDeviation(pq[k][j], mean_pq[k][j], nsamples);
      fprintf(pq_file, "%le\t%le\t%le\n", q[j], mean_pq[k][j], err_pq[k][j]);

    }

    fclose(pq_file);
    plot_disorder_average_ifo(location, k);
    char* command = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(command);
    sprintf(command, "mv %s DIS_AVE/IFOs/", location);
    system(command);
    free(command);
  }

  free(mean_pq);
  free(pq);
  free(q);
  free(err_pq);
  free(location);
}


void dynamics_ifo_hist_disorder_average(int nsamples, int block, int cumulative) {

  char* location;
  double*** pq = (double***)calloc(NPT, sizeof(double**));
  checkHypMat(pq);
  double* q = (double*)calloc(HISTBARS, sizeof(double));
  checkDouble(q);


  for(int k = 0; k < NPT; k++) {
    pq[k] = (double**)calloc(nsamples, sizeof(double*));
    checkMat(pq[k]);

    for(int n = 0; n < nsamples; n++){
      pq[k][n] = (double*)calloc(HISTBARS, sizeof(double));
      checkDouble(pq[k][n]);

      location = (char*)calloc(STRING_SIZE, sizeof(char));
      checkString(location);
      if(cumulative == 0) {
        sprintf(location, "sample%d/IFOs/histograms/experimental/blocks_indep/IFO_exp_block%d_temp%d.dat", n+1, block, k);
      }else if(cumulative == 1){
        sprintf(location, "sample%d/IFOs/histograms/experimental/cumulative/IFO_cumulative_exp_block%d_temp%d.dat", n+1, block, k);

      }else if(cumulative == 2) {
              sprintf(location, "sample%d/IFOs/histograms/experimental/equilibrium/IFO_exp_temp%d.dat", n+1, k);
      }
      FILE* pq_file = fopen(location, "r");
      checkFile(pq_file, location);
      //fscanf(pq_file, "%*[^\n]\n");
      for(int j = 0; j < HISTBARS; j++) {
        fscanf(pq_file, "%le\t%le\n", &q[j], &pq[k][n][j]);
      }
      fclose(pq_file);

    }
  }



  double** mean_pq = (double**)calloc(NPT, sizeof(double*));
  checkMat(mean_pq);
  double** err_pq = (double**)calloc(NPT, sizeof(double*));
  checkMat(err_pq);

  for(int k = 0; k < NPT; k++) {

    pq[k] = TransposeMat(pq[k], nsamples, HISTBARS);
    mean_pq[k] = (double*)calloc(HISTBARS, sizeof(double));
    checkDouble(mean_pq[k]);
    err_pq[k] = (double*)calloc(HISTBARS, sizeof(double));
    checkDouble(err_pq[k]);
    location = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(location);
    if(cumulative == 0) {
      sprintf(location, "dis_ave_exp_IFO_block%d_temp%d.dat", block, k);
    }else if(cumulative == 1) {
      sprintf(location, "dis_ave_exp_IFO_cumulative_block%d_temp%d.dat", block, k);

    }else if(cumulative == 2) {
          sprintf(location, "dis_ave_exp_IFO_eq_temp%d.dat", k);
    }
    FILE* pq_file = fopen(location, "w");
    checkFile(pq_file, location);

    for(int j = 0; j < HISTBARS; j++) {

      for(int n = 0; n < nsamples; n++) {
        mean_pq[k][j] += pq[k][j][n];
      }

      mean_pq[k][j] /= nsamples;
      err_pq[k][j] = StandardDeviation(pq[k][j], mean_pq[k][j], nsamples);
      fprintf(pq_file, "%le\t%le\t%le\n", q[j], mean_pq[k][j], err_pq[k][j]);

    }

    fclose(pq_file);
    plot_exp_disorder_average_ifo(location, k, block, cumulative);
    char* command = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(command);
    if(cumulative == 0) {
      sprintf(command, "mv %s DIS_AVE/IFOs/experimental/blocks_indep", location);
    }else if(cumulative == 1){
      sprintf(command, "mv %s DIS_AVE/IFOs/experimental/cumulative", location);

    }else if(cumulative == 2) {
      sprintf(command, "mv %s DIS_AVE/IFOs/experimental/equilibrium", location);
    }
    system(command);
    free(command);
  }

  free(mean_pq);
  free(pq);
  free(q);
  free(err_pq);
  free(location);
}




void plot_disorder_average_plaqs(char* filename, int k) {

  char* tempfile_name = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(tempfile_name);
  sprintf(tempfile_name, "plot_script.p");
	//int max = 3;
  FILE* script = fopen(tempfile_name, "w");
  checkFile(script,tempfile_name);

  fprintf(script, "set terminal png\n");

  fprintf(script, "set output 'dis_ave_plaqs_temp%d.png'\n", k);


  fprintf(script, "set log y 10\n");
  fprintf(script, "set yrange[0.001:100]\n");
  fprintf(script, "set xlabel 'Q'\n");
  fprintf(script, "set ylabel 'P(Q)'\n");
  fprintf(script, "set title 'Plaquette Overlaps Distribution'\n");

  fprintf(script, "plot '%s' u 1:2:3 w ye lt 7 notitle, '%s' w lp lt 7 title 'T %d'\n", filename, filename, k);

  fclose(script);
  system("gnuplot plot_script.p");
  system("rm plot_script.p");
}

void plaq_hist_disorder_average(int nsamples) {

  char* location;
  double*** pq = (double***)calloc(NPT, sizeof(double**));
  checkHypMat(pq);
   int max = 3;

  double* q = (double*)calloc(max*HISTBARS, sizeof(double));
  checkDouble(q);
 
  for(int k = 0; k < NPT; k++) {
    pq[k] = (double**)calloc(nsamples, sizeof(double*));
    checkMat(pq[k]);

    for(int n = 0; n < nsamples; n++){
      pq[k][n] = (double*)calloc(max*HISTBARS, sizeof(double));
      checkDouble(pq[k][n]);

      location = (char*)calloc(STRING_SIZE, sizeof(char));
      checkString(location);
      sprintf(location, "sample%d/PQ_plaqs/histograms/EPS_temp%d.dat", n+1, k);
      FILE* pq_file = fopen(location, "r");
      checkFile(pq_file, location);
      //fscanf(pq_file, "%*[^\n]\n");
      for(int j = 0; j < max*HISTBARS; j++) {
        fscanf(pq_file, "%le\t%le\n", &q[j], &pq[k][n][j]);
      }
      fclose(pq_file);

    }
  }



  double** mean_pq = (double**)calloc(NPT, sizeof(double*));
  checkMat(mean_pq);
  double** err_pq = (double**)calloc(NPT, sizeof(double*));
  checkMat(err_pq);

  for(int k = 0; k < NPT; k++) {

    pq[k] = TransposeMat(pq[k], nsamples, max*HISTBARS);
    mean_pq[k] = (double*)calloc(max*HISTBARS, sizeof(double));
    checkDouble(mean_pq[k]);
    err_pq[k] = (double*)calloc(max*HISTBARS, sizeof(double));
    checkDouble(err_pq[k]);
    location = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(location);
    sprintf(location, "dis_ave_EPS_temp%d.dat", k);
    FILE* pq_file = fopen(location, "w");
    checkFile(pq_file, location);

    for(int j = 0; j < max*HISTBARS; j++) {

      for(int n = 0; n < nsamples; n++) {
        mean_pq[k][j] += pq[k][j][n];
      }

      mean_pq[k][j] /= nsamples;
      err_pq[k][j] = StandardDeviation(pq[k][j], mean_pq[k][j], nsamples);
      fprintf(pq_file, "%le\t%le\t%le\n", q[j], mean_pq[k][j], err_pq[k][j]);

    }

    fclose(pq_file);
    plot_disorder_average_plaqs(location, k);
    char* command = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(command);
    sprintf(command, "mv %s DIS_AVE/PQ_plaqs/", location);
    system(command);
    free(command);
  }

  free(mean_pq);
  free(pq);
  free(q);
  free(err_pq);
  free(location);
}

//NEW MACRO FUNCTIONS
/*
double*** read_parallel_tempering(int iter) {
	
	double*** erg = (double***)calloc(NREPLICAS, sizeof(double***));
	checkHypMat(erg);
	for(int r = 0; r < NREPLICAS; r++) {
		erg[r] = (double**)calloc(NPT-1, sizeof(double*));
		checkMat(erg[r]);
		for(int k = 0; k < NPT-1; k++) {
			erg[r][k] = (double*)calloc(iter, sizeof(double));
			checkDouble(erg[r][k]);
		}
	}
	
	double* temperature = (double*)calloc(NPT-1, sizeof(double));
	checkDouble(temperature);
	
	//FACCIAMO UN CICLO SUI FILE, CHE CORRISPONDONO ALLE VARIE REPLICHE
	for(int r = 0; r < NREPLICAS; r++) {
		char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
		checkString(filename);
		sprintf(filename, "parallel_tempering%d.dat", r);
		
		FILE* file_parallel = fopen(filename, "r");
		checkFile(file_parallel, filename);
		
		//CICLO SULLE RIGHE DEL FILE, CHE SONO I PASSI MONTE CARLO
		
		for(int t = 0; t < iter; t++) {
			
			int mcs;
			fscanf(file_parallel, " %d", &mcs);
			
			//CICLO SULLE COLONNE DEL FILE, CHE SONO LE TEMPERATURE
			
			for(int k = NPT-2; k >= 0; k--) {
					double temp = 0.;
					double t_erg = 0.;
					
					double acc = 0.;
					if(k == 0) {
						fscanf(file_parallel, " \t %lf %le %le\n", &temp, &acc, &t_erg);
					}else{
						fscanf(file_parallel, " \t %lf %le %le", &temp, &acc, &t_erg);
					}
					//printf("erg[%d][%d][%d] = %le\n", r, k, t, t_erg);
					erg[r][k][t] = t_erg * Size;     //Perchè i valori nei files sono delle energie per volume
					if(r == 0 && t == 0) {
						temperature[k] = temp;
					}
			}//FINE FOR RIGHE
				
		}//FINE FOR COLONNE
		
		fclose(file_parallel);
		
	}//FINE CICLO REPLICHE
	
	//ADESSO VEDIAMO L'ENERGIA NEL TEMPO. FACCIAMO UNA MEDIA SULLE FINESTRE LOGARITMICHE IN BASE 2
	printf("Energy vs Time....\n");
	
	//erg[r][k][t]
	
	for(int r = 0; r < NREPLICAS; r++) {

		//file per le energie
		char* filename_erg = (char*)calloc(STRING_SIZE, sizeof(char));
		checkString(filename_erg);
		sprintf(filename_erg, "energy_nrep%d.dat", r);
		FILE* file_erg = fopen(filename_erg, "w");
		checkFile(file_erg, filename_erg);
		free(filename_erg);
		
		fprintf(file_erg, "#time\tenergy\terror\ttemp\n"); 
		
		for(int k = 0; k < NPT-1; k++) {
		
			int b_max = (int)log2(iter);
			
			for(int b = 0; b < b_max; b++) {
				
			
				int t_min = intPow(2, b) - 1;
				int t_max = intPow(2, b + 1) - 2;
				
				double* tempErg = (double*)calloc(t_max-t_min+1, sizeof(double));
				checkDouble(tempErg);
				int l = 0;
				
				for(int t = t_min; t <= t_max; t++) {
					tempErg[l] = erg[r][k][t];
					l++;
				}

				double mean_erg = MeanVec(tempErg, t_max-t_min+1);
				double err_erg = StandardDeviation(tempErg, mean_erg, t_max-t_min+1);
				fprintf(file_erg, "%d\t%le\t%le\t%lf\n", intPow(2,b)*NSTEP, mean_erg, err_erg, temperature[k]);
				free(tempErg);
			}
			
			fprintf(file_erg, "\n");
		}
		
		fclose(file_erg);
		
	}
	
	return erg;
}
*/
void freeMemory5(double***** data, int R, int L, int K, int C) {
	
	//data[R][L][K]
	for(int r = 0; r < R; r++) {
		for(int l = 0; l < L; l++) {
			for(int k = 0; k < K; k++) {
				for(int c = 0; c < C; c++) {
					free(data[r][l][k][c]);
				}
				free(data[r][l][k]);
			}
			free(data[r][l]);
		}
		free(data[r]);
	}//wholeDynamics[r][1:erg/0:acc][k][t]
	
	free(data);

}

void freeMemory5C(complex***** data, int R, int L, int K, int C) {
	
	//data[R][L][K]
	for(int r = 0; r < R; r++) {
		for(int l = 0; l < L; l++) {
			for(int k = 0; k < K; k++) {
				for(int c = 0; c < C; c++) {
					free(data[r][l][k][c]);
				}
				free(data[r][l][k]);
			}
			free(data[r][l]);
		}
		free(data[r]);
	}//wholeDynamics[r][1:erg/0:acc][k][t]
	
	free(data);

}

void freeMemory4(double**** data, int R, int L, int K) {
	
	//data[R][L][K]
	for(int r = 0; r < R; r++) {
		for(int l = 0; l < L; l++) {
			for(int k = 0; k < K; k++) {
				free(data[r][l][k]);
			}
			free(data[r][l]);
		}
		free(data[r]);
	}//wholeDynamics[r][1:erg/0:acc][k][t]
	
	free(data);

}

void freeMemory4C(complex**** data, int R, int L, int K) {
	
	//data[R][L][K]
	for(int r = 0; r < R; r++) {
		for(int l = 0; l < L; l++) {
			for(int k = 0; k < K; k++) {
				free(data[r][l][k]);
			}
			free(data[r][l]);
		}
		free(data[r]);
	}//wholeDynamics[r][1:erg/0:acc][k][t]
	
	free(data);

}

void freeMemory3(double*** data, int R, int L) {
	
	//data[R][L][K]
	for(int r = 0; r < R; r++) {
		for(int l = 0; l < L; l++) {
			free(data[r][l]);
		}
		free(data[r]);
	}
	
	free(data);

}

void freeMemory2(double** data, int R) {
	
	//data[R][L][K]
	for(int r = 0; r < R; r++) {
		free(data[r]);
	}
	free(data);
}

void compute_specific_heat_V3(double*** erg, int iter, double* temp) {

	//erg[r][k][t]
	char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(filename);
		
	sprintf(filename, "specific_heat_N%d.dat", Size);
		
	FILE* sh_file = fopen(filename, "w");
	checkFile(sh_file, filename);
	free(filename);
		
	fprintf(sh_file, "#temp\tspec_heat\terrors\n");
	
	double** sh = (double**)calloc(NPT, sizeof(double*));
	checkMat(sh);
	
	double** err_sh = (double**)calloc(NPT, sizeof(double*));
	checkMat(err_sh);
	
	for(int k = 0; k < NPT - 1; k++) {
		
		sh[k] = (double*)calloc(NREPLICAS, sizeof(double));
		checkDouble(sh[k]);
		err_sh[k] = (double*)calloc(NREPLICAS, sizeof(double));
		checkDouble(err_sh[k]);
	}
	
	for(int r = 0; r < NREPLICAS; r++) {
		
		
		for(int k = 0; k < NPT -1; k++) {
			
			double* temp_erg = (double*)calloc(iter - (int)(iter/2) +1, sizeof(double));
			checkDouble(temp_erg);
			double* erg2 = (double*)calloc(iter- (int)(iter/2) +1, sizeof(double));
			checkDouble(temp_erg);
			int l = 0; 
			for(int t = (int)(iter/2); t < iter; t++) {                             //per il calcolo del calore specifico mi servono solo i dati 
																					// all'equilibrio, in quanto deve valere il teorema fluttuazione-diss.
				temp_erg[l] = erg[r][k][t];											// suppongo che la seconda metà della simulazione sia all'equilibrio.
				erg2[l] = erg[r][k][t] * erg[r][k][t];
				l++;
			}
			
			double e = MeanVec(temp_erg, iter - (int)(iter/2)+1);
			double e2 = MeanVec(erg2, iter- (int)(iter/2)+1);
			double se = StandardDeviation(temp_erg, e, iter - (int)(iter/2) + 1);
			double se2 = StandardDeviation(erg2, e2, iter-(int)(iter/2) +1);
			
			sh[k][r] = (e2 - e*e)/(Size*temp[k]*temp[k]);
			err_sh[k][r] = (1.0/(Size*temp[k]*temp[k])) * sqrt(se2*se2 + 4*e*e*se*se);
			
			//fprintf(sh_file, "%lf\t%le\t%le\t%d\n", temp[k], sh, err_sh, r);
			
			free(temp_erg);
			free(erg2);
		}
		//fprintf(sh_file, "\n");
	}
	
	for(int k = 0; k < NPT - 1; k++) {
		
		double sum = 0.;
		double err = 0.;
		for(int r = 0; r < NREPLICAS; r++) {
			sum += sh[k][r]/NREPLICAS;
			err += err_sh[k][r]/NREPLICAS;
		}
		
		fprintf(sh_file, "%lf\t%le\t%le\n", temp[k], sum, err);
		
	}
	
	freeMemory2(sh, NPT);
	freeMemory2(err_sh, NPT);
	fclose(sh_file);
	

}

complex** read_configuration_file_V3(int r, int iter) {

	complex** conf = (complex**)calloc(NPT, sizeof(complex*));
	checkComplexMat(conf);
	char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(filename);
	sprintf(filename, "config_nrep%d_iter_%d.dat", r, iter);
	FILE* fr = fopen(filename, "r");
	checkFile(fr, filename);
	free(filename);
	for(int k = 0; k < NPT; k++) {
		conf[k] = (complex*)calloc(Size, sizeof(complex));
		checkComplex(conf[k]);
		for(int j = 0; j < Size; j++) {
			double temp1 = 0.;
			int temp2 = 0;
			fscanf(fr, "%lf %d %le %le\n", &temp1, &temp2, &conf[k][j].re, &conf[k][j].im);
			
		}
	}
	fclose(fr);
	
	//printf("%le\t%le\n", conf[NPT-1][Size-1].re, conf[NPT-1][Size-1].im);
	
	return conf;

}

complex**** load_configurations_PT(int* iter, int nFILE) {						//Questa funzione carica tutti i file di configurazione come fossero
																				// un unico blocco (tutti dati all'equilibrio termodinamico)
	complex**** spin = (complex****)calloc(NREPLICAS, sizeof(complex***));
	checkHyp2ComplexMat(spin);
	
	for(int r = 0; r < NREPLICAS; r++) {
		
		spin[r] = (complex***)calloc(nFILE, sizeof(complex**));
		checkHypComplexMat(spin[r]);
		
		for(int i = 0; i < nFILE; i++) {
			
			spin[r][i] = read_configuration_file_V3(r, iter[i]);
		}
	}
	
	complex**** new_spin = (complex****)calloc(NREPLICAS, sizeof(complex***));
	checkHyp2ComplexMat(new_spin);
	
	for(int r = 0; r < NREPLICAS; r++) {
		
		new_spin[r] = (complex***)calloc(NPT, sizeof(complex**));
		checkHypComplexMat(new_spin[r]);
		
		for(int k = 0; k < NPT; k++) {
			
			new_spin[r][k] = (complex**)calloc(Size, sizeof(complex*));
			checkComplexMat(new_spin[r][k]);
			
			for(int j = 0; j < Size; j++) {
				
				new_spin[r][k][j] = (complex*)calloc(nFILE, sizeof(complex));
				checkComplex(new_spin[r][k][j]);
				
				for(int i = 0; i < nFILE; i++) {
					
					new_spin[r][k][j][i].re = spin[r][i][k][j].re;
					new_spin[r][k][j][i].im = spin[r][i][k][j].im;
					
				}
			}
		}
	}
	
	freeMemory4C(spin, NREPLICAS, nFILE, NPT);
	
	return new_spin;
} 

/*
complex***** load_configurations_NOPT(int* iter, int nFILE) {
	
	complex**** spin = (complex****)calloc(NREPLICAS, sizeof(complex***));
	checkHyp2ComplexMat(spin);
	
	int b_max = (int)log2(nFILE);
	
	for(int r = 0; r < NREPLICAS; r++) {
	
		spin[r] = (complex***)calloc(nFILE, sizeof(complex**));
		checkHypComplexMat(spin[r]);
		
		for(int i = 0; i < nFILE; i++) {
			
			spin[r][i] = read_configuration_file_V3(r, iter[i]);
		}
	}
	
	//spin[r][i][k][j] ---- > new_spin[r][k][j][b][i]
	
	complex***** new_spin = (complex*****)calloc(NREPLICAS, sizeof(complex****));
	
	for(int r = 0; r < NREPLICAS; r++) {
		
		new_spin[r] = (complex****)calloc(NPT, sizeof(complex***));
		checkHyp2ComplexMat(new_spin[r]);
		
		for(int k = 0; k < NPT; k++) {
			
			new_spin[r][k] = (complex***)calloc(Size, sizeof(complex**));
			checkHypComplexMat(new_spin[r][k]);
			
			for(int j = 0; j < Size; j++) {
				
				new_spin[r][k][j] = (complex**)calloc(b_max, sizeof(complex*));
				checkComplexMat(new_spin[r][k][j]);
				
				for(int b = 0; b < b_max; b++) {
				
					int t_min = intPow(2, b) - 1;
					int t_max = intPow(2, b + 1) - 2;
					
					new_spin[r][k][j][b] = (complex*)calloc(t_max-t_min+1, sizeof(complex));
					checkComplex(new_spin[r][k][j][b]);
					
					int l = 0;
					for(int t = t_min; t <= t_max; t++) {
						
						new_spin[r][k][j][b][t-t_min].re = spin[r][t][k][j].re;
						new_spin[r][k][j][b][t-t_min].im = spin[r][t][k][j].im;
						l++;
		
					}
				}
			}
		}
	}
	
	freeMemory4C(spin, NREPLICAS, nFILE, NPT);
	//new_spin[r][k][j][b][i]
	
	return new_spin;

}
*/
double**** instant_spectrum_PT(complex**** spin, int nFILE) {

	double**** intensity = (double****)calloc(NREPLICAS, sizeof(double***));
	checkHyp2Mat(intensity);
	
	for(int r = 0; r < NREPLICAS; r++) {
		
		intensity[r] = (double***)calloc(NPT, sizeof(double**));
		checkHypMat(intensity[r]);
		
		for(int k = 0; k < NPT; k++) {
			
			intensity[r][k] = (double**)calloc(Size, sizeof(double*));
			checkMat(intensity[r][k]);
			
			for(int j = 0; j < Size; j++) {
				
				intensity[r][k][j] = (double*)calloc(nFILE, sizeof(double));
				checkDouble(intensity[r][k][j]);
				
				for(int i = 0; i < nFILE; i++) {
					
					intensity[r][k][j][i] = Mod2Value(spin[r][k][j][i]);
				}
			}
		}
	}
	
	return intensity;
}

/*
double***** instant_spectrum_NOPT(complex***** spin, int nFILE) {

	int b_max = (int)log2(nFILE);
	
	double***** intensity = (double*****)calloc(NREPLICAS, sizeof(double****));
	
	for(int r = 0; r < NREPLICAS; r++) {
		
		intensity[r] = (double****)calloc(NPT, sizeof(double***));
		checkHyp2Mat(intensity[r]);
		
		for(int k = 0; k < NPT; k++) {
			
			intensity[r][k] = (double***)calloc(Size, sizeof(double**));
			checkHypMat(intensity[r][k]);
			
			for(int j = 0; j < Size; j++) {
				
				intensity[r][k][j] = (double**)calloc(b_max, sizeof(double*));
				checkMat(intensity[r][k][j]);
				
				for(int b = 0; b < b_max; b++) {
					
					intensity[r][k][j][b] = (double*)calloc(intPow(2, b), sizeof(double));
					checkDouble(intensity[r][k][j][b]);
					
					for(int i = 0; i < intPow(2, b); i++) {
						
						intensity[r][k][j][b][i] = Mod2Value(spin[r][k][j][b][i]);
					}
				}
			}
		}
	}


	return intensity;

}

*/
double*** mean_spectrum_PT(double**** intensity, int nFILE, double* freq, double* temp) {

	double*** spectrum = (double***)calloc(NREPLICAS, sizeof(double));
	checkHypMat(spectrum);
	double*** error = (double***)calloc(NREPLICAS, sizeof(double));
	checkHypMat(error);
	
	printf("Printing spectrum files...\n");
	for(int r = 0; r < NREPLICAS; r++) {
		
		spectrum[r] = (double**)calloc(NPT, sizeof(double*));
		checkMat(spectrum[r]);
		error[r] = (double**)calloc(NPT, sizeof(double*));
		checkMat(error[r]);
		
		for(int k = 0; k < NPT; k++) {
			
			spectrum[r][k] = (double*)calloc(Size, sizeof(double));
			checkDouble(spectrum[r][k]);
			
			error[r][k] = (double*)calloc(Size, sizeof(double));
			checkDouble(error[r][k]);
			for(int j = 0; j < Size; j++) {
				
				spectrum[r][k][j] = MeanVec(intensity[r][k][j], nFILE);
				error[r][k][j] = StandardDeviation(intensity[r][k][j], spectrum[r][k][j], nFILE);
				
				//fprintf(spectrum_file, "%le\t%le\t%le\t%lf\n", freq[j], spectrum[r][k][j]/sqrt(temp[k]), error, temp[k]);	
			}
			
			//fprintf(spectrum_file, "\n");
		}
		//fclose(spectrum_file);
	}
	
	char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(filename);
		
	sprintf(filename, "spectrum_N%d.dat", Size);
		
	FILE* spectrum_file = fopen(filename, "w");
	checkFile(spectrum_file, filename);
	free(filename);
		
	fprintf(spectrum_file, "#frequency\tintensity\terror\ttemperature\n");
	
	for(int k = 0; k < NPT; k++) {
		
		for(int j = 0; j < Size; j++) {
			
			double sum = 0.;
			double err = 0.;
			
			for(int r = 0; r < NREPLICAS; r++) {
				
				sum += spectrum[r][k][j] / NREPLICAS;
				err += error[r][k][j] / NREPLICAS;
					
			}
			
			fprintf(spectrum_file, "%le\t%le\t%le\t%lf\n", freq[j], sum/sqrt(temp[k]), err/sqrt(temp[k]), temp[k]);	
			
		}
		
		fprintf(spectrum_file, "\n");
	}	
	
	fclose(spectrum_file);	
	return spectrum;
	
}

/*
double**** mean_spectrum_NOPT(double***** intensity, int nFILE, double* freq, double* temp) {
	
	
	int b_max = (int)log2(nFILE);
	double**** spectrum = (double****)calloc(NREPLICAS, sizeof(double***));
	checkHyp2Mat(spectrum);
	
	double**** error = (double****)calloc(NREPLICAS, sizeof(double***));
	checkHyp2Mat(error);
	
	for(int r = 0; r < NREPLICAS; r++) {
		
		spectrum[r] = (double***)calloc(NPT, sizeof(double**));
		checkHypMat(spectrum[r]);
		
		error[r] = (double***)calloc(NPT, sizeof(double**));
		checkHypMat(error[r]);
		
		for(int k = 0; k < NPT; k++) {
			
			spectrum[r][k] = (double**)calloc(Size, sizeof(double*));
			checkMat(spectrum[r][k]);
			
			error[r][k] = (double**)calloc(Size, sizeof(double*));
			checkMat(error[r][k]);
			
			for(int j = 0; j < Size; j++) {
				
				spectrum[r][k][j] = (double*)calloc(b_max, sizeof(double));
				checkDouble(spectrum[r][k][j]);
				
				error[r][k][j] = (double*)calloc(b_max, sizeof(double));
				checkDouble(error[r][k][j]);
				
				for(int b = 0; b < b_max; b++) {
					
					spectrum[r][k][j][b] = MeanVec(intensity[r][k][j][b], intPow(2,b));
					error[r][k][j][b] = StandardDeviation(intensity[r][k][j][b], spectrum[r][k][j][b], intPow(2, b));
				}
			}
		}
	}
	
	printf("Printing spectrum files...\n");
	
	for(int b = 0; b < b_max; b++) {
		char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
		checkString(filename);
			
		sprintf(filename, "spectrum_N%d_block%d.dat",Size, b);
			
		FILE* spectrum_file = fopen(filename, "w");
		checkFile(spectrum_file, filename);
		free(filename);
			
		fprintf(spectrum_file, "#frequency\tintensity\terror\ttemperature\n");
		
		for(int k = 0; k < NPT; k++) {
			
			for(int j = 0; j < Size; j++) {
				
				double sum = 0.;
				double err = 0.;
				
				for(int r = 0; r < NREPLICAS; r++) {
					
					sum += spectrum[r][k][j][b]/ NREPLICAS;
					err += error[r][k][j][b] / NREPLICAS;
				}
				
				fprintf(spectrum_file, "%le\t%le\t%le\t%lf\n", freq[j], sum/sqrt(temp[k]), err/sqrt(temp[k]), temp[k]);
			}
			
			fprintf(spectrum_file, "\n");
		}
		
		fclose(spectrum_file);
		
	}
	freeMemory4(error, NREPLICAS, NPT, Size);
	
	return spectrum;
}

*/
void histogram_PT(double** overlap, int nFILE, double* temp, char* type) {

    double min = -1.00;
    double max = 1.00;

    int nq = (int)(NREPLICAS*(NREPLICAS-1)/2)*nFILE;
    
    for(int k = 0; k < NPT; k++) {
    	
    	for(int l = 0; l < nq; l++) {
    		
    		if(overlap[k][l] < min || overlap[k][l] > max) {
    			printf("OUT OF RANGE!\n");
    		}
    	}
    }

    double bin_size = (max - min)/HISTBARS;
	
	// if(vec[j] > min + i*bin_size && vec[j] <= min + (i+1)*bin_size)
	// qMean[i] = min + bin_size/2 + i*bin_size;

	double* q = (double*)calloc(HISTBARS, sizeof(double));
	checkDouble(q);
	
	for(int i = 0; i < HISTBARS; i++) {
		q[i] = min + bin_size/2 + i* bin_size;
	}
	
	char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(filename);
	
	sprintf(filename, "%s_histogram.dat", type);
	
	FILE* histogram_file = fopen(filename, "w");
	checkFile(histogram_file, filename);
	
	free(filename);
	
	fprintf(histogram_file, "#overlap\tdist\ttemperature\n");
	
	for(int k = 0; k < NPT; k++) {
	
		double* p = (double*)calloc(HISTBARS, sizeof(double));
		checkDouble(p);
		
		for(int i = 0; i < HISTBARS; i++) {
			
			for(int l = 0; l < nq; l++) {
				
				if((overlap[k][l] > min + i* bin_size) && (overlap[k][l] <= min + (i+1)*bin_size)) {
					p[i]++;
				} 
			}
		}
		
		double norm = ArraySum(p, HISTBARS) * bin_size;
		
		for(int i = 0; i < HISTBARS; i++) {
			fprintf(histogram_file, "%le\t%le\t%lf\n", q[i], p[i]/norm, temp[k]);
		}
		
		fprintf(histogram_file, "\n");
		free(p);
	}
	
    fclose(histogram_file);
    free(q);
}


void compute_theo_ifo_PT(double**** spectrum, double*** mean_spectrum, int nFILE, double* temp) {
	//riscrivo anche la parte che calcola il singolo elemento di matrice per l'overlap
	
	char* overlap_filename = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(overlap_filename);
	
	sprintf(overlap_filename, "ifo_theo_N%d.dat", Size);
	
	FILE* overlap_file = fopen(overlap_filename, "w");
	checkFile(overlap_file, overlap_filename);
	free(overlap_filename);
	
	fprintf(overlap_file, "#time\tifo\ttemperature\talpha\tbeta\n");
	
	double** ifo = (double**)calloc(NPT, sizeof(double*));
	checkMat(ifo);
	
	
	for(int k = 0; k < NPT; k++) {
		
		ifo[k] = (double*)calloc((int)(NREPLICAS*(NREPLICAS-1)/2) * nFILE, sizeof(double));
		checkDouble(ifo[k]);
		
		int l = 0;
		for(int i = 0; i < nFILE; i++) {
		
			for(int r1 = 0; r1 < NREPLICAS; r1++) {
			
				for(int r2 = 0; r2 < r1; r2++) {
				
					double sum1 = 0.;
					double sum2 = 0.;
					double overlap = 0.;
				
					for(int j = 0; j < Size; j++) {
					
						double delta1 = spectrum[r1][k][j][i] - mean_spectrum[r1][k][j];
						double delta2 = spectrum[r2][k][j][i] - mean_spectrum[r2][k][j];
						
						overlap += delta1*delta2;
						sum1 += delta1*delta1;
						sum2 += delta2*delta2;
						
					}
					
					ifo[k][l] = overlap / (sqrt(sum1 * sum2));
					fprintf(overlap_file, "%d\t%le\t%lf\t%d\t%d\n", i, ifo[k][l], temp[k], r1, r2);
					
					l++;
					
				}
			}
		}
		
		fprintf(overlap_file, "\n");
	}
	
	fclose(overlap_file);
	
	printf("Making histograms...\n");
	char* type = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(type);
	sprintf(type, "theo_ifo");
	
	histogram_PT(ifo, nFILE, temp, type);
	
	free(type);
	freeMemory2(ifo, NPT);
	
}

void histogram_NOPT(double** overlap, int nFILE, double* temp, char* type, int b) {

    double min = -1.00;
    double max = 1.00;

    int nq = (int)(NREPLICAS*(NREPLICAS-1)/2)*nFILE;
    
    for(int k = 0; k < NPT; k++) {
    	
    	for(int l = 0; l < nq; l++) {
    		
    		if(overlap[k][l] < min || overlap[k][l] > max) {
    			printf("OUT OF RANGE!\n");
    		}
    	}
    }

    double bin_size = (max - min)/HISTBARS;
	
	// if(vec[j] > min + i*bin_size && vec[j] <= min + (i+1)*bin_size)
	// qMean[i] = min + bin_size/2 + i*bin_size;

	double* q = (double*)calloc(HISTBARS, sizeof(double));
	checkDouble(q);
	
	for(int i = 0; i < HISTBARS; i++) {
		q[i] = min + bin_size/2 + i* bin_size;
	}
	
	char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(filename);
	
	sprintf(filename, "%s_histogram_block%d.dat", type, b);
	
	FILE* histogram_file = fopen(filename, "w");
	checkFile(histogram_file, filename);
	
	free(filename);
	
	fprintf(histogram_file, "#overlap\tdist\ttemperature\n");
	
	for(int k = 0; k < NPT; k++) {
	
		double* p = (double*)calloc(HISTBARS, sizeof(double));
		checkDouble(p);
		
		for(int i = 0; i < HISTBARS; i++) {
			
			for(int l = 0; l < nq; l++) {
				
				if((overlap[k][l] > min + i* bin_size) && (overlap[k][l] <= min + (i+1)*bin_size)) {
					p[i]++;
				} 
			}
		}
		
		double norm = ArraySum(p, HISTBARS) * bin_size;
		
		for(int i = 0; i < HISTBARS; i++) {
			fprintf(histogram_file, "%le\t%le\t%lf\n", q[i], p[i]/norm, temp[k]);
		}
		
		fprintf(histogram_file, "\n");
		free(p);
	}
	
    fclose(histogram_file);
    free(q);
}



/*
void compute_theo_ifo_NOPT(double***** spectrum, double**** mean_spectrum, int nFILE, double* temp) {


	int b_max = (int)log2(nFILE);
	double*** ifo = (double***)calloc(b_max, sizeof(double**));
	checkHypMat(ifo);
	
	
	for(int b = 0; b < b_max; b++) {
	
		char* overlap_filename = (char*)calloc(STRING_SIZE, sizeof(char));
		checkString(overlap_filename);
	
		sprintf(overlap_filename, "ifo_theo_N%d_block%d.dat", Size, b);
	
		FILE* overlap_file = fopen(overlap_filename, "w");
		checkFile(overlap_file, overlap_filename);
		free(overlap_filename);
	
		fprintf(overlap_file, "#time\tifo\ttemperature\talpha\tbeta\n");

		ifo[b] = (double**)calloc(NPT, sizeof(double*));
		checkMat(ifo[b]);
		
		for(int k = 0; k < NPT; k++) {
			
			ifo[b][k] = (double*)calloc((int)(NREPLICAS*(NREPLICAS-1)/2) * intPow(2, b), sizeof(double));
			checkDouble(ifo[b][k]);
			
			int l = 0;
			
			for(int i = 0; i < intPow(2, b); i++) {
				
				for(int r1 = 0; r1 < NREPLICAS; r1++) {
					
					for(int r2 = 0; r2 < r1; r2++) {
					
						double sum1 = 0.;
						double sum2 = 0.;
						double overlap = 0.;
						
						for(int j = 0; j < Size; j++) {
						
							double delta1 = spectrum[r1][k][j][b][i] - mean_spectrum[r1][k][j][b];
							double delta2 = spectrum[r2][k][j][b][i] - mean_spectrum[r2][k][j][b];
						
							overlap += delta1*delta2;
							sum1 += delta1*delta1;
							sum2 += delta2*delta2;
							
						}
						
						ifo[b][k][l] = overlap / (sqrt(sum1 * sum2));
						fprintf(overlap_file, "%d\t%le\t%lf\t%d\t%d\n", i, ifo[b][k][l], temp[k], r1, r2);
						l++;
					}
				}
			}
			
			fprintf(overlap_file, "\n");
		}
		
		fclose(overlap_file);
	}
	
	char* type = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(type);
	sprintf(type, "theo_ifo");
	
	for(int b = 1; b < b_max; b++) {
		histogram_NOPT(ifo[b], intPow(2,b), temp, type, b);
	}
	free(type);

	freeMemory3(ifo, b_max, NPT);

}

*/
void compute_exp_ifo_PT(double*** mean_spectrum, int nFILE, double* temp) {
	
	
	//mean_spectrum[r][k][j]
	char* overlap_filename = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(overlap_filename);
	
	
	double** mean_replica = (double**)calloc(NPT, sizeof(double*));
	checkMat(mean_replica);
	
	
	for(int k = 0; k < NPT; k++) {
		
		mean_replica[k] = (double*)calloc(Size, sizeof(double));
		checkDouble(mean_replica[k]);
		
		for(int j = 0; j < Size; j++) {
			
			double sum = 0.;
			
			for(int r = 0; r < NREPLICAS; r++) {
				
				sum += mean_spectrum[r][k][j];
			}
			
			mean_replica[k][j] = sum/NREPLICAS;
		}
	}
	
	
	
	sprintf(overlap_filename, "ifo_exp_N%d.dat", Size);
	
	FILE* overlap_file = fopen(overlap_filename, "w");
	checkFile(overlap_file, overlap_filename);
	free(overlap_filename);
	
	fprintf(overlap_file, "#time\tifo\ttemperature\talpha\tbeta\n");
	
	double** ifo = (double**)calloc(NPT, sizeof(double*));
	checkMat(ifo);
	
	
	for(int k = 0; k < NPT; k++) {
		
		ifo[k] = (double*)calloc((int)(NREPLICAS*(NREPLICAS-1)/2), sizeof(double));
		checkDouble(ifo[k]);
		
		int l = 0;
		
		for(int r1 = 0; r1 < NREPLICAS; r1++) {
			
			for(int r2 = 0; r2 < r1; r2++) {
				
				ifo[k][l] = OverlapIFO_Element(mean_spectrum[r1][k], mean_spectrum[r2][k], mean_replica[k], mean_replica[k], Size);
				fprintf(overlap_file, "%d\t%le\t%lf\t%d\t%d\n", l, ifo[k][l], temp[k], r1, r2);	
				l++;
				
			}
		}
		
		
		fprintf(overlap_file, "\n");
	}
	
	fclose(overlap_file);
	
	printf("Making histograms...\n");
	char* type = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(type);
	sprintf(type, "exp_ifo");
	
	histogram_PT(ifo, 1, temp, type);
	
	free(type);
	freeMemory2(ifo, NPT);
	
	
}


void compute_parisi_overlap_PT(complex**** spin, int nFILE, double* temp) {
	
	
	double** overlap = (double**)calloc(NPT, sizeof(double*));
	checkMat(overlap);
	
	char* overlap_filename = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(overlap_filename);
	sprintf(overlap_filename, "parisi_overlap_N%d.dat", Size);
	
	FILE* overlap_file = fopen(overlap_filename, "w");
	checkFile(overlap_file, overlap_filename);
	free(overlap_filename);
	
	fprintf(overlap_file, "#time\toverlap\ttemperature\talpha\tbeta\n");
	
	for(int k = 0; k < NPT; k++) {
		
		overlap[k] = (double*)calloc((int)(NREPLICAS*(NREPLICAS-1)/2)*nFILE, sizeof(double));
		checkDouble(overlap[k]);
		
		int l = 0;
		for(int i = 0; i < nFILE; i++) {
			
			for(int r1 = 0; r1 < NREPLICAS; r1++) {
				
				for(int r2 = 0; r2 < r1; r2++) {
					
					double sum = 0.;
					
					for(int j = 0; j < Size; j++) {
						
						sum += (spin[r1][k][j][i].re * spin[r2][k][j][i].re + spin[r1][k][j][i].im * spin[r2][k][j][i].im) / Size;
					}
					
					overlap[k][l] = sum;
					fprintf(overlap_file, "%d\t%le\t%lf\t%d\t%d\n", i, overlap[k][l], temp[k], r1, r2);	
					l++;
					
				}
			}
		}
		fprintf(overlap_file, "\n");
	}


	fclose(overlap_file);
	
	printf("Making histograms...\n");
	char* type = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(type);
	sprintf(type, "parisi");
	
	histogram_PT(overlap, nFILE, temp, type);
	
	free(type);
	freeMemory2(overlap, NPT);

}

/*

void compute_parisi_overlap_NOPT(complex***** spin, int nFILE, double* temp) {
	
	int b_max = (int)log2(nFILE);
	double*** overlap = (double***)calloc(b_max, sizeof(double**));
	checkHypMat(overlap);
	
	for(int b = 0; b < b_max; b++) {
	
		overlap[b] = (double**)calloc(NPT, sizeof(double*));
		checkMat(overlap[b]);
	
		char* overlap_filename = (char*)calloc(STRING_SIZE, sizeof(char));
		checkString(overlap_filename);
		sprintf(overlap_filename, "parisi_overlap_N%d.dat", Size);
	
		FILE* overlap_file = fopen(overlap_filename, "w");
		checkFile(overlap_file, overlap_filename);
		free(overlap_filename);
	
		fprintf(overlap_file, "#time\toverlap\ttemperature\talpha\tbeta\n");
	
		for(int k = 0; k < NPT; k++) {
		
			overlap[b][k] = (double*)calloc((int)(NREPLICAS*(NREPLICAS-1)/2)*intPow(2, b), sizeof(double));
			checkDouble(overlap[b][k]);
		
			int l = 0;
			for(int i = 0; i < intPow(2, b); i++) {
			
				for(int r1 = 0; r1 < NREPLICAS; r1++) {
			
					for(int r2 = 0; r2 < r1; r2++) {
					
						double sum = 0.;
					
						for(int j = 0; j < Size; j++) {
						
							sum += (spin[r1][k][j][b][i].re * spin[r2][k][j][b][i].re + spin[r1][k][j][b][i].im * spin[r2][k][j][b][i].im) / Size;
						}
					
						overlap[b][k][l] = sum;
						fprintf(overlap_file, "%d\t%le\t%lf\t%d\t%d\n", i, overlap[b][k][l], temp[k], r1, r2);	
						l++;
					
					}
				}
			}
			fprintf(overlap_file, "\n");
		}
		
		fclose(overlap_file);
	}

	
	
	printf("Making histograms...\n");
	char* type = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(type);
	sprintf(type, "parisi");
	
	for(int b = 0; b < b_max; b++) {
		histogram_NOPT(overlap[b], intPow(2, b), temp, type, b);
	}
	
	free(type);
	freeMemory3(overlap, b_max, NPT);

}



void compute_exp_ifo_NOPT(double**** mean_spectrum, int nFILE, double* temp) {
	
	
	//mean_spectrum[r][k][j][b]
	int b_max = (int)log2(nFILE);
	
	double*** mean_replica = (double***)calloc(b_max, sizeof(double**));
	checkHypMat(mean_replica);
	
	for(int b = 0; b < b_max; b++) {
		
		mean_replica[b] = (double**)calloc(NPT, sizeof(double*));
		checkMat(mean_replica[b]);
		
		for(int k = 0; k < NPT; k++) {
			
			mean_replica[b][k] = (double*)calloc(Size, sizeof(double));
			checkDouble(mean_replica[b][k]);
			
			for(int j = 0; j < Size; j++) {
				
				for(int r = 0; r < NREPLICAS; r++) {
					
					mean_replica[b][k][j] += mean_spectrum[r][k][j][b] / NREPLICAS;
				}
			}
		}
	}
	
	double*** ifo = (double***)calloc(b_max, sizeof(double**));
	checkHypMat(ifo);
	
	for(int b = 0; b < b_max; b++) {
		
		char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
		checkString(filename);
		
		sprintf(filename, "ifo_exp_N%d_block%d.dat", Size, b);
		ifo[b] = (double**)calloc(NPT, sizeof(double*));
		checkMat(ifo[b]);
		
		FILE* overlap_file = fopen(filename, "w");
		checkFile(overlap_file, filename);
		free(filename);
		
		fprintf(overlap_file, "#time\tifo\ttemperature\talpha\tbeta\n");
		
		for(int k = 0; k < NPT; k++) {
			
			ifo[b][k] = (double*)calloc((int)(NREPLICAS * (NREPLICAS - 1)/2), sizeof(double));
			checkDouble(ifo[b][k]);
			
			int l = 0;
			
			for(int r1 = 0; r1 < NREPLICAS; r1++) {
				
				for(int r2 = 0; r2 < r1; r2++) {
					
					double overlap = 0.;
					double sum1 = 0.;
					double sum2 = 0.;
					
					for(int j = 0; j < Size; j++) {
						
						double delta1 = mean_spectrum[r1][k][j][b] - mean_replica[b][k][j];
						double delta2 = mean_spectrum[r2][k][j][b] - mean_replica[b][k][j];
						
						overlap += delta1*delta2;
						sum1 += delta1*delta1;
						sum2 += delta2*delta2;
						
					}
					
					ifo[b][k][l] = overlap/sqrt(sum1*sum2);
					fprintf(overlap_file, "%d\t%le\t%lf\t%d\t%d\n", l, ifo[b][k][l], temp[k], r1, r2);	
					
					l++;
				}
			}
			
			fprintf(overlap_file, "\n");
		}
		
		fclose(overlap_file);
	}
	
	printf("Making histograms...\n");
	char* type = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(type);
	sprintf(type, "exp_ifo");
	
	for(int b = 0; b < b_max; b++) {
		histogram_NOPT(ifo[b], 1, temp, type, b);
	}
	
	free(type);
	freeMemory3(ifo, b_max, NPT);
	
}

*/
double KLD(double* p1, double* p2) {
	
	double d = 0.;
	
	for(int i = 0; i < HISTBARS; i++) {
		
		
		if(p1[i] != 0 && p2[i] != 0) {
			d += p1[i]*log(p1[i]/p2[i]);
		}else{
			d += 0;
		}
		//printf("KLD:\t%le\n", d);
	}
	
	//printf("KLD:\t%le\n", d);
	return d;
	
}

dist* ReadDist(char* type) {

	char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(filename);

	sprintf(filename, "%s_dis_ave_N%d.dat", type, Size);

	FILE* fr = fopen(filename, "r");
	checkFile(fr, filename);
	free(filename);

	skipLine(fr);

	dist* overlap = (dist*)calloc(NPT, sizeof(dist));


	for(int k = 0; k < NPT; k++) {

		overlap[k].p = (double*)calloc(HISTBARS, sizeof(double));
		checkDouble(overlap[k].p);
		overlap[k].q = (double*)calloc(HISTBARS, sizeof(double));
		checkDouble(overlap[k].q);
		overlap[k].e = (double*)calloc(HISTBARS, sizeof(double));
		checkDouble(overlap[k].e);

		for(int i = 0; i < HISTBARS; i++) {
			double t_temp = 0.;
			fscanf(fr, "%le\t%le\t%le\t%lf\n", &overlap[k].q[i], &overlap[k].p[i], &overlap[k].e[i], &t_temp);
			//printf("%le\t%le\t%le\t%lf\n", overlap[k].q[i], overlap[k].p[i], overlap[k].e[i], t_temp);
			if(i == 0) {
				overlap[k].T = t_temp;
			}
		}		
	}

	return overlap;
}

dist* ReadDistBlock(char* type, int b) {

	char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(filename);

	sprintf(filename, "%s_dis_ave_N%d_block%d.dat", type, Size, b);

	FILE* fr = fopen(filename, "r");
	checkFile(fr, filename);
	free(filename);

	skipLine(fr);

	dist* overlap = (dist*)calloc(NPT, sizeof(dist));


	for(int k = 0; k < NPT; k++) {

		overlap[k].p = (double*)calloc(HISTBARS, sizeof(double));
		checkDouble(overlap[k].p);
		overlap[k].q = (double*)calloc(HISTBARS, sizeof(double));
		checkDouble(overlap[k].q);
		overlap[k].e = (double*)calloc(HISTBARS, sizeof(double));
		checkDouble(overlap[k].e);

		for(int i = 0; i < HISTBARS; i++) {
			double t_temp = 0.;
			fscanf(fr, "%le\t%le\t%le\t%lf\n", &overlap[k].q[i], &overlap[k].p[i], &overlap[k].e[i], &t_temp);
			//printf("%le\t%le\t%le\t%lf\n", overlap[k].q[i], overlap[k].p[i], overlap[k].e[i], t_temp);
			if(i == 0) {
				overlap[k].T = t_temp;
			}
		}		
	}

	return overlap;
}


double errKLD(double* p1, double* p2, double* e1, double*e2) {
	double err = 0.;

	for(int i = 0; i < HISTBARS; i++) {

		if(p1[i] != 0 && p2[i] != 0) {
			err += (log(p1[i]/p2[i]) + 1)*(log(p1[i]/p2[i]) + 1)* e1[i]*e1[i] + (p1[i]*p1[i])/(p2[i]*p2[i])*e2[i]*e2[i];
		}

	}
	return sqrt(err);
}

void compute_KLD_PT(const char* name1, const char* name2) {

	char* filename1 = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(filename1);
	char* filename2 = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(filename2);
	printf("KULLBACK-LEIBLER DIVERGENCE FOR %s VS %s.\nReading distribution Files...\n", name1, name2);

	sprintf(filename1, "%s", name1);
	sprintf(filename2, "%s", name2);

	
	dist* overlap1 = ReadDist(filename1);
	dist* overlap2 = ReadDist(filename2);
	
	char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(filename);

	sprintf(filename, "KLD_%s_%s_N%d.dat", filename1, filename2, Size);
	FILE* fw = fopen(filename, "w");
	checkFile(fw, filename);

	fprintf(fw, "#temperature\tKLD\terror\n");

	for(int k = 0; k < NPT; k++) {

		double divergence = KLD(overlap1[k].p, overlap2[k].p);
		double err = errKLD(overlap1[k].p, overlap2[k].p, overlap1[k].e, overlap2[k].e);
		fprintf(fw, "%le\t%le\t%le\n", overlap1[k].T, divergence, err);

	}
	fclose(fw);



}

void compute_KLD_NOPT(const char* name1, const char* name2, int b) {

	char* filename1 = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(filename1);
	char* filename2 = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(filename2);
	printf("KULLBACK-LEIBLER DIVERGENCE FOR %s VS %s.\nReading distribution Files...\n", name1, name2);

	sprintf(filename1, "%s", name1);
	sprintf(filename2, "%s", name2);

	
	dist* overlap1 = ReadDistBlock(filename1, b);
	dist* overlap2 = ReadDistBlock(filename2, b);
	
	char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(filename);

	sprintf(filename, "KLD_%s_%s_N%d_block%d.dat", filename1, filename2, Size, b);
	FILE* fw = fopen(filename, "w");
	checkFile(fw, filename);

	fprintf(fw, "#temperature\tKLD\terror\n");

	for(int k = 0; k < NPT; k++) {

		double divergence = KLD(overlap1[k].p, overlap2[k].p);
		double err = errKLD(overlap1[k].p, overlap2[k].p, overlap1[k].e, overlap2[k].e);
		fprintf(fw, "%le\t%le\t%le\n", overlap1[k].T, divergence, err);

	}
	fclose(fw);



}

void disorder_average_specific_heat(int nsamples) {
	
	
	double** sh = (double**)calloc(NPT, sizeof(double*));
	checkMat(sh);
	double** err_sh = (double**)calloc(NPT, sizeof(double*));
	checkMat(err_sh);
	
	for(int k = 0; k < NPT-1; k++) {
		
		sh[k] = (double*)calloc(nsamples, sizeof(double));
		checkDouble(sh[k]);
		err_sh[k] = (double*)calloc(nsamples, sizeof(double));
		checkDouble(err_sh[k]);
		
	}
	
	double* temp = (double*)calloc(NPT-1, sizeof(double));
	checkDouble(temp);
	
	for(int i = 0; i < nsamples; i++) {
	
		char* filename = (char*)calloc(2*STRING_SIZE, sizeof(char));
		checkString(filename);
		
		sprintf(filename, "sample%d/specific_heat_N%d.dat", i+1, Size);
		
		FILE* fr = fopen(filename, "r");
		checkFile(fr, filename);
		free(filename);
		
		fscanf(fr, "%*[^\n]\n");
		
		for(int k = 0; k < NPT-1; k++) {
			double s = 0.;
			double e = 0.;
			double t = 0.;
			
			fscanf(fr, "%lf\t%le\t%le\n", &t, &s, &e);
			
			if(i == 0) {
				temp[k] = t;
			}
			
			sh[k][i] = s;
			err_sh[k][i] = e;
		}
		
		fclose(fr);
		
	
	}
	
	char* filename_sh = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(filename_sh);
	
	sprintf(filename_sh, "specific_heat_dis_ave_N%d.dat", Size);
	
	FILE* file_sh = fopen(filename_sh, "w");
	checkFile(file_sh, filename_sh);
	free(filename_sh);
	
	fprintf(file_sh, "#temp\tspec_heat\terror\n");
	
	for(int k = 0; k < NPT-1; k++) {
		
		double s = MeanVec(sh[k], nsamples);
		double e = MeanVec(err_sh[k], nsamples);
		
		fprintf(file_sh, "%lf\t%le\t%le\n", temp[k], s, e);
		
	}
	
	fclose(file_sh);
	free(temp);
	freeMemory2(sh, NPT);
	freeMemory2(err_sh, NPT);
	
}

void disorder_average_spectrum_PT(int nsamples) {
	
	double** freq = (double**)calloc(nsamples, sizeof(double*));
	checkMat(freq);
	
	double*** spectrum = (double***)calloc(NPT, sizeof(double**));
	checkHypMat(spectrum);
	
	double*** error = (double***)calloc(NPT, sizeof(double**));
	checkHypMat(error);
	
	double* temp = (double*)calloc(NPT, sizeof(double));
	checkDouble(temp);
	
	for(int k = 0; k < NPT; k++) {
		
		spectrum[k] = (double**)calloc(nsamples, sizeof(double*));
		checkMat(spectrum[k]);
		
		error[k] = (double**)calloc(nsamples, sizeof(double*));
		checkMat(error[k]);
		
		for(int i = 0; i < nsamples; i++) {
			
			spectrum[k][i] = (double*)calloc(Size, sizeof(double));
			checkDouble(spectrum[k][i]);
			
			error[k][i] = (double*)calloc(Size, sizeof(double));
			checkDouble(error[k][i]);
		}
	}
	
	printf("Acquiring spectrum data...\n");
	
	for(int i = 0; i < nsamples; i++) {
		
		freq[i] = (double*)calloc(Size, sizeof(double));
		checkDouble(freq[i]);
		
		char* filename = (char*)calloc(2*STRING_SIZE, sizeof(char));
		checkString(filename);
		
		sprintf(filename, "sample%d/spectrum_N%d.dat", i+1, Size);
		
		FILE* fr = fopen(filename, "r");
		checkFile(fr, filename);
		free(filename);
		
		fscanf(fr, "%*[^\n]\n");
		
		for(int k = 0; k < NPT; k++) {
			
			for(int j = 0; j < Size; j++) {
				
				double f = 0.;
				double s = 0.;
				double e = 0.;
				double t = 0.;
				
				fscanf(fr, "%le\t%le\t%le\t%lf\n", &f, &s, &e, &t);
				
				if(i == 0 && j == 0) {
					temp[k] = t;
				}
				if(k == 0) {
					freq[i][j] = f;
				}
				
				spectrum[k][i][j] = s;
				error[k][i][j] = e;
			}
		}
		
		fclose(fr);
	
	}
	
	double* all_freq = (double*)calloc(nsamples*Size, sizeof(double));
	double** all_spec = (double**)calloc(NPT, sizeof(double*));
	double** all_err = (double**)calloc(NPT, sizeof(double*));
	checkMat(all_spec);
	checkMat(all_err);
	
	printf("Ordering frequencies...\n");
	
	for(int k = 0; k < NPT; k++) {
		
		all_spec[k] = (double*)calloc(nsamples*Size, sizeof(double));
		all_err[k] = (double*)calloc(nsamples*Size, sizeof(double));
		checkDouble(all_spec[k]);
		checkDouble(all_err[k]);
		
	}
	
	
	int l = 0;
	for(int i = 0; i < nsamples; i++) {
		
		for(int j = 0; j < Size; j++) {
			
			all_freq[l] = freq[i][j];
			l++;
		}
	}

	for(int k = 0; k < NPT; k++) {
		
		l = 0;
		
		for(int i = 0; i < nsamples; i++) {
			
			for(int j = 0; j < Size; j++) {
				
				all_spec[k][l] = spectrum[k][i][j];
				all_err[k][l] = error[k][i][j];
				l++;
				
			}
		}
	}
	
	freeMemory3(spectrum, NPT, nsamples);
	freeMemory3(error, NPT, nsamples);
	
	double a = 0.;
	double b = 0.;
	double c = 0.;
	
	for(int i = 0; i < nsamples*Size; i++) {
		
		for(int j = i+1; j < nsamples*Size; j++) {
			
			if(all_freq[i]>all_freq[j]) {
				a = all_freq[i];
				all_freq[i] = all_freq[j];
				all_freq[j] = a;
				
				for(int k = 0; k < NPT; k++) {
					
					b = all_spec[k][i];
					all_spec[k][i] = all_spec[k][j];
					all_spec[k][j] = b;
					
					c = all_err[k][i];
					all_err[k][i] = all_err[k][j];
					all_err[k][j] = c;
				}
				
			}
		}
	}
	
	char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(filename);
	sprintf(filename, "spectrum_all_N%d.dat", Size);
	
	FILE* fw = fopen(filename, "w");
	checkFile(fw, filename);
	free(filename);
	
	fprintf(fw, "#frequencies\tintensities\terror\ttemperature\n");
	
	for(int k = 0; k < NPT; k++) {
		
		for(int i = 0; i < nsamples*Size; i++) {
			
			fprintf(fw, "%le\t%le\t%le\t%lf\n", all_freq[i], all_spec[k][i]/temp[k], all_err[k][i]/temp[k], temp[k]);
		}
		
		fprintf(fw, "\n");
	}
	
	fclose(fw);
	
	printf("Rebinning spectrum...\n");
	
	double min = 0.;
	double max = 1.;
	
	double bin_size = (max - min)/HISTBARS; 
	
	double* f = (double*)calloc(HISTBARS, sizeof(double));
	checkDouble(f);
	
	for(int i = 0; i < HISTBARS; i++) {
		f[i] = min + bin_size/2 + i* bin_size;
	}
	
	int** counter = (int**)calloc(NPT, sizeof(int*));
	//all_spec[k][i]
	for(int k = 0; k < NPT; k++) {
		
		counter[k] = (int*)calloc(HISTBARS, sizeof(int));
		checkInt(counter[k]);
		
		for(int j = 0; j < HISTBARS; j++) {
			
			for(int i = 0; i < nsamples*Size; i++) {
				
				if((all_freq[i] > min + j* bin_size) && (all_freq[i] <= min + (j+1)*bin_size)) {
					
					counter[k][j]++;
				}
			}
		}
	}
	
	char* reb_filename = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(reb_filename);
	
	sprintf(reb_filename, "spectrum_rebinned_N%d.dat", Size);
	FILE* reb_file = fopen(reb_filename, "w");
	checkFile(reb_file, reb_filename);
	free(reb_filename);
	
	fprintf(reb_file, "#frequencies\tintensities\terror\ttemperature\n");
	for(int k = 0; k < NPT; k++) {
		
		int i = 0;
		for(int j = 0; j < HISTBARS; j++) {
			
			double s = 0.;
			double e = 0.;
			
			for(int u = 0; u < counter[k][j]; u++) {
				
				s += all_spec[k][i]/counter[k][j];
				e += all_err[k][i]/counter[k][j];
				
				i++;
			}
			
			fprintf(reb_file, "%le\t%le\t%le\t%lf\n", f[j], s, e, temp[k]);
			
		}
		
		fprintf(reb_file, "\n");
	}
	
	fclose(reb_file);
	free(all_freq);
	freeMemory2(all_spec, NPT);
	freeMemory2(all_err, NPT);
	free(temp);
	
} 

void disorder_average_theo_ifo_PT(int nsamples) {
	
	double* q = (double*)calloc(HISTBARS, sizeof(double));
	checkDouble(q);
	
	double* temp = (double*)calloc(NPT, sizeof(double));
	checkDouble(temp);
	
	double*** p = (double***)calloc(nsamples, sizeof(double**));
	checkHypMat(p);
	
	
	
	for(int n = 0; n < nsamples; n++) {
		
		p[n] = (double**)calloc(NPT, sizeof(double*));
		checkMat(p[n]);
		
		char* filename = (char*)calloc(2*STRING_SIZE, sizeof(char));
		checkString(filename);
		
		sprintf(filename, "sample%d/theo_ifo_histogram.dat", n+1);
		
		FILE* fr = fopen(filename, "r");
		checkFile(fr, filename);
		free(filename);
		
		fscanf(fr, "%*[^\n]\n");
		
		for(int k = 0; k < NPT; k++) {
			
			p[n][k] = (double*)calloc(HISTBARS, sizeof(double));
			checkDouble(p[n][k]);
			
			for(int i = 0; i < HISTBARS; i++) {
				
				double ov = 0.;
				double pr = 0.;
				double t = 0.;
				
				fscanf(fr, "%le\t%le\t%lf\n", &ov, &pr, &t);
				
				if(n == 0 && i == 0) {
					temp[k] = t;
				}
				if(n == 0 && k == 0) {
					q[i] = ov;
				}
				
				p[n][k][i] = pr;
				
			}
		}
		
		fclose(fr);
	
	}
	
	double*** new_pr = (double***)calloc(NPT, sizeof(double**));
	checkHypMat(new_pr);
	
	char* ifo_filename = (char*)calloc(STRING_SIZE, sizeof(char*));
	checkString(ifo_filename);
	sprintf(ifo_filename, "theo_ifo_dis_ave_N%d.dat", Size);
	
	FILE* file_ifo = fopen(ifo_filename, "w");
	checkFile(file_ifo, ifo_filename);
	free(ifo_filename);
	
	fprintf(file_ifo, "#overlap\tdist\terror\ttemperature\n");
	
	for(int k = 0; k < NPT; k++) {
		
		new_pr[k] = (double**)calloc(HISTBARS, sizeof(double*));
		checkMat(new_pr[k]);
		
		for(int i = 0; i < HISTBARS; i++) {
			
			new_pr[k][i] = (double*)calloc(nsamples, sizeof(double));
			checkDouble(new_pr[k][i]);
			
			for(int n = 0; n < nsamples; n++) {
				
				new_pr[k][i][n] = p[n][k][i];
			}
			
			double mean = MeanVec(new_pr[k][i], nsamples);
			double err = StandardDeviation(new_pr[k][i], mean, nsamples);
			
			fprintf(file_ifo, "%le\t%le\t%le\t%lf\n", q[i], mean, err, temp[k]);
			
		}
		
		fprintf(file_ifo, "\n");
	}
	
	fclose(file_ifo);
	
	freeMemory3(p, nsamples, NPT);
	freeMemory3(new_pr, NPT, HISTBARS);
	free(q);
	free(temp);
	

}


void disorder_average_exp_ifo_PT(int nsamples) {

	double* q = (double*)calloc(HISTBARS, sizeof(double));
	checkDouble(q);
	
	double* temp = (double*)calloc(NPT, sizeof(double));
	checkDouble(temp);
	
	double*** p = (double***)calloc(nsamples, sizeof(double**));
	checkHypMat(p);
	
	
	
	for(int n = 0; n < nsamples; n++) {
		
		p[n] = (double**)calloc(NPT, sizeof(double*));
		checkMat(p[n]);
		
		char* filename = (char*)calloc(2*STRING_SIZE, sizeof(char));
		checkString(filename);
		
		sprintf(filename, "sample%d/exp_ifo_histogram.dat", n+1);
		
		FILE* fr = fopen(filename, "r");
		checkFile(fr, filename);
		free(filename);
		
		fscanf(fr, "%*[^\n]\n");
		
		for(int k = 0; k < NPT; k++) {
			
			p[n][k] = (double*)calloc(HISTBARS, sizeof(double));
			checkDouble(p[n][k]);
			
			for(int i = 0; i < HISTBARS; i++) {
				
				double ov = 0.;
				double pr = 0.;
				double t = 0.;
				
				fscanf(fr, "%le\t%le\t%lf\n", &ov, &pr, &t);
				
				if(n == 0 && i == 0) {
					temp[k] = t;
				}
				if(n == 0 && k == 0) {
					q[i] = ov;
				}
				
				p[n][k][i] = pr;
				
			}
		}
		
		fclose(fr);
	
	}
	
	double*** new_pr = (double***)calloc(NPT, sizeof(double**));
	checkHypMat(new_pr);
	
	char* ifo_filename = (char*)calloc(STRING_SIZE, sizeof(char*));
	checkString(ifo_filename);
	sprintf(ifo_filename, "exp_ifo_dis_ave_N%d.dat", Size);
	
	FILE* file_ifo = fopen(ifo_filename, "w");
	checkFile(file_ifo, ifo_filename);
	free(ifo_filename);
	
	fprintf(file_ifo, "#overlap\tdist\terror\ttemperature\n");
	
	for(int k = 0; k < NPT; k++) {
		
		new_pr[k] = (double**)calloc(HISTBARS, sizeof(double*));
		checkMat(new_pr[k]);
		
		for(int i = 0; i < HISTBARS; i++) {
			
			new_pr[k][i] = (double*)calloc(nsamples, sizeof(double));
			checkDouble(new_pr[k][i]);
			
			for(int n = 0; n < nsamples; n++) {
				
				new_pr[k][i][n] = p[n][k][i];
			}
			
			double mean = MeanVec(new_pr[k][i], nsamples);
			double err = StandardDeviation(new_pr[k][i], mean, nsamples);
			
			fprintf(file_ifo, "%le\t%le\t%le\t%lf\n", q[i], mean, err, temp[k]);
			
		}
		
		fprintf(file_ifo, "\n");
	}
	
	fclose(file_ifo);
	
	freeMemory3(p, nsamples, NPT);
	freeMemory3(new_pr, NPT, HISTBARS);
	free(q);
	free(temp);
	


}


void disorder_average_parisi_PT(int nsamples) {

	double* q = (double*)calloc(HISTBARS, sizeof(double));
	checkDouble(q);
	
	double* temp = (double*)calloc(NPT, sizeof(double));
	checkDouble(temp);
	
	double*** p = (double***)calloc(nsamples, sizeof(double**));
	checkHypMat(p);
	
	
	
	for(int n = 0; n < nsamples; n++) {
		
		p[n] = (double**)calloc(NPT, sizeof(double*));
		checkMat(p[n]);
		
		char* filename = (char*)calloc(2*STRING_SIZE, sizeof(char));
		checkString(filename);
		
		sprintf(filename, "sample%d/parisi_histogram.dat", n+1);
		
		FILE* fr = fopen(filename, "r");
		checkFile(fr, filename);
		free(filename);
		
		fscanf(fr, "%*[^\n]\n");
		
		for(int k = 0; k < NPT; k++) {
			
			p[n][k] = (double*)calloc(HISTBARS, sizeof(double));
			checkDouble(p[n][k]);
			
			for(int i = 0; i < HISTBARS; i++) {
				
				double ov = 0.;
				double pr = 0.;
				double t = 0.;
				
				fscanf(fr, "%le\t%le\t%lf\n", &ov, &pr, &t);
				
				if(n == 0 && i == 0) {
					temp[k] = t;
				}
				if(n == 0 && k == 0) {
					q[i] = ov;
				}
				
				p[n][k][i] = pr;
				
			}
		}
		
		fclose(fr);
	
	}
	
	double*** new_pr = (double***)calloc(NPT, sizeof(double**));
	checkHypMat(new_pr);
	
	char* ifo_filename = (char*)calloc(STRING_SIZE, sizeof(char*));
	checkString(ifo_filename);
	sprintf(ifo_filename, "parisi_dis_ave_N%d.dat", Size);
	
	FILE* file_ifo = fopen(ifo_filename, "w");
	checkFile(file_ifo, ifo_filename);
	free(ifo_filename);
	
	fprintf(file_ifo, "#overlap\tdist\terror\ttemperature\n");
	
	for(int k = 0; k < NPT; k++) {
		
		new_pr[k] = (double**)calloc(HISTBARS, sizeof(double*));
		checkMat(new_pr[k]);
		
		for(int i = 0; i < HISTBARS; i++) {
			
			new_pr[k][i] = (double*)calloc(nsamples, sizeof(double));
			checkDouble(new_pr[k][i]);
			
			for(int n = 0; n < nsamples; n++) {
				
				new_pr[k][i][n] = p[n][k][i];
			}
			
			double mean = MeanVec(new_pr[k][i], nsamples);
			double err = StandardDeviation(new_pr[k][i], mean, nsamples);
			
			fprintf(file_ifo, "%le\t%le\t%le\t%lf\n", q[i], mean, err, temp[k]);
			
		}
		
		fprintf(file_ifo, "\n");
	}
	
	fclose(file_ifo);
	
	freeMemory3(p, nsamples, NPT);
	freeMemory3(new_pr, NPT, HISTBARS);
	free(q);
	free(temp);
	


}


void disorder_average_spectrum_NOPT(int nsamples, int b) {
	
	double** freq = (double**)calloc(nsamples, sizeof(double*));
	checkMat(freq);
	
	double*** spectrum = (double***)calloc(NPT, sizeof(double**));
	checkHypMat(spectrum);
	
	double*** error = (double***)calloc(NPT, sizeof(double**));
	checkHypMat(error);
	
	double* temp = (double*)calloc(NPT, sizeof(double));
	checkDouble(temp);
	
	for(int k = 0; k < NPT; k++) {
		
		spectrum[k] = (double**)calloc(nsamples, sizeof(double*));
		checkMat(spectrum[k]);
		
		error[k] = (double**)calloc(nsamples, sizeof(double*));
		checkMat(error[k]);
		
		for(int i = 0; i < nsamples; i++) {
			
			spectrum[k][i] = (double*)calloc(Size, sizeof(double));
			checkDouble(spectrum[k][i]);
			
			error[k][i] = (double*)calloc(Size, sizeof(double));
			checkDouble(error[k][i]);
		}
	}
	
	//printf("Acquiring spectrum data...\n");
	
	for(int i = 0; i < nsamples; i++) {
		
		freq[i] = (double*)calloc(Size, sizeof(double));
		checkDouble(freq[i]);
		
		char* filename = (char*)calloc(2*STRING_SIZE, sizeof(char));
		checkString(filename);
		
		sprintf(filename, "sample%d/spectrum_N%d_block%d.dat", i+1, Size, b);
		
		FILE* fr = fopen(filename, "r");
		checkFile(fr, filename);
		free(filename);
		
		fscanf(fr, "%*[^\n]\n");
		
		for(int k = 0; k < NPT; k++) {
			
			for(int j = 0; j < Size; j++) {
				
				double f = 0.;
				double s = 0.;
				double e = 0.;
				double t = 0.;
				
				fscanf(fr, "%le\t%le\t%le\t%lf\n", &f, &s, &e, &t);
				
				if(i == 0 && j == 0) {
					temp[k] = t;
				}
				if(k == 0) {
					freq[i][j] = f;
				}
				
				spectrum[k][i][j] = s;
				error[k][i][j] = e;
			}
		}
		
		fclose(fr);
	
	}
	
	double* all_freq = (double*)calloc(nsamples*Size, sizeof(double));
	double** all_spec = (double**)calloc(NPT, sizeof(double*));
	double** all_err = (double**)calloc(NPT, sizeof(double*));
	checkMat(all_spec);
	checkMat(all_err);
	
	//printf("Ordering frequencies...\n");
	
	for(int k = 0; k < NPT; k++) {
		
		all_spec[k] = (double*)calloc(nsamples*Size, sizeof(double));
		all_err[k] = (double*)calloc(nsamples*Size, sizeof(double));
		checkDouble(all_spec[k]);
		checkDouble(all_err[k]);
		
	}
	
	
	int l = 0;
	for(int i = 0; i < nsamples; i++) {
		
		for(int j = 0; j < Size; j++) {
			
			all_freq[l] = freq[i][j];
			l++;
		}
	}

	for(int k = 0; k < NPT; k++) {
		
		l = 0;
		
		for(int i = 0; i < nsamples; i++) {
			
			for(int j = 0; j < Size; j++) {
				
				all_spec[k][l] = spectrum[k][i][j];
				all_err[k][l] = error[k][i][j];
				l++;
				
			}
		}
	}
	
	freeMemory3(spectrum, NPT, nsamples);
	freeMemory3(error, NPT, nsamples);
	
	double a = 0.;
	double bb = 0.;
	double c = 0.;
	
	for(int i = 0; i < nsamples*Size; i++) {
		
		for(int j = i+1; j < nsamples*Size; j++) {
			
			if(all_freq[i]>all_freq[j]) {
				a = all_freq[i];
				all_freq[i] = all_freq[j];
				all_freq[j] = a;
				
				for(int k = 0; k < NPT; k++) {
					
					bb = all_spec[k][i];
					all_spec[k][i] = all_spec[k][j];
					all_spec[k][j] = bb;
					
					c = all_err[k][i];
					all_err[k][i] = all_err[k][j];
					all_err[k][j] = c;
				}
				
			}
		}
	}
	
	char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(filename);
	sprintf(filename, "spectrum_all_N%d_block%d.dat", Size, b);
	
	FILE* fw = fopen(filename, "w");
	checkFile(fw, filename);
	free(filename);
	
	fprintf(fw, "#frequencies\tintensities\terror\ttemperature\n");
	
	for(int k = 0; k < NPT; k++) {
		
		for(int i = 0; i < nsamples*Size; i++) {
			
			fprintf(fw, "%le\t%le\t%le\t%lf\n", all_freq[i], all_spec[k][i]/temp[k], all_err[k][i]/temp[k], temp[k]);
		}
		
		fprintf(fw, "\n");
	}
	
	fclose(fw);
	
	//printf("Rebinning spectrum...\n");
	
	double min = 0.;
	double max = 1.;
	
	double bin_size = (max - min)/HISTBARS; 
	
	double* f = (double*)calloc(HISTBARS, sizeof(double));
	checkDouble(f);
	
	for(int i = 0; i < HISTBARS; i++) {
		f[i] = min + bin_size/2 + i* bin_size;
	}
	
	int** counter = (int**)calloc(NPT, sizeof(int*));
	//all_spec[k][i]
	for(int k = 0; k < NPT; k++) {
		
		counter[k] = (int*)calloc(HISTBARS, sizeof(int));
		checkInt(counter[k]);
		
		for(int j = 0; j < HISTBARS; j++) {
			
			for(int i = 0; i < nsamples*Size; i++) {
				
				if((all_freq[i] > min + j* bin_size) && (all_freq[i] <= min + (j+1)*bin_size)) {
					
					counter[k][j]++;
				}
			}
		}
	}
	
	char* reb_filename = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(reb_filename);
	
	sprintf(reb_filename, "spectrum_rebinned_N%d_block%d.dat", Size, b);
	FILE* reb_file = fopen(reb_filename, "w");
	checkFile(reb_file, reb_filename);
	free(reb_filename);
	
	fprintf(reb_file, "#frequencies\tintensities\terror\ttemperature\n");
	for(int k = 0; k < NPT; k++) {
		
		int i = 0;
		for(int j = 0; j < HISTBARS; j++) {
			
			double s = 0.;
			double e = 0.;
			
			for(int u = 0; u < counter[k][j]; u++) {
				
				s += all_spec[k][i]/counter[k][j];
				e += all_err[k][i]/counter[k][j];
				
				i++;
			}
			
			fprintf(reb_file, "%le\t%le\t%le\t%lf\n", f[j], s, e, temp[k]);
			
		}
		
		fprintf(reb_file, "\n");
	}
	
	fclose(reb_file);
	free(all_freq);
	freeMemory2(all_spec, NPT);
	freeMemory2(all_err, NPT);
	free(temp);
	
} 



void disorder_average_theo_ifo_NOPT(int nsamples, int b) {
	
	double* q = (double*)calloc(HISTBARS, sizeof(double));
	checkDouble(q);
	
	double* temp = (double*)calloc(NPT, sizeof(double));
	checkDouble(temp);
	
	double*** p = (double***)calloc(nsamples, sizeof(double**));
	checkHypMat(p);
	
	
	
	for(int n = 0; n < nsamples; n++) {
		
		p[n] = (double**)calloc(NPT, sizeof(double*));
		checkMat(p[n]);
		
		char* filename = (char*)calloc(2*STRING_SIZE, sizeof(char));
		checkString(filename);
		
		sprintf(filename, "sample%d/theo_ifo_histogram_block%d.dat", n+1, b);
		
		FILE* fr = fopen(filename, "r");
		checkFile(fr, filename);
		free(filename);
		
		fscanf(fr, "%*[^\n]\n");
		
		for(int k = 0; k < NPT; k++) {
			
			p[n][k] = (double*)calloc(HISTBARS, sizeof(double));
			checkDouble(p[n][k]);
			
			for(int i = 0; i < HISTBARS; i++) {
				
				double ov = 0.;
				double pr = 0.;
				double t = 0.;
				
				fscanf(fr, "%le\t%le\t%lf\n", &ov, &pr, &t);
				
				if(n == 0 && i == 0) {
					temp[k] = t;
				}
				if(n == 0 && k == 0) {
					q[i] = ov;
				}
				
				p[n][k][i] = pr;
				
			}
		}
		
		fclose(fr);
	
	}
	
	double*** new_pr = (double***)calloc(NPT, sizeof(double**));
	checkHypMat(new_pr);
	
	char* ifo_filename = (char*)calloc(STRING_SIZE, sizeof(char*));
	checkString(ifo_filename);
	sprintf(ifo_filename, "theo_ifo_dis_ave_N%d_block%d.dat", Size, b);
	
	FILE* file_ifo = fopen(ifo_filename, "w");
	checkFile(file_ifo, ifo_filename);
	free(ifo_filename);
	
	fprintf(file_ifo, "#overlap\tdist\terror\ttemperature\n");
	
	for(int k = 0; k < NPT; k++) {
		
		new_pr[k] = (double**)calloc(HISTBARS, sizeof(double*));
		checkMat(new_pr[k]);
		
		for(int i = 0; i < HISTBARS; i++) {
			
			new_pr[k][i] = (double*)calloc(nsamples, sizeof(double));
			checkDouble(new_pr[k][i]);
			
			for(int n = 0; n < nsamples; n++) {
				
				new_pr[k][i][n] = p[n][k][i];
			}
			
			double mean = MeanVec(new_pr[k][i], nsamples);
			double err = StandardDeviation(new_pr[k][i], mean, nsamples);
			
			fprintf(file_ifo, "%le\t%le\t%le\t%lf\n", q[i], mean, err, temp[k]);
			
		}
		
		fprintf(file_ifo, "\n");
	}
	
	fclose(file_ifo);
	
	freeMemory3(p, nsamples, NPT);
	freeMemory3(new_pr, NPT, HISTBARS);
	free(q);
	free(temp);
	

}




void disorder_average_exp_ifo_NOPT(int nsamples, int b) {

	double* q = (double*)calloc(HISTBARS, sizeof(double));
	checkDouble(q);
	
	double* temp = (double*)calloc(NPT, sizeof(double));
	checkDouble(temp);
	
	double*** p = (double***)calloc(nsamples, sizeof(double**));
	checkHypMat(p);
	
	
	
	for(int n = 0; n < nsamples; n++) {
		
		p[n] = (double**)calloc(NPT, sizeof(double*));
		checkMat(p[n]);
		
		char* filename = (char*)calloc(2*STRING_SIZE, sizeof(char));
		checkString(filename);
		
		sprintf(filename, "sample%d/exp_ifo_histogram_block%d.dat", n+1, b);
		
		FILE* fr = fopen(filename, "r");
		checkFile(fr, filename);
		free(filename);
		
		fscanf(fr, "%*[^\n]\n");
		
		for(int k = 0; k < NPT; k++) {
			
			p[n][k] = (double*)calloc(HISTBARS, sizeof(double));
			checkDouble(p[n][k]);
			
			for(int i = 0; i < HISTBARS; i++) {
				
				double ov = 0.;
				double pr = 0.;
				double t = 0.;
				
				fscanf(fr, "%le\t%le\t%lf\n", &ov, &pr, &t);
				
				if(n == 0 && i == 0) {
					temp[k] = t;
				}
				if(n == 0 && k == 0) {
					q[i] = ov;
				}
				
				p[n][k][i] = pr;
				
			}
		}
		
		fclose(fr);
	
	}
	
	double*** new_pr = (double***)calloc(NPT, sizeof(double**));
	checkHypMat(new_pr);
	
	char* ifo_filename = (char*)calloc(STRING_SIZE, sizeof(char*));
	checkString(ifo_filename);
	sprintf(ifo_filename, "exp_ifo_dis_ave_N%d_block%d.dat", Size, b);
	
	FILE* file_ifo = fopen(ifo_filename, "w");
	checkFile(file_ifo, ifo_filename);
	free(ifo_filename);
	
	fprintf(file_ifo, "#overlap\tdist\terror\ttemperature\n");
	
	for(int k = 0; k < NPT; k++) {
		
		new_pr[k] = (double**)calloc(HISTBARS, sizeof(double*));
		checkMat(new_pr[k]);
		
		for(int i = 0; i < HISTBARS; i++) {
			
			new_pr[k][i] = (double*)calloc(nsamples, sizeof(double));
			checkDouble(new_pr[k][i]);
			
			for(int n = 0; n < nsamples; n++) {
				
				new_pr[k][i][n] = p[n][k][i];
			}
			
			double mean = MeanVec(new_pr[k][i], nsamples);
			double err = StandardDeviation(new_pr[k][i], mean, nsamples);
			
			fprintf(file_ifo, "%le\t%le\t%le\t%lf\n", q[i], mean, err, temp[k]);
			
		}
		
		fprintf(file_ifo, "\n");
	}
	
	fclose(file_ifo);
	
	freeMemory3(p, nsamples, NPT);
	freeMemory3(new_pr, NPT, HISTBARS);
	free(q);
	free(temp);
	


}




void disorder_average_parisi_NOPT(int nsamples, int b) {

	double* q = (double*)calloc(HISTBARS, sizeof(double));
	checkDouble(q);
	
	double* temp = (double*)calloc(NPT, sizeof(double));
	checkDouble(temp);
	
	double*** p = (double***)calloc(nsamples, sizeof(double**));
	checkHypMat(p);
	
	
	
	for(int n = 0; n < nsamples; n++) {
		
		p[n] = (double**)calloc(NPT, sizeof(double*));
		checkMat(p[n]);
		
		char* filename = (char*)calloc(2*STRING_SIZE, sizeof(char));
		checkString(filename);
		
		sprintf(filename, "sample%d/parisi_histogram_block%d.dat", n+1, b);
		
		FILE* fr = fopen(filename, "r");
		checkFile(fr, filename);
		free(filename);
		
		fscanf(fr, "%*[^\n]\n");
		
		for(int k = 0; k < NPT; k++) {
			
			p[n][k] = (double*)calloc(HISTBARS, sizeof(double));
			checkDouble(p[n][k]);
			
			for(int i = 0; i < HISTBARS; i++) {
				
				double ov = 0.;
				double pr = 0.;
				double t = 0.;
				
				fscanf(fr, "%le\t%le\t%lf\n", &ov, &pr, &t);
				
				if(n == 0 && i == 0) {
					temp[k] = t;
				}
				if(n == 0 && k == 0) {
					q[i] = ov;
				}
				
				p[n][k][i] = pr;
				
			}
		}
		
		fclose(fr);
	
	}
	
	double*** new_pr = (double***)calloc(NPT, sizeof(double**));
	checkHypMat(new_pr);
	
	char* ifo_filename = (char*)calloc(STRING_SIZE, sizeof(char*));
	checkString(ifo_filename);
	sprintf(ifo_filename, "parisi_dis_ave_N%d_block%d.dat", Size, b);
	
	FILE* file_ifo = fopen(ifo_filename, "w");
	checkFile(file_ifo, ifo_filename);
	free(ifo_filename);
	
	fprintf(file_ifo, "#overlap\tdist\terror\ttemperature\n");
	
	for(int k = 0; k < NPT; k++) {
		
		new_pr[k] = (double**)calloc(HISTBARS, sizeof(double*));
		checkMat(new_pr[k]);
		
		for(int i = 0; i < HISTBARS; i++) {
			
			new_pr[k][i] = (double*)calloc(nsamples, sizeof(double));
			checkDouble(new_pr[k][i]);
			
			for(int n = 0; n < nsamples; n++) {
				
				new_pr[k][i][n] = p[n][k][i];
			}
			
			double mean = MeanVec(new_pr[k][i], nsamples);
			double err = StandardDeviation(new_pr[k][i], mean, nsamples);
			
			fprintf(file_ifo, "%le\t%le\t%le\t%lf\n", q[i], mean, err, temp[k]);
			
		}
		
		fprintf(file_ifo, "\n");
	}
	
	fclose(file_ifo);
	
	freeMemory3(p, nsamples, NPT);
	freeMemory3(new_pr, NPT, HISTBARS);
	free(q);
	free(temp);
	


}

