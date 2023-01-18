int getMasterSeed(int resume, int bash_seed) {
  int master_seed;

  char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(filename);
  FILE* seed_file;
  int temp_seed = 0;
  if(resume == 0) {
    master_seed = bash_seed;

  }else if(resume == 1) {
    sprintf(filename, "master_seed.txt");

    seed_file = fopen(filename, "r");
    checkFile(seed_file, filename);
    free(filename);
    fscanf(seed_file, "%d\n", &temp_seed);
    master_seed = temp_seed;
    fclose(seed_file);
  }else{
    printf("Resume value not valid! Exit.\n");
    exit(BAD_SIMULATION);
  }
  printf("MASTER SEED: %d\n", master_seed);
  return master_seed;
}

void printMasterSeed(int seed) {
  char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(filename);
  sprintf(filename, "master_seed.txt");
  FILE* seed_file = fopen(filename, "w");
  checkFile(seed_file, filename);
  free(filename);

  fprintf(seed_file,"%d\n", seed);
  fclose(seed_file);
}


int getGraphSeed(int resume) {
  int seed_int;
  char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(filename);
  FILE* seed_file;
  int temp_seed = 0;

  if(resume == 0) {
    seed_int = lrand48();
  }else if(resume == 1) {
    sprintf(filename, "graph_seed.txt");

    seed_file = fopen(filename, "r");
    checkFile(seed_file, filename);
    free(filename);
    fscanf(seed_file, "%d\n", &temp_seed);
    seed_int = temp_seed;
    fclose(seed_file);
  }else {
    printf("Resume value not valid! Exit.\n");
    exit(BAD_SIMULATION);
  }
  printf("INTERACTION SEED: %d\n", seed_int);
  return seed_int;


}

void printGraphSeed(int seed) {
  char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(filename);
  sprintf(filename, "graph_seed.txt");
  FILE* seed_file = fopen(filename, "w");
  checkFile(seed_file, filename);
  free(filename);

  fprintf(seed_file,"%d\n", seed);
  fclose(seed_file);
}


int* getReplicaSeed(int resume) {
  int* seed = (int*)calloc(NREPLICAS, sizeof(int));
  checkInt(seed);

  char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(filename);
  FILE* seed_file;



  if(resume == 0) {
    for(int r = 0; r < NREPLICAS; r++) {
      seed[r] = lrand48();
    }
  }else if(resume == 1) {
    sprintf(filename, "replica_seed.txt");

    seed_file = fopen(filename, "r");
    checkFile(seed_file, filename);
    free(filename);

    for(int r = 0; r < NREPLICAS; r++) {
      fscanf(seed_file, "%d\n", &seed[r]);
    }
    fclose(seed_file);
  }else {
    printf("Resume value not valid! Exit.\n");
    exit(BAD_SIMULATION);
  }

  for(int r = 0; r < NREPLICAS; r++) {
    printf("SEED REPLICA %d: %d\n", r, seed[r]);
  }
  return seed;

}

void printReplicaSeed(int* seed) {
  char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(filename);
  sprintf(filename, "replica_seed.txt");
  FILE* seed_file = fopen(filename, "w");
  checkFile(seed_file, filename);
  free(filename);
  for(int r = 0; r < NREPLICAS; r++) {
    fprintf(seed_file, "%d\n", seed[r]);
  }

  fclose(seed_file);
}


frequency genFrequencies(int resume) {
  frequency f;

  f.size = Size;
  f.gain = (double*)calloc(f.size, sizeof(double));
  checkDouble(f.gain);
  f.index = (int*)calloc(f.size, sizeof(int));
  checkInt(f.index);
  f.freq = (double*)calloc(f.size, sizeof(double));
  checkDouble(f.freq);

  char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(filename);
  FILE* freq_file;
  if(resume == 0) {

    #if FREQ_ENABLE == 1
      int nbFreq = generateDiscreteWss(N,ws,wsindices,seed2);                 // generazione frequenze con input da file
      #if _GainMax_ < 1.e-8
        for(int i=0;i<f.size;i++){
          f.gain[i] = 0.;
        }
        #else
        generateGain(N,wsindices,gain,_GainMax_ / sqrt(_Tref_));                   // generazione gain dallo stesso file di input con massimo _GainMax_
        #endif

    #elif FREQ_ENABLE == 2

      int nbFreq = f.size;
      double mySigma = f.size/4.;          // gaussian gain
      double myMean = f.size/2.-0.5;

      f.freq = GenCombLikeFrequencies(0, DW, f.size);
      for(int i=0;i<f.size;i++){
        //ws[i] = (double)i;

        f.index[i] = i ;
        f.gain[i] = _GainMax_ * ( exp(-(i-myMean)*(i-myMean)/(2.*mySigma*mySigma) )  ) / ( sqrt(twopi) * mySigma * sqrt(_Tref_) );   // gain gaussiano con sigma = mySigma
      }

    #elif FREQ_ENABLE == 3
      f.freq = GenUniformFrequencies(0, DW, f.size);
      for(int i=0;i<f.size;i++){
        //ws[i] = i;
        f.index[i] = i ;
        f.gain[i] = 0.;
      }

    #else
      int nbFreq = 1;
      for(int i=0;i<N;i++){
        f.freq[i] = 1.;
        f.index[i] = 0;
        f.gain[i] = 0.;
      }
      #endif

    #ifdef MAXNFREQ
      if(nbFreq > MAXNFREQ) {
        cout << "\n#main: ERROR: number of frequencies exceding the constant variable NBFREQ. exiting" << endl << endl;
        exit(1);
      }
    #endif

    #if I_WANT_GAIN
      FILE *fileg = fopen("final_gain.dat","w");
      for(int i=0;i<f.size;i++){
        fprintf(fileg,"%d %f \n", i, f.gain[i]);
      }
      fclose(fileg);

    #endif
  }else if(resume == 1){
    sprintf(filename, "frequencies.dat");

    freq_file = fopen(filename, "r");
    checkFile(freq_file, filename);
    free(filename);
    for(int i = 0; i < f.size; i++) {
      fscanf(freq_file, "%d\t%le\t%le\n", &f.index[i], &f.freq[i], &f.gain[i]);
    }
    fclose(freq_file);
  }else{
    printf("Resume value not valid! Exit.\n");
    exit(BAD_SIMULATION);
  }
  return f;
}


void printFrequencies(frequency f) {
  char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(filename);
  sprintf(filename, "frequencies.dat");
  FILE* freq_file = fopen(filename, "w");
  checkFile(freq_file, filename);
  free(filename);

  for(int i = 0; i < f.size; i++) {
    fprintf(freq_file, "%d\t%le\t%le\n", f.index[i], f.freq[i], f.gain[i]);
  }
  fclose(freq_file);
}

void SaveSysCriteria(Conf_type*** sys, int iter) {

  for(int r = 0; r < NREPLICAS; r++) {
    char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(filename);

    sprintf(filename, "sys_criteria_it%d_rep%d.dat",iter, r);

    FILE* conf_file = fopen(filename, "w");
    checkFile(conf_file, filename);
    //free(filename);
    fprintf(conf_file, "%d\n", iter);

    for(int k = 0; k < NPT; k++) {
      if(sys[r][k]->N) {
        fprintf(conf_file, "%d\n",1);
      }else{
        fprintf(conf_file, "%d\n",0);
      }

      if(sys[r][k]->T) {
        fprintf(conf_file, "%d\n",1);
      }else{
        fprintf(conf_file, "%d\n",0);
      }

      if(sys[r][k]->identity) {
        fprintf(conf_file, "%d\n",1);
      }else{
        fprintf(conf_file, "%d\n",0);
      }

      if(sys[r][k]->gain) {
        fprintf(conf_file, "%d\n",1);
      }else{
        fprintf(conf_file, "%d\n",0);
      }

      if(sys[r][k]->pl_ene) {
        fprintf(conf_file, "%d\n",1);
      }else{
        fprintf(conf_file, "%d\n",0);
      }


      if(sys[r][k]->pl_ene_new) {
        fprintf(conf_file, "%d\n",1);
      }else{
        fprintf(conf_file, "%d\n",0);
      }


      if(sys[r][k]->pl_de) {
        fprintf(conf_file, "%d\n",1);
      }else{
        fprintf(conf_file, "%d\n",0);
      }


      if(sys[r][k]->pl_de_block) {
        fprintf(conf_file, "%d\n",1);
      }else{
        fprintf(conf_file, "%d\n",0);
      }

      if(sys[r][k]->xs) {
        fprintf(conf_file, "%d\n",1);
      }else{
        fprintf(conf_file, "%d\n",0);
      }

      if(sys[r][k]->ys) {
        fprintf(conf_file, "%d\n",1);
      }else{
        fprintf(conf_file, "%d\n",0);
      }
    }

    fclose(conf_file);
    char* command = (char*)calloc(STRING_SIZE, sizeof(char));
    sprintf(command, "mv %s backup", filename);
    system(command);
    free(command);
    free(filename);
  }

}

void SaveConf(Conf_type*** sys, int Nplaqs, int N, int iter) {

  for(int r = 0; r < NREPLICAS; r++) {
    char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(filename);

    sprintf(filename, "conf_file_backup_it%d_rep%d.dat",iter, r);

    FILE* conf_file = fopen(filename, "w");
    checkFile(conf_file, filename);
    //free(filename);
    fprintf(conf_file, "%d\n", iter);

    for(int k = 0; k < NPT; k++) {
      if(sys[r][k]->N) {

        fprintf(conf_file, "N\n");
        fprintf(conf_file, "%d\n", sys[r][k]->N);
      }

      if(sys[r][k]->T) {

        fprintf(conf_file, "T\n");
        fprintf(conf_file, "%.20le\n", sys[r][k]->T);
      }

      if(sys[r][k]->identity) {

        fprintf(conf_file, "IDENTITY\n");
        fprintf(conf_file, "%s\n", sys[r][k]->identity);
      }

      if(sys[r][k]->gain) {

        fprintf(conf_file, "GAIN\n");
        for(int j = 0; j < N; j++) {

          fprintf(conf_file, "%.20le\n", sys[r][k]->gain[j]);
        }
      }

      if(sys[r][k]->pl_ene) {

        fprintf(conf_file, "PL_ENE\n");
        for(int n = 0; n < Nplaqs; n++) {
          fprintf(conf_file, "%.20le\n", sys[r][k]->pl_ene[n]);
        }
      }


      if(sys[r][k]->pl_ene_new) {

        fprintf(conf_file, "PL_ENE_NEW\n");
        for(int n = 0; n < Nplaqs; n++) {
          fprintf(conf_file, "%.20le\n", sys[r][k]->pl_ene_new[n]);
        }
      }


      if(sys[r][k]->pl_de) {

        fprintf(conf_file, "PL_DE\n");
        for(int n = 0; n < Nplaqs; n++) {
          fprintf(conf_file, "%.20le\n", sys[r][k]->pl_de[n]);
        }
      }


      if(sys[r][k]->pl_de_block) {

        fprintf(conf_file, "PL_DE_BLOCK\n");
        for(int n = 0; n < Nplaqs/N_THREADS_1_BLOCK; n++) {
          fprintf(conf_file, "%.20le\n", sys[r][k]->pl_de_block[n]);
        }
      }

      if(sys[r][k]->xs) {

        fprintf(conf_file, "XS\n");
        for(int j = 0; j < N; j++) {
          fprintf(conf_file, "%.20le\n", sys[r][k]->xs[j]);
        }
      }

      if(sys[r][k]->ys) {

        fprintf(conf_file, "YS\n");
        for(int j = 0; j < N; j++) {
          fprintf(conf_file, "%.20le\n", sys[r][k]->ys[j]);
        }
      }
    }

    fclose(conf_file);
    char* command = (char*)calloc(STRING_SIZE, sizeof(char));
    sprintf(command, "mv %s backup", filename);
    system(command);
    free(command);
    free(filename);
  }

}



void CheckConf(Conf_type*** sys, int Nplaqs, int N, int iter) {

  for(int r = 0; r < NREPLICAS; r++) {
    char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(filename);

    sprintf(filename, "conf_file_backup_it%d_rep%d-check.dat",iter, r);

    FILE* conf_file = fopen(filename, "w");
    checkFile(conf_file, filename);
    free(filename);
    fprintf(conf_file, "%d\n", iter);

    for(int k = 0; k < NPT; k++) {
      if(sys[r][k]->N) {

        fprintf(conf_file, "N\n");
        fprintf(conf_file, "%d\n", sys[r][k]->N);
      }

      if(sys[r][k]->T) {

        fprintf(conf_file, "T\n");
        fprintf(conf_file, "%.20le\n", sys[r][k]->T);
      }

      if(sys[r][k]->identity) {

        fprintf(conf_file, "IDENTITY\n");
        fprintf(conf_file, "%s\n", sys[r][k]->identity);
      }

      if(sys[r][k]->gain) {

        fprintf(conf_file, "GAIN\n");
        for(int j = 0; j < N; j++) {

          fprintf(conf_file, "%.20le\n", sys[r][k]->gain[j]);
        }
      }

      if(sys[r][k]->pl_ene) {

        fprintf(conf_file, "PL_ENE\n");
        for(int n = 0; n < Nplaqs; n++) {
          fprintf(conf_file, "%.20le\n", sys[r][k]->pl_ene[n]);
        }
      }


      if(sys[r][k]->pl_ene_new) {

        fprintf(conf_file, "PL_ENE_NEW\n");
        for(int n = 0; n < Nplaqs; n++) {
          fprintf(conf_file, "%.20le\n", sys[r][k]->pl_ene_new[n]);
        }
      }


      if(sys[r][k]->pl_de) {

        fprintf(conf_file, "PL_DE\n");
        for(int n = 0; n < Nplaqs; n++) {
          fprintf(conf_file, "%.20le\n", sys[r][k]->pl_de[n]);
        }
      }


      if(sys[r][k]->pl_de_block) {

        fprintf(conf_file, "PL_DE_BLOCK\n");
        for(int n = 0; n < Nplaqs/N_THREADS_1_BLOCK; n++) {
          fprintf(conf_file, "%.20le\n", sys[r][k]->pl_de_block[n]);
        }
      }

      if(sys[r][k]->xs) {

        fprintf(conf_file, "XS\n");
        for(int j = 0; j < N; j++) {
          fprintf(conf_file, "%.20le\n", sys[r][k]->xs[j]);
        }
      }

      if(sys[r][k]->ys) {

        fprintf(conf_file, "YS\n");
        for(int j = 0; j < N; j++) {
          fprintf(conf_file, "%.20le\n", sys[r][k]->ys[j]);
        }
      }
    }

    fclose(conf_file);
  }

}










void SaveTime(double t, int iter) {
  char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(filename);

  sprintf(filename, "time_file_it%d.dat",iter);

  FILE* time_file = fopen(filename, "w");
  checkFile(time_file, filename);
  free(filename);
  fprintf(time_file, "%d\n", iter);



  fprintf(time_file, "%.20le\n", t);


  fclose(time_file);



}

void SaveMcCriteria(MC_type*** mc_step, int iter) {

  for(int r = 0; r < NREPLICAS; r++) {
    char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(filename);
    sprintf(filename, "mc_criteria_it%d_rep%d.dat",iter, r);


    FILE* mc_file = fopen(filename, "w");
    checkFile(mc_file, filename);
    //free(filename);
    fprintf(mc_file, "%d\n", iter);


    for(int k = 0; k < NPT; k++) {
      if(mc_step[r][k]->T) {
        fprintf(mc_file, "%d\n",1);
      }else{
        fprintf(mc_file, "%d\n",0);
      }
      if(mc_step[r][k]->flag) {
        fprintf(mc_file, "%d\n",1);
      }else{
        fprintf(mc_file, "%d\n",0);
      }
      if(mc_step[r][k]->icoppia) {
        fprintf(mc_file, "%d\n",1);
      }else{
        fprintf(mc_file, "%d\n",0);
      }
      if(mc_step[r][k]->Ncoppie) {
        fprintf(mc_file, "%d\n",1);
      }else{
        fprintf(mc_file, "%d\n",0);
      }

      if(mc_step[r][k]->coppie) {
        fprintf(mc_file, "%d\n",1);
      }else{
        fprintf(mc_file, "%d\n",0);
      }

      if(mc_step[r][k]->rnumbers_coppie) {
        fprintf(mc_file, "%d\n",1);
      }else{
        fprintf(mc_file, "%d\n",0);
      }

      if(mc_step[r][k]->nx1) {
        fprintf(mc_file, "%d\n",1);
      }else{
        fprintf(mc_file, "%d\n",0);
      }

      if(mc_step[r][k]->nx2) {
        fprintf(mc_file, "%d\n",1);
      }else{
        fprintf(mc_file, "%d\n",0);
      }

      if(mc_step[r][k]->ny1) {
        fprintf(mc_file, "%d\n",1);
      }else{
        fprintf(mc_file, "%d\n",0);
      }

      if(mc_step[r][k]->ny2) {
        fprintf(mc_file, "%d\n",1);
      }else{
        fprintf(mc_file, "%d\n",0);
      }

      if(mc_step[r][k]->alpha_rand) {
        fprintf(mc_file, "%d\n",1);
      }else{
        fprintf(mc_file, "%d\n",0);
      }

      if(mc_step[r][k]->phi1_rand) {
        fprintf(mc_file, "%d\n",1);
      }else{
        fprintf(mc_file, "%d\n",0);
      }

      if(mc_step[r][k]->phi2_rand) {
        fprintf(mc_file, "%d\n",1);
      }else{
        fprintf(mc_file, "%d\n",0);
      }

    }


    fclose(mc_file);

    char* command = (char*)calloc(STRING_SIZE, sizeof(char));
    sprintf(command, "mv %s backup", filename);
    system(command);
    free(command);
    free(filename);
  }

}


void SaveMC(MC_type*** mc_step, int N, int Ncoppie, int iter) {

  for(int r = 0; r < NREPLICAS; r++) {
    char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(filename);
    sprintf(filename, "mc_file_backup_it%d_rep%d.dat", iter, r);

    FILE* mc_file = fopen(filename, "w");
    checkFile(mc_file, filename);

    fprintf(mc_file, "%d\n", iter);


    for(int k = 0; k < NPT; k++) {
      if(mc_step[r][k]->T) {
        fprintf(mc_file, "T\n");
        fprintf(mc_file, "%.20le\n", mc_step[r][k]->T);
      }
      if(mc_step[r][k]->flag) {
        fprintf(mc_file, "FLAG\n");
        fprintf(mc_file, "%d\n", mc_step[r][k]->flag);
      }
      if(mc_step[r][k]->icoppia) {
        fprintf(mc_file, "ICOPPIA\n");
        fprintf(mc_file, "%d\n", mc_step[r][k]->icoppia);
      }
      if(mc_step[r][k]->Ncoppie) {
        fprintf(mc_file, "NCOPPIE\n");
        fprintf(mc_file, "%d\n", mc_step[r][k]->Ncoppie);
      }

      if(mc_step[r][k]->coppie) {
        fprintf(mc_file, "NCOPPIE\n");
        for(int j = 0; j < N; j++) {
          fprintf(mc_file, "%d\n", mc_step[r][k]->coppie[j]);
        }
      }

      if(mc_step[r][k]->rnumbers_coppie) {
        fprintf(mc_file, "RNUMBERS_COPPIE\n");
        for(int j = 0; j < N; j++) {
          fprintf(mc_file, "%.20e\n", mc_step[r][k]->rnumbers_coppie[j]);
        }
      }

      if(mc_step[r][k]->nx1) {
        fprintf(mc_file, "NX1\n");
        for(int nc = 0; nc < Ncoppie; nc++) {
          fprintf(mc_file, "%.20le\n", mc_step[r][k]->nx1[nc]);
        }
      }

      if(mc_step[r][k]->nx2) {
        fprintf(mc_file, "NX2\n");
        for(int nc = 0; nc < Ncoppie; nc++) {
          fprintf(mc_file, "%.20le\n", mc_step[r][k]->nx2[nc]);
        }
      }

      if(mc_step[r][k]->ny1) {
        fprintf(mc_file, "NY1\n");


        for(int nc = 0; nc < Ncoppie; nc++) {
          fprintf(mc_file, "%.20le\n", mc_step[r][k]->ny1[nc]);
        }
      }

      if(mc_step[r][k]->ny2) {

        fprintf(mc_file, "NY2\n");

        for(int nc = 0; nc < Ncoppie; nc++) {
          fprintf(mc_file, "%.20le\n", mc_step[r][k]->ny2[nc]);
        }
      }

      if(mc_step[r][k]->alpha_rand) {

        fprintf(mc_file, "ALPHA RAND\n");

        for(int nc = 0; nc < Ncoppie; nc++) {
          fprintf(mc_file, "%.20le\n", mc_step[r][k]->alpha_rand[nc]);
        }
      }

      if(mc_step[r][k]->phi1_rand) {
        fprintf(mc_file, "PHI1 RAND\n");


        for(int nc = 0; nc < Ncoppie; nc++) {
          fprintf(mc_file, "%.20le\n", mc_step[r][k]->phi1_rand[nc]);
        }
      }

      if(mc_step[r][k]->phi2_rand) {

        fprintf(mc_file, "PHI2 RAND\n");

        for(int nc = 0; nc < Ncoppie; nc++) {
          fprintf(mc_file, "%.20le\n", mc_step[r][k]->phi2_rand[nc]);
        }
      }

    }
    fclose(mc_file);

    char* command = (char*)calloc(STRING_SIZE, sizeof(char));
    sprintf(command, "mv %s backup", filename);
    system(command);
    free(command);
    free(filename);

  }

}

void CheckMC(MC_type*** mc_step, int N, int Ncoppie, int iter) {

  for(int r = 0; r < NREPLICAS; r++) {
    char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(filename);
    sprintf(filename, "mc_file_backup_it%d_rep%d-check.dat", iter, r);

    FILE* mc_file = fopen(filename, "w");
    checkFile(mc_file, filename);
    free(filename);
    fprintf(mc_file, "%d\n", iter);


    for(int k = 0; k < NPT; k++) {
      if(mc_step[r][k]->T) {
        fprintf(mc_file, "T\n");
        fprintf(mc_file, "%.20le\n", mc_step[r][k]->T);
      }
      if(mc_step[r][k]->flag) {
        fprintf(mc_file, "FLAG\n");
        fprintf(mc_file, "%d\n", mc_step[r][k]->flag);
      }
      if(mc_step[r][k]->icoppia) {
        fprintf(mc_file, "ICOPPIA\n");
        fprintf(mc_file, "%d\n", mc_step[r][k]->icoppia);
      }
      if(mc_step[r][k]->Ncoppie) {
        fprintf(mc_file, "NCOPPIE\n");
        fprintf(mc_file, "%d\n", mc_step[r][k]->Ncoppie);
      }

      if(mc_step[r][k]->coppie) {
        fprintf(mc_file, "NCOPPIE\n");
        for(int j = 0; j < N; j++) {
          fprintf(mc_file, "%d\n", mc_step[r][k]->coppie[j]);
        }
      }

      if(mc_step[r][k]->rnumbers_coppie) {
        fprintf(mc_file, "RNUMBERS_COPPIE\n");
        for(int j = 0; j < N; j++) {
          fprintf(mc_file, "%.20e\n", mc_step[r][k]->rnumbers_coppie[j]);
        }
      }

      if(mc_step[r][k]->nx1) {
        fprintf(mc_file, "NX1\n");
        for(int nc = 0; nc < Ncoppie; nc++) {
          fprintf(mc_file, "%.20le\n", mc_step[r][k]->nx1[nc]);
        }
      }

      if(mc_step[r][k]->nx2) {
        fprintf(mc_file, "NX2\n");
        for(int nc = 0; nc < Ncoppie; nc++) {
          fprintf(mc_file, "%.20le\n", mc_step[r][k]->nx2[nc]);
        }
      }

      if(mc_step[r][k]->ny1) {
        fprintf(mc_file, "NY1\n");


        for(int nc = 0; nc < Ncoppie; nc++) {
          fprintf(mc_file, "%.20le\n", mc_step[r][k]->ny1[nc]);
        }
      }

      if(mc_step[r][k]->ny2) {

        fprintf(mc_file, "NY2\n");

        for(int nc = 0; nc < Ncoppie; nc++) {
          fprintf(mc_file, "%.20le\n", mc_step[r][k]->ny2[nc]);
        }
      }

      if(mc_step[r][k]->alpha_rand) {

        fprintf(mc_file, "ALPHA RAND\n");

        for(int nc = 0; nc < Ncoppie; nc++) {
          fprintf(mc_file, "%.20le\n", mc_step[r][k]->alpha_rand[nc]);
        }
      }

      if(mc_step[r][k]->phi1_rand) {
        fprintf(mc_file, "PHI1 RAND\n");


        for(int nc = 0; nc < Ncoppie; nc++) {
          fprintf(mc_file, "%.20le\n", mc_step[r][k]->phi1_rand[nc]);
        }
      }

      if(mc_step[r][k]->phi2_rand) {

        fprintf(mc_file, "PHI2 RAND\n");

        for(int nc = 0; nc < Ncoppie; nc++) {
          fprintf(mc_file, "%.20le\n", mc_step[r][k]->phi2_rand[nc]);
        }
      }

    }


    fclose(mc_file);
  }

}



void SaveClockCriteria(Clock_type** clock, int iter) {


  for(int r = 0; r < NREPLICAS; r++) {
    char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(filename);
    sprintf(filename, "clock_criteria_it%d_rep%d.dat", iter, r);

    FILE* clock_file = fopen(filename, "w");
    checkFile(clock_file, filename);
    //free(filename);
    fprintf(clock_file, "%d\n", iter);



      if(clock[r]->prof_time) {
        fprintf(clock_file, "%d\n", 1);
      }else{
        fprintf(clock_file, "%d\n", 0);
      }

      if(clock[r]->nrg) {
        fprintf(clock_file, "%d\n", 1);
      }else{
        fprintf(clock_file, "%d\n", 0);
      }

      if(clock[r]->n_attemp_exchange) {
        fprintf(clock_file, "%d\n", 1);
      }else{
        fprintf(clock_file, "%d\n", 0);
      }
      if(clock[r]->acc_rate_exchange) {
        fprintf(clock_file, "%d\n", 1);
      }else{
        fprintf(clock_file, "%d\n", 0);
      }

      if(clock[r]->acc_rate) {
        fprintf(clock_file, "%d\n", 1);
      }else{
        fprintf(clock_file, "%d\n", 0);
      }

      if(clock[r]->n_attemp) {
        fprintf(clock_file, "%d\n", 1);
      }else{
        fprintf(clock_file, "%d\n", 0);
      }



      fclose(clock_file);

      char* command = (char*)calloc(STRING_SIZE, sizeof(char));
      sprintf(command, "mv %s backup", filename);
      system(command);
      free(command);
      free(filename);
  }

}

void SaveClock(Clock_type** clock, int iter) {



  for(int r = 0; r < NREPLICAS; r++) {
    char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(filename);
    sprintf(filename, "clock_file_backup_it%d_rep%d.dat", iter, r);

    FILE* clock_file = fopen(filename, "w");
    checkFile(clock_file, filename);

    fprintf(clock_file, "%d\n", iter);

      if(clock[r]->prof_time) {
        fprintf(clock_file, "PROF_TIME\n");
        for(int n = 0; n < NTIMES_PROFILING; n++) {
          fprintf(clock_file, "%.20le\n", clock[r]->prof_time[n]);
        }
      }

      if(clock[r]->nrg) {
        fprintf(clock_file, "NRG\n");
        for(int l = 0; l < NPT; l++) {
          fprintf(clock_file, "%.20le\n", clock[r]->nrg[l]);
        }
      }

      if(clock[r]->n_attemp_exchange) {
        fprintf(clock_file, "N_ATTEMP_EXCHANGE\n");
        for(int n = 0; n < NTJUMPS; n++) {
          fprintf(clock_file, "%d\n", clock[r]->n_attemp_exchange[n]);
        }
      }
      if(clock[r]->acc_rate_exchange) {
        fprintf(clock_file, "ACC_RATE_EXCHANGE\n");
        for(int n = 0; n < NTJUMPS; n++) {
          fprintf(clock_file, "%d\n", clock[r]->acc_rate_exchange[n]);
        }
      }

      if(clock[r]->acc_rate) {
        fprintf(clock_file, "ACC_RATE\n");
        for(int n = 0; n < NPT; n++) {
          fprintf(clock_file, "%d\n", clock[r]->acc_rate[n]);
        }
      }

      if(clock[r]->n_attemp) {
        fprintf(clock_file, "N_ATTEMP\n");
        for(int n = 0; n < NPT; n++) {
          fprintf(clock_file, "%d\n", clock[r]->n_attemp[n]);
        }
      }

      fclose(clock_file);

      char* command = (char*)calloc(STRING_SIZE, sizeof(char));
      sprintf(command, "mv %s backup", filename);
      system(command);
      free(command);
      free(filename);

    }


}




void CheckClock(Clock_type** clock, int iter) {



  for(int r = 0; r < NREPLICAS; r++) {
    char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(filename);
    sprintf(filename, "clock_file_backup_it%d_rep%d-check.dat", iter, r);

    FILE* clock_file = fopen(filename, "w");
    checkFile(clock_file, filename);
    free(filename);
    fprintf(clock_file, "%d\n", iter);

      if(clock[r]->prof_time) {
        fprintf(clock_file, "PROF_TIME\n");
        for(int n = 0; n < NTIMES_PROFILING; n++) {
          fprintf(clock_file, "%.20le\n", clock[r]->prof_time[n]);
        }
      }

      if(clock[r]->nrg) {
        fprintf(clock_file, "NRG\n");
        for(int l = 0; l < NPT; l++) {
          fprintf(clock_file, "%.20le\n", clock[r]->nrg[l]);
        }
      }

      if(clock[r]->n_attemp_exchange) {
        fprintf(clock_file, "N_ATTEMP_EXCHANGE\n");
        for(int n = 0; n < NTJUMPS; n++) {
          fprintf(clock_file, "%d\n", clock[r]->n_attemp_exchange[n]);
        }
      }
      if(clock[r]->acc_rate_exchange) {
        fprintf(clock_file, "ACC_RATE_EXCHANGE\n");
        for(int n = 0; n < NTJUMPS; n++) {
          fprintf(clock_file, "%d\n", clock[r]->acc_rate_exchange[n]);
        }
      }

      if(clock[r]->acc_rate) {
        fprintf(clock_file, "ACC_RATE\n");
        for(int n = 0; n < NPT; n++) {
          fprintf(clock_file, "%d\n", clock[r]->acc_rate[n]);
        }
      }

      if(clock[r]->n_attemp) {
        fprintf(clock_file, "N_ATTEMP\n");
        for(int n = 0; n < NPT; n++) {
          fprintf(clock_file, "%d\n", clock[r]->n_attemp[n]);
        }
      }

      fclose(clock_file);


    }


}
// DA RIFARE DACCAPO
last_iteration lastConf(int iter) {

  last_iteration choise;
  char** loadconf_filename = (char**)calloc(NREPLICAS, sizeof(char*));
  char** loadmc_filename = (char**)calloc(NREPLICAS, sizeof(char*));
  char** loadclock_filename = (char**)calloc(NREPLICAS, sizeof(char*));
  char** criteriaconf_filename = (char**)calloc(NREPLICAS, sizeof(char*));
  char** criteriamc_filename = (char**)calloc(NREPLICAS, sizeof(char*));
  char** criteriaclock_filename = (char**)calloc(NREPLICAS, sizeof(char*));
  choise.loadconf = (char**)calloc(NREPLICAS, sizeof(char*));
  choise.loadmc = (char**)calloc(NREPLICAS, sizeof(char*));
  choise.loadclock = (char**)calloc(NREPLICAS, sizeof(char*));
  choise.criteriaconf = (char**)calloc(NREPLICAS, sizeof(char*));
  choise.criteriamc = (char**)calloc(NREPLICAS, sizeof(char*));
  choise.criteriaclock = (char**)calloc(NREPLICAS, sizeof(char*));
  for(int r = 0; r < NREPLICAS; r++) {
    loadconf_filename[r] = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(loadconf_filename[r]);
	  loadmc_filename[r] = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(loadmc_filename[r]);
	  loadclock_filename[r] = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(loadclock_filename[r]);
	  criteriaconf_filename[r] = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(criteriaconf_filename[r]);
	  criteriamc_filename[r] = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(criteriamc_filename[r]);
	  criteriaclock_filename[r] = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(criteriaclock_filename[r]);


    sprintf(loadconf_filename[r], "backup/conf_file_backup_it%d_rep%d.dat",iter, r);
    sprintf(loadmc_filename[r],"backup/mc_file_backup_it%d_rep%d.dat", iter, r);
    sprintf(loadclock_filename[r], "backup/clock_file_backup_it%d_rep%d.dat", iter, r);
    sprintf(criteriaconf_filename[r], "backup/sys_criteria_it%d_rep%d.dat",iter, r);
    sprintf(criteriamc_filename[r], "backup/mc_criteria_it%d_rep%d.dat",iter, r);
    sprintf(criteriaclock_filename[r], "backup/clock_criteria_it%d_rep%d.dat", iter, r);

    FILE* f = fopen(loadconf_filename[r], "r");
    checkFile(f, loadconf_filename[r]);
    int temp_iter = 0;
    fscanf(f, "%d\n", &temp_iter);
    printf("%s -> %d\n", loadconf_filename[r], NSTEP*temp_iter);
    fclose(f);

    f = fopen(loadmc_filename[r], "r");
    checkFile(f, loadmc_filename[r]);
    fscanf(f, "%d\n", &temp_iter);
    printf("%s -> %d\n", loadmc_filename[r], NSTEP*temp_iter);
    fclose(f);

    f = fopen(loadclock_filename[r], "r");
    checkFile(f, loadclock_filename[r]);
    fscanf(f, "%d\n", &temp_iter);
    printf("%s -> %d\n", loadclock_filename[r], NSTEP*temp_iter);
    fclose(f);

    f = fopen(criteriaconf_filename[r], "r");
    checkFile(f, criteriaconf_filename[r]);
    fscanf(f, "%d\n", &temp_iter);
    printf("%s -> %d\n", criteriaconf_filename[r], NSTEP*temp_iter);
    fclose(f);

    f = fopen(criteriamc_filename[r], "r");
    checkFile(f, criteriamc_filename[r]);
    fscanf(f, "%d\n", &temp_iter);
    printf("%s -> %d\n", criteriamc_filename[r], NSTEP*temp_iter);
    fclose(f);

    f = fopen(criteriaclock_filename[r], "r");
    checkFile(f, criteriaclock_filename[r]);
    fscanf(f, "%d\n", &temp_iter);
    printf("%s -> %d\n", criteriaclock_filename[r], NSTEP*temp_iter);
    fclose(f);
    /*
    f = fopen(loadtime_filename[r], "r");
    checkFile(f, loadtime_filename[r]);
    fscanf(f, "%d\n", &temp_iter);
    printf("%s -> %d\n", loadtime_filename[r], NSTEP*temp_iter);
    fclose(f);
    */
    choise.loadconf[r] = loadconf_filename[r];
    choise.loadmc[r] = loadmc_filename[r];
    choise.loadclock[r] = loadclock_filename[r];
    choise.criteriaconf[r] = criteriaconf_filename[r];
    choise.criteriamc[r] = criteriamc_filename[r];
    choise.criteriaclock[r] = criteriaclock_filename[r];
    //choise[r].loadtime = loadtime_filename[r];
  }

  free(loadconf_filename);
  free(loadmc_filename);
  free(loadclock_filename);
  free(criteriaconf_filename);
  free(criteriamc_filename);
  free(criteriaclock_filename);

  return choise;

}

int*** readCriteria(int size, char** filename) {

  int*** criteria = (int***)calloc(NREPLICAS, sizeof(int**));

  for(int r = 0; r < NREPLICAS; r++) {

    FILE* fread = fopen(filename[r], "r");
    checkFile(fread, filename[r]);


    int temp_iter;

    fscanf(fread, "%d\n", &temp_iter);


    criteria[r] = (int**)calloc(NPT, sizeof(int*));
    for(int k = 0; k < NPT; k++) {
      criteria[r][k] = (int*)calloc(size, sizeof(int));
      for(int n = 0; n < size; n++) {
        fscanf(fread, "%d\n", &criteria[r][k][n]);

      }

    }

  }

  return criteria;
}

double LoadTime(int iter) {

  char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(filename);
  sprintf(filename, "time_file_it%d.dat", iter);
  FILE* fread=fopen(filename, "r");
  checkFile(fread, filename);

  double t =0.;

  int temp_iter;

  fscanf(fread, "%d\n", &temp_iter);

    fscanf(fread, "%le\n", &t);
    printf("TIME SWEEP: %le\n", t);


  return t;
}

Conf_type*** LoadConf(char** filename, int Nplaqs, int N, int*** criteria) {



  Conf_type*** sys = (Conf_type***)calloc(NREPLICAS, sizeof(Conf_type**));

  for(int r = 0; r < NREPLICAS; r++) {

    printf("filename: %s\n", filename[r]);
    FILE* fread = fopen(filename[r], "r");
    checkFile(fread, filename[r]);
    int temp_iter = 0;
    fscanf(fread, "%d\n", &temp_iter);

    sys[r] = (Conf_type**)calloc(NPT, sizeof(Conf_type*));
    for(int k = 0; k < NPT; k++) {
      char* temp_str = (char*)calloc(STRING_SIZE, sizeof(char));
      sys[r][k] = (Conf_type*)calloc(1, sizeof(Conf_type));

      sys[r][k]->identity = (char*)calloc(NCHAR_IDENTITY, sizeof(char));
      sys[r][k]->gain = (double*)calloc(N, sizeof(double));
      sys[r][k]->pl_ene = (double*)calloc(Nplaqs, sizeof(double));
      sys[r][k]->pl_ene_new = (double*)calloc(Nplaqs, sizeof(double));
      sys[r][k]->pl_de = (double*)calloc(Nplaqs, sizeof(double));
      sys[r][k]->pl_de_block = (double*)calloc(int(Nplaqs/N_THREADS_1_BLOCK), sizeof(double));
      sys[r][k]->xs = (spin_t*)calloc(N, sizeof(spin_t));
      sys[r][k]->ys = (spin_t*)calloc(N, sizeof(spin_t));


      if(criteria[r][k][0] == 1) {
        fscanf(fread, "%s\n", temp_str);
        fscanf(fread, "%d\n", &sys[r][k]->N);
        //printf("%d\n", sys[r][k][0].N);
      }
      if(criteria[r][k][1] == 1) {
        fscanf(fread, "%s\n", temp_str);

        fscanf(fread, "%le\n", &sys[r][k]->T);
        if(PRINTLOADING==1) printf("T:\n");
        if(PRINTLOADING==1) printf("%le\n", sys[r][k]->T);

      }
      if(criteria[r][k][2] == 1) {
        checkString(sys[r][k][0].identity);
        fscanf(fread, "%s\n", temp_str);

        fscanf(fread, "%s\n", sys[r][k]->identity);
        if(PRINTLOADING==1) printf("IDENTITY:\n");
        if(PRINTLOADING==1) printf("%s\n", sys[r][k]->identity);

      }

      if(criteria[r][k][3] == 1) {
        fscanf(fread, "%s\n", temp_str);
        if(PRINTLOADING==1) printf("GAIN:\n");
        for(int j = 0; j < N; j++) {
          fscanf(fread, "%le\n",&sys[r][k]->gain[j]);
          if(PRINTLOADING==1) printf("%le\n", sys[r][k]->gain[j]);
        }
      }


      if(criteria[r][k][4] == 1) {
        fscanf(fread, "%s\n", temp_str);
        if(PRINTLOADING==1) printf("PL_ENE:\n");
        for(int n = 0; n < Nplaqs; n++) {
          fscanf(fread, "%le\n", &sys[r][k]->pl_ene[n]);
          if(PRINTLOADING==1) printf("%le\n", sys[r][k]->pl_ene[n]);

        }
      }

      if(criteria[r][k][5] == 1) {
        fscanf(fread, "%s\n", temp_str);
        if(PRINTLOADING==1) printf("PL_ENE_EW:\n");

        for(int n = 0; n < Nplaqs; n++) {
          fscanf(fread, "%le\n", &sys[r][k]->pl_ene_new[n]);
          if(PRINTLOADING==1) printf("%le\n", sys[r][k]->pl_ene_new[n]);
        }
      }
      if(criteria[r][k][6] == 1) {
        fscanf(fread, "%s\n", temp_str);
        if(PRINTLOADING==1) printf("PL_DE\n");

        for(int n = 0; n < Nplaqs; n++) {
          fscanf(fread, "%le\n", &sys[r][k]->pl_de[n]);
          if(PRINTLOADING==1) printf("%le\n", sys[r][k]->pl_de[n]);
        }
      }
      if(criteria[r][k][7] == 1) {
        fscanf(fread, "%s\n", temp_str);
        if(PRINTLOADING==1) printf("PL_DE_BLOCK\n");

        for(int n = 0; n < Nplaqs/N_THREADS_1_BLOCK; n++) {
          fscanf(fread, "%le\n", &sys[r][k]->pl_de_block[n]);
          if(PRINTLOADING==1) printf("%le\n", sys[r][k]->pl_de_block[n]);
        }
      }


      if(criteria[r][k][8] == 1) {
        fscanf(fread, "%s\n", temp_str);
        if(PRINTLOADING==1) printf("XS\n");

        for(int n = 0; n < N; n++) {
          fscanf(fread, "%le\n", &sys[r][k]->xs[n]);
          if(PRINTLOADING==1) printf("%le\n", sys[r][k]->xs[n]);
        }
      }
      if(criteria[r][k][9] == 1) {
        fscanf(fread, "%s\n", temp_str);
        if(PRINTLOADING==1) printf("YS:\n");

        for(int n = 0; n < N; n++) {
          fscanf(fread, "%le\n", &sys[r][k]->ys[n]);
          if(PRINTLOADING==1) printf("%le\n", sys[r][k]->ys[n]);
        }
      }
      free(temp_str);

    }
    fclose(fread);

  }

  return sys;


}

MC_type*** LoadMCS(char** filename, int Ncoppie, int N, int*** criteria) {


  MC_type*** mc = (MC_type***)calloc(NREPLICAS, sizeof(MC_type**));
  //SIZE STRUCT 13
  for(int r = 0; r < NREPLICAS; r++) {

    printf("filename: %s\n", filename[r]);
    FILE* fread = fopen(filename[r], "r");
    checkFile(fread, filename[r]);
    int temp_iter = 0;
    fscanf(fread, "%d\n", &temp_iter);

    mc[r] = (MC_type**)calloc(NPT, sizeof(MC_type*));
    for(int k = 0; k < NPT; k++) {
      char* temp_str = (char*)calloc(STRING_SIZE, sizeof(char));

      mc[r][k] = (MC_type*)calloc(1, sizeof(MC_type));
      mc[r][k]->coppie = (int*)calloc(N, sizeof(int));
      mc[r][k]->rnumbers_coppie = (float*)calloc(N, sizeof(float));
      mc[r][k]->nx1 = (spin_t*)calloc(Ncoppie, sizeof(spin_t));
      mc[r][k]->nx2 = (spin_t*)calloc(Ncoppie, sizeof(spin_t));
      mc[r][k]->ny1 = (spin_t*)calloc(Ncoppie, sizeof(spin_t));
      mc[r][k]->ny2 = (spin_t*)calloc(Ncoppie, sizeof(spin_t));
      mc[r][k]->alpha_rand = (double*)calloc(Ncoppie, sizeof(double));
      mc[r][k]->phi1_rand = (double*)calloc(Ncoppie, sizeof(double));
      mc[r][k]->phi2_rand = (double*)calloc(Ncoppie, sizeof(double));

      if(criteria[r][k][0] == 1) {
        fscanf(fread, "%s\n", temp_str);
        if(PRINTLOADING==1) printf("T:\n");

        fscanf(fread, "%le\n", &mc[r][k]->T);
        if(PRINTLOADING==1) printf("%le\n",mc[r][k]->T);
      }
      if(criteria[r][k][1] == 1) {
        fscanf(fread, "%s\n", temp_str);
        if(PRINTLOADING==1) printf("FLAG:\n");
        fscanf(fread, "%d\n", &mc[r][k]->flag);
        if(PRINTLOADING==1) printf("%d\n", mc[r][k]->flag);
      }
      if(criteria[r][k][2] == 1) {
        fscanf(fread, "%s\n", temp_str);
        if(PRINTLOADING==1) printf("ICOPPIA:\n");
        fscanf(fread, "%d\n", &mc[r][k]->icoppia);
        if(PRINTLOADING==1) printf("%d\n", mc[r][k]->icoppia);
      }
      if(criteria[r][k][3] == 1) {
        fscanf(fread, "%s\n", temp_str);
        if(PRINTLOADING==1) printf("NCOPPIE:\n");
        fscanf(fread, "%d\n", &mc[r][k]->Ncoppie);
        if(PRINTLOADING==1) printf("%d\n", mc[r][k]->Ncoppie);
      }
      if(criteria[r][k][4] == 1) {
        fscanf(fread, "%s\n", temp_str);
        if(PRINTLOADING==1) printf("COPPIE:\n");

        for(int i = 0; i < N; i++) {
          fscanf(fread, "%d ", &mc[r][k]->coppie[i]);
          if(PRINTLOADING==1) printf("%d\n", mc[r][k]->coppie[i]);
        }
      }
      if(criteria[r][k][5] == 1) {
        fscanf(fread, "%s\n", temp_str);
        if(PRINTLOADING==1) printf("RNUMBERS_COPPIE:\n");

        for(int i = 0; i < N; i++) {
          fscanf(fread, "%e\n", &mc[r][k]->rnumbers_coppie[i]);
          if(PRINTLOADING==1) printf("%e\n", mc[r][k]->rnumbers_coppie[i]);
        }
      }
      if(criteria[r][k][6] == 1) {
        fscanf(fread, "%s\n", temp_str);
        if(PRINTLOADING==1) printf("NX1\n");

        for(int i = 0; i < Ncoppie; i++) {
          fscanf(fread, "%le\n", &mc[r][k]->nx1[i]);
          if(PRINTLOADING==1) printf("%le\n", mc[r][k]->nx1[i]);
        }
      }
      if(criteria[r][k][7] == 1) {
        fscanf(fread, "%s\n", temp_str);
        if(PRINTLOADING==1) printf("NX2\n");

        for(int i = 0; i < Ncoppie; i++) {
          fscanf(fread, "%le\n", &mc[r][k]->nx2[i]);
          if(PRINTLOADING==1) printf("%le\n", mc[r][k]->nx2[i]);
        }
      }
      if(criteria[r][k][8] == 1) {
        fscanf(fread, "%s\n", temp_str);
        if(PRINTLOADING==1) printf("NY1\n");

        for(int i = 0; i < Ncoppie; i++) {
          fscanf(fread, "%le\n", &mc[r][k]->ny1[i]);
          if(PRINTLOADING==1) printf("%le\n", mc[r][k]->ny1[i]);
        }
      }

      if(criteria[r][k][9] == 1) {
        fscanf(fread, "%s\n", temp_str);
        if(PRINTLOADING==1) printf("NY2\n");

        for(int i = 0; i < Ncoppie; i++) {
          fscanf(fread, "%le\n", &mc[r][k]->ny2[i]);
          if(PRINTLOADING==1) printf("%le\n", mc[r][k]->ny2[i]);
        }
      }

      if(criteria[r][k][10] == 1) {
        fscanf(fread, "%s\n", temp_str);
        if(PRINTLOADING==1) printf("ALPHA_RAND:\n");

        for(int i = 0; i < Ncoppie; i++) {
          fscanf(fread, "%le\n", &mc[r][k]->alpha_rand[i]);
          if(PRINTLOADING==1) printf("%le\n", mc[r][k]->alpha_rand[i]);
        }
      }
      if(criteria[r][k][11] == 1) {
        fscanf(fread, "%s\n", temp_str);
        if(PRINTLOADING==1) printf("PHI1_RAND:\n");

        for(int i = 0; i < Ncoppie; i++) {
          fscanf(fread, "%le\n", &mc[r][k]->phi1_rand[i]);
          if(PRINTLOADING==1) printf("%le\n", mc[r][k]->phi1_rand[i]);
        }
      }
      if(criteria[r][k][12] == 1) {
        fscanf(fread, "%s\n", temp_str);
        if(PRINTLOADING==1) printf("PHI2_RAND:\n");

        for(int i = 0; i < Ncoppie; i++) {
          fscanf(fread, "%le\n", &mc[r][k]->phi2_rand[i]);
          if(PRINTLOADING==1) printf("%le\n", mc[r][k]->phi2_rand[i]);
        }

      }
      free(temp_str);
    }
    fclose(fread);

  }
  return mc;
}

int** clockCritera(int size, char** filename) {



  int** criteria = (int**)calloc(NREPLICAS, sizeof(int*));
  for(int r = 0; r < NREPLICAS; r++) {
    FILE* fread = fopen(filename[r], "r");
    checkFile(fread, filename[r]);
    int temp_iter;

    fscanf(fread, "%d\n", &temp_iter);
    criteria[r] = (int*)calloc(size, sizeof(int));
    for(int n = 0; n < size; n++) {
      fscanf(fread, "%d\n", &criteria[r][n]);
    }
    fclose(fread);

  }


  return criteria;
}

Clock_type** LoadClock(char** filename, int** criteria) {


  Clock_type** clock=(Clock_type**)calloc(NREPLICAS, sizeof(Clock_type*));
  for(int r = 0; r < NREPLICAS; r++) {

    char* temp_str = (char*)calloc(STRING_SIZE, sizeof(char));

    printf("filename: %s\n", filename[r]);
    FILE* fread = fopen(filename[r], "r");
    checkFile(fread, filename[r]);
    int temp_iter = 0;
    fscanf(fread, "%d\n", &temp_iter);

    clock[r] = (Clock_type*)calloc(1, sizeof(Clock_type));
    clock[r]->prof_time = (double*)calloc(NTIMES_PROFILING, sizeof(double));
    clock[r]->nrg = (double*)calloc(NPT, sizeof(double));
    clock[r]->acc_rate = (int*)calloc(NPT, sizeof(double));
    clock[r]->n_attemp = (int*)calloc(NPT, sizeof(int));


    if(criteria[r][0] == 1) {
      fscanf(fread, "%s\n", temp_str);
      if(PRINTLOADING==1) printf("PROF_TIME:\n");

      for(int i = 0; i < NTIMES_PROFILING; i++) {
        fscanf(fread, "%le\n", &clock[r]->prof_time[i]);
        if(PRINTLOADING==1) printf("%le\n", clock[r]->prof_time[i]);

      }

    }
    if(criteria[r][1] == 1) {
      fscanf(fread, "%s\n", temp_str);
      if(PRINTLOADING==1) printf("NRG:\n");

      for(int i = 0; i < NPT; i++) {
        fscanf(fread, "%le\n", &clock[r]->nrg[i]);
        if(PRINTLOADING==1) printf("%le\n", clock[r]->nrg[i]);

      }
    }
    if(criteria[r][2] == 1) {
      fscanf(fread, "%s\n", temp_str);
      if(PRINTLOADING==1) printf("N_ATTEMP_EXCHANGE:\n");
      for(int i = 0; i < NTJUMPS; i++) {
        clock[r]->n_attemp_exchange[i] = 0;
        fscanf(fread, "%d\n", &clock[r]->n_attemp_exchange[i]);
        if(PRINTLOADING==1) printf("%d\n", clock[r]->n_attemp_exchange[i]);
      }
    }
    if(criteria[r][3] == 1) {
      fscanf(fread, "%s\n", temp_str);
      if(PRINTLOADING==1) printf("ACC_RATE_EXCHANGE:\n");
      for(int i = 0; i < NTJUMPS; i++) {
        clock[r]->acc_rate_exchange[i] = 0;
        fscanf(fread, "%d\n",&clock[r]->acc_rate_exchange[i]);
        if(PRINTLOADING==1) printf("%d\n", clock[r]->acc_rate_exchange[i]);
      }
    }

    if(criteria[r][4] == 1) {
      fscanf(fread, "%s\n", temp_str);
      if(PRINTLOADING==1) printf("ACC_RATE:\n");

      for(int i = 0; i < NPT; i++) {
        fscanf(fread, "%d\n", &clock[r]->acc_rate[i]);
        if(PRINTLOADING==1) printf("%d\n", clock[r]->acc_rate[i]);
      }
    }
    if(criteria[r][5] == 1) {
      fscanf(fread, "%s\n", temp_str);
      if(PRINTLOADING==1) printf("N_ATTEMP:\n");

      for(int i = 0; i < NPT; i++) {
        fscanf(fread, "%d\n", &clock[r]->n_attemp[i]);
        if(PRINTLOADING==1) printf("%d\n", clock[r]->n_attemp[i]);
      }
    }
    free(temp_str);
    fclose(fread);

  }

  return clock;

}




void Resume_counters(Clock_type * clock, Clock_type * d_clock){


  double* point_prof_time;
  double* point_nrg;
  int* point_acc_rate;
  int* point_n_attemp;
  int* point_n_attemp_exchange;
  int* point_acc_rate_exchange;

  double* temp_prof_time = (double*)calloc(NTIMES_PROFILING, sizeof(double));
  for(int i = 0; i < NTIMES_PROFILING; i++) {
    temp_prof_time[i] = clock->prof_time[i];
  }
  double* temp_nrg = (double*)calloc(NPT, sizeof(double));
  for(int i = 0; i < NPT; i++) {
    temp_nrg[i] = clock->nrg[i];
  }
  int* temp_acc_rate = (int*)calloc(NPT, sizeof(double));
  for(int i = 0; i < NPT; i++) {
    temp_acc_rate[i] = clock->acc_rate[i];
  }
  int* temp_n_attemp = (int*)calloc(NPT, sizeof(double));
  for(int i = 0; i < NPT; i++) {
    temp_n_attemp[i] = clock->n_attemp[i];
  }
  int* temp_n_attemp_exchange = (int*)calloc(NTJUMPS, sizeof(int));
  for(int i = 0; i < NTJUMPS; i++) {
    temp_n_attemp_exchange[i] = clock->n_attemp_exchange[i];
  }
  int* temp_acc_rate_exchange = (int*)calloc(NTJUMPS, sizeof(int));
  for(int i = 0; i < NTJUMPS; i++) {
    temp_acc_rate_exchange[i] = clock->acc_rate_exchange[i];
  }

  cudaMalloc((double **) &point_prof_time,NTIMES_PROFILING*sizeof(double));
  cudaMalloc((double **) &point_nrg,NPT*sizeof(double));
  cudaMalloc((int **) &point_acc_rate,NPT*sizeof(int));
  cudaMalloc((int **) &point_n_attemp,NPT*sizeof(int));
  cudaMalloc((int**)&point_n_attemp_exchange, NTJUMPS*sizeof(int));
  cudaMalloc((int**)&point_acc_rate_exchange, NTJUMPS*sizeof(int));

  cudaMemcpy(&(d_clock->prof_time),&(point_prof_time),sizeof(double *),cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_clock->nrg),&(point_nrg),sizeof(double *),cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_clock->acc_rate),&(point_acc_rate),sizeof(int *),cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_clock->n_attemp),&(point_n_attemp),sizeof(int *),cudaMemcpyHostToDevice);
  //cudaMemcpy(&(d_clock->n_attemp_exchange), &(point_n_attemp_exchange, sizeof(int*)), cudaMemcpyHostToDevice);
  //cudaMemcpy(&(d_clock->acc_rate_exchange), &(point_acc_rate_exchange, sizeof(int*)), cudaMemcpyHostToDevice);

  cudaMemcpy(point_prof_time,temp_prof_time,NTIMES_PROFILING*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy(point_nrg,temp_nrg,NPT*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy(point_acc_rate,temp_acc_rate,NPT*sizeof(int),cudaMemcpyHostToDevice);
  cudaMemcpy(point_n_attemp,temp_n_attemp,NPT*sizeof(int),cudaMemcpyHostToDevice);
  //cudaMemcpy(point_n_attemp_exchange, NTJUMPS*sizeof(int), cudaMemcpyHostToDevice);
  //cudaMemcpy(point_acc_rate_exchange, NTJUMPS*sizeof(int), cudaMemcpyHostToDevice);
  //cudaDeviceSynchronize();


  free(temp_prof_time);
  free(temp_nrg);
  free(temp_acc_rate);
  free(temp_n_attemp);
  free(temp_n_attemp_exchange);
  free(temp_acc_rate_exchange);


}



void Resume_System_device(int ireplica, int seed, int N, int Nplaqs, double temperature, Conf_type * sys, Conf_type * d_sys){

  int N_temp = N;
  double T_temp = temperature;
  char* identity_temp = (char*)calloc(NCHAR_IDENTITY, sizeof(char));
  sprintf(identity_temp, "%s", sys->identity);
  double* gain_temp = (double*)calloc(N, sizeof(double));
  for(int i = 0; i < N; i++) {
    gain_temp[i] = sys->gain[i];
  }
  double* pl_ene_temp = (double*)calloc(Nplaqs, sizeof(double));
  for(int i = 0; i < Nplaqs; i++) {
    pl_ene_temp[i] = sys->pl_ene[i];
  }

  double* pl_ene_new_temp = (double*)calloc(Nplaqs, sizeof(double));
  for(int i = 0; i < Nplaqs; i++) {
    pl_ene_new_temp[i] = sys->pl_ene_new[i];
  }
  double* pl_de_temp = (double*)calloc(Nplaqs, sizeof(double));
  for(int i = 0; i < Nplaqs; i++) {
    pl_de_temp[i] = sys->pl_de[i];
  }
  double* pl_de_block_temp = (double*)calloc(Nplaqs/N_THREADS_1_BLOCK, sizeof(double));
  for(int i = 0; i < Nplaqs/N_THREADS_1_BLOCK; i++) {
    pl_de_block_temp[i] = sys->pl_de_block[i];
  }
  spin_t* xs_temp = (spin_t*)calloc(N, sizeof(spin_t));
  for(int i = 0; i < N; i++) {
    xs_temp[i] = sys->xs[i];
  }
  spin_t* ys_temp = (spin_t*)calloc(N, sizeof(spin_t));
  for(int i = 0; i < N; i++) {
    ys_temp[i] = sys->ys[i];
  }

  char* point_identity;
  double* point_gain;
  double* point_pl_ene;
  double* point_pl_ene_new;
  double* point_pl_de;
  double* point_pl_de_block;
  spin_t* point_xs;
  spin_t* point_ys;
//ALLOCO NELLA GPU
  cudaMalloc((char**)&point_identity, NCHAR_IDENTITY*sizeof(char));
  cudaMalloc((double**)&point_gain, N*sizeof(double));
  cudaMalloc((double**)&point_pl_ene, Nplaqs*sizeof(double));
  cudaMalloc((double**)&point_pl_ene_new, Nplaqs*sizeof(double));
  cudaMalloc((double**)&point_pl_de, Nplaqs*sizeof(double));
  cudaMalloc((double**)&point_pl_de_block, (Nplaqs/N_THREADS_1_BLOCK)*sizeof(double));
  cudaMalloc((spin_t**)&point_xs, N*sizeof(spin_t));
  cudaMalloc((spin_t**)&point_ys, N*sizeof(spin_t));

  //COPIO GLI INDIRIZZI

  cudaMemcpy(&(d_sys->identity), &(point_identity), sizeof(char*), cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_sys->gain), &(point_gain), sizeof(double*), cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_sys->pl_ene), &(point_pl_ene), sizeof(double*), cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_sys->pl_ene_new), &(point_pl_ene_new), sizeof(double*), cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_sys->pl_de), &(point_pl_de), sizeof(double*), cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_sys->pl_de_block), &(point_pl_de_block), sizeof(double*), cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_sys->xs), &(point_xs), sizeof(spin_t*), cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_sys->ys), &(point_ys), sizeof(spin_t*), cudaMemcpyHostToDevice);

  //COPIO I VALORI

  cudaMemcpy(point_identity, identity_temp, NCHAR_IDENTITY*sizeof(char), cudaMemcpyHostToDevice);
  cudaMemcpy(point_gain, gain_temp, N*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(point_pl_ene, pl_ene_temp, Nplaqs*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(point_pl_ene_new, pl_ene_new_temp, Nplaqs*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(point_pl_de, pl_de_temp, Nplaqs*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(point_pl_de_block, pl_de_block_temp, (Nplaqs/N_THREADS_1_BLOCK)*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(point_xs, xs_temp, N*sizeof(spin_t), cudaMemcpyHostToDevice);
  cudaMemcpy(point_ys, ys_temp, N*sizeof(spin_t), cudaMemcpyHostToDevice);

  cudaMemcpy(&(d_sys->T),&(T_temp),sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_sys->N),&(N_temp),sizeof(int),cudaMemcpyHostToDevice);

  //cudaDeviceSynchronize();

  free(identity_temp);
  free(gain_temp);
  free(pl_ene_temp);
  free(pl_ene_new_temp);
  free(pl_de_temp);
  free(pl_de_block_temp);
  free(xs_temp);
  free(ys_temp);



}


void Resume_MC_step_variables(int ireplica, double T, int N, int Ncoppie, MC_type * mc_step, MC_type * d_mc_step){

  double temp_T = T;
  int temp_flag = mc_step->flag;
  int temp_icoppia = mc_step->icoppia;
  int temp_Ncoppie = mc_step->Ncoppie;
  int* temp_coppie = (int*)calloc(N, sizeof(int));
  for(int i = 0; i < N; i++) {
    temp_coppie[i] = mc_step->coppie[i];
  }
  float* temp_rnumbers_coppie = (float*)calloc(N, sizeof(float));
  for(int i = 0; i < N; i++) {
    temp_rnumbers_coppie[i] = mc_step->rnumbers_coppie[i];
  }
  spin_t* temp_nx1 = (spin_t*)calloc(Ncoppie, sizeof(spin_t));
  spin_t* temp_nx2 = (spin_t*)calloc(Ncoppie, sizeof(spin_t));
  spin_t* temp_ny1 = (spin_t*)calloc(Ncoppie, sizeof(spin_t));
  spin_t* temp_ny2 = (spin_t*)calloc(Ncoppie, sizeof(spin_t));

  for(int i = 0; i < Ncoppie; i++) {
    temp_nx1[i] = mc_step->nx1[i];
    temp_nx2[i] = mc_step->nx2[i];
    temp_ny1[i] = mc_step->ny1[i];
    temp_ny2[i] = mc_step->ny2[i];

  }

  double* temp_alpha_rand = (double*)calloc(Ncoppie, sizeof(double));
  double* temp_phi1_rand = (double*)calloc(Ncoppie, sizeof(double));
  double* temp_phi2_rand = (double*)calloc(Ncoppie, sizeof(double));

  for(int i = 0; i < Ncoppie; i++) {
    temp_alpha_rand[i] = mc_step->alpha_rand[i];
    temp_phi1_rand[i] = mc_step->phi1_rand[i];
    temp_phi2_rand[i] = mc_step->phi2_rand[i];
  }

  int* point_coppie;
  float* point_rnumbers_coppie;
  spin_t* point_nx1;
  spin_t* point_nx2;
  spin_t* point_ny1;
  spin_t* point_ny2;
  double* point_alpha_rand;
  double* point_phi1_rand;
  double* point_phi2_rand;


  //cudaMalloc((char**)&point_identity, NCHAR_IDENTITY*sizeof(char));

  cudaMalloc((int**)&point_coppie, N*sizeof(int));
  cudaMalloc((float**)&point_rnumbers_coppie, N*sizeof(float));
  cudaMalloc((spin_t**)&point_nx1, Ncoppie*sizeof(spin_t));
  cudaMalloc((spin_t**)&point_nx2, Ncoppie*sizeof(spin_t));
  cudaMalloc((spin_t**)&point_ny1, Ncoppie*sizeof(spin_t));
  cudaMalloc((spin_t**)&point_ny2, Ncoppie*sizeof(spin_t));

  cudaMalloc((double**)&point_alpha_rand, Ncoppie*sizeof(double));
  cudaMalloc((double**)&point_phi1_rand, Ncoppie*sizeof(double));
  cudaMalloc((double**)&point_phi2_rand, Ncoppie*sizeof(double));

  //cudaMemcpy(&(d_sys->idendity), &(point_identity), sizeof(char*), cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_mc_step->coppie), &(point_coppie), sizeof(int*), cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_mc_step->rnumbers_coppie), &(point_rnumbers_coppie), sizeof(float*), cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_mc_step->nx1), &(point_nx1), sizeof(spin_t*), cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_mc_step->nx2), &(point_nx2), sizeof(spin_t*), cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_mc_step->ny1), &(point_ny1), sizeof(spin_t*), cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_mc_step->ny2), &(point_ny2), sizeof(spin_t*), cudaMemcpyHostToDevice);

  cudaMemcpy(&(d_mc_step->alpha_rand), &(point_alpha_rand), sizeof(double*), cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_mc_step->phi1_rand), &(point_phi1_rand), sizeof(double*), cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_mc_step->phi2_rand), &(point_phi2_rand), sizeof(double*), cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_mc_step->T),&(temp_T),sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_mc_step->flag),&(temp_flag),sizeof(int),cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_mc_step->icoppia),&(temp_icoppia),sizeof(int),cudaMemcpyHostToDevice);
  cudaMemcpy(&(d_mc_step->Ncoppie),&(temp_Ncoppie),sizeof(int),cudaMemcpyHostToDevice);


//  cudaMemcpy(point_identity, identity_temp, NCHAR_IDENTITY*sizeof(char), cudaMemcpyHostToDevice);

  cudaMemcpy(point_coppie, temp_coppie, N*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(point_rnumbers_coppie, temp_rnumbers_coppie, N*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(point_nx1, temp_nx1, Ncoppie*sizeof(spin_t), cudaMemcpyHostToDevice);
  cudaMemcpy(point_nx2, temp_nx2, Ncoppie*sizeof(spin_t), cudaMemcpyHostToDevice);
  cudaMemcpy(point_ny1, temp_ny1, Ncoppie*sizeof(spin_t), cudaMemcpyHostToDevice);
  cudaMemcpy(point_ny2, temp_ny2, Ncoppie*sizeof(spin_t), cudaMemcpyHostToDevice);
  cudaMemcpy(point_alpha_rand, temp_alpha_rand, Ncoppie*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(point_phi1_rand, temp_phi1_rand, Ncoppie*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(point_phi2_rand, temp_phi2_rand, Ncoppie*sizeof(double), cudaMemcpyHostToDevice);
  //cudaDeviceSynchronize();

  free(temp_coppie);
  free(temp_rnumbers_coppie);
  free(temp_nx1);
  free(temp_nx2);
  free(temp_ny1);
  free(temp_ny2);
  free(temp_alpha_rand);
  free(temp_phi1_rand);
  free(temp_phi2_rand);


}

void resume_energy_plaquettes(long long int Nplaqs, int nreplica, Int_type * d_inter, Conf_type ** sys, Conf_type ** d_sys, Clock_type * clock, Clock_type * d_clock){

  double * point_pl_ene;

  for(int irep=0; irep<NPT; irep++){

    total_energy_parallel(d_sys[irep],d_inter); // initialize energy of each plaquette

    cudaMemcpy(&(point_pl_ene),&(d_sys[irep]->pl_ene),sizeof(double *),cudaMemcpyDeviceToHost);

    cudaMemcpy(sys[irep]->pl_ene,point_pl_ene,Nplaqs*sizeof(double),cudaMemcpyDeviceToHost);

    clock->nrg[irep]=0;

    for(int iplaq=0; iplaq<Nplaqs; iplaq++) clock->nrg[irep]+=sys[irep]->pl_ene[iplaq];

    printf("Energia totale: %g\n", clock->nrg[irep]/Size);

  }

  cudaMemcpy(&(point_pl_ene),&(d_clock->nrg),sizeof(double *),cudaMemcpyDeviceToHost);
  cudaMemcpy(point_pl_ene,clock->nrg,NPT*sizeof(double),cudaMemcpyDeviceToHost);
  cudaDeviceSynchronize();
}





void resume_replica(int seed, int N, int Nplaqs, int Ncoppie, double * temp, Conf_type ** sys, Conf_type ** d_sys, MC_type ** mc_step, MC_type ** d_mc_step){

  // ALLOCA LA STRUTTURA DATI "conf"

  for(int i=0;i<NPT;i++){
    cudaMalloc((Conf_type **) &(d_sys[i]),sizeof(Conf_type));
    Resume_System_device(i,seed,N,Nplaqs,temp[i],sys[i],d_sys[i]);
  }

  for(int i=0; i<NPT; i++){
    cudaMalloc((MC_type **) &(d_mc_step[i]),sizeof(MC_type));
    Resume_MC_step_variables(i,temp[i],N,Ncoppie,mc_step[i],d_mc_step[i]);
  }

}




void resumeParallelTempering(int r, int nPT, int last) {//THIS FUNCTION READ PARALLEL_TEMPERING%d.dat

	char* fileName = (char*)calloc(STRING_SIZE, sizeof(char));
	checkString(fileName);
  sprintf(fileName, "parallel_tempering%d.dat", r);

  char* fileName_restored = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(fileName_restored);
  sprintf(fileName_restored, "parallel_tempering%d_restored.dat", r);
	FILE* fr = fopen(fileName, "r");
	checkFile(fr, fileName);
	//free(fileName);
  FILE* fw = fopen(fileName_restored, "w");
  checkFile(fw, fileName_restored);
	//TEMP VARIABLES
	int mcs;

	int nCol = nPT - 1;


	for(int j = 0; j < last; j++) { // CICLO SULLE RIGHE
		fscanf(fr, " %d \t", &mcs);
		//printf("\n%d\n\n", mcs);
    fprintf(fw, " %d \t", mcs);

		for(int i = 0; i < nCol; i++) { //CICLO SULLE COLONNE
      double T = 0.;
      double a = 0.;
      double e = 0.;
			fscanf(fr, " %lf %lg %le \t", &T, &a, &e);
      fprintf(fw," %g %g %8.8e \t", T, a, e);

		}//FINE CICLO SULLE COLONNE
    //fscanf(fr, "\n");
		fprintf(fw, "\n");
	} //FINE CICLO SULLE RIGHE

	fclose(fr);
  fclose(fw);
  char* command = (char*)calloc(10*STRING_SIZE, sizeof(char));
  sprintf(command, "mv %s %s", fileName_restored, fileName);

  system(command);
  free(command);
  free(fileName);
  free(fileName_restored);

}



void RemoveBackup(int iter) {

  for(int r = 0; r < NREPLICAS; r++) {
    char* command = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(command);
    sprintf(command, "rm sys_criteria_it%d_rep%d.dat", iter, r);
    system(command);
    command = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(command);
    sprintf(command, "rm conf_file_backup_it%d_rep%d.dat", iter, r);
    system(command);
    command = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(command);
    sprintf(command, "rm mc_criteria_it%d_rep%d.dat", iter, r);
    system(command);
    command = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(command);
    sprintf(command, "rm mc_file_backup_it%d_rep%d.dat", iter, r);
    system(command);
    command = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(command);
    sprintf(command, "rm clock_criteria_it%d_rep%d.dat", iter, r);
    system(command);
    command = (char*)calloc(STRING_SIZE, sizeof(char));
    checkString(command);
    sprintf(command, "rm clock_file_backup_it%d_rep%d.dat", iter, r);
    system(command);
    free(command);
  }
  char* command_time = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(command_time);
  sprintf(command_time, "rm time_file_it%d.dat", iter);
  system(command_time);
  free(command_time);
}

void SaveIter(int iter) {
  char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(filename);
  sprintf(filename, "last_iter.dat");
  FILE* fw = fopen(filename, "a");
  checkFile(fw, filename);
  fprintf(fw, "%d\n", iter);
  fclose(fw);
  free(filename);
}

int getIter() {

  char* filename = (char*)calloc(STRING_SIZE, sizeof(char));
  checkString(filename);
  sprintf(filename, "last_iter.dat");
  FILE* fr = fopen(filename, "r");
  checkFile(fr, filename);
  int temp_iter = 0;
  int counter = 0;
  printf("GETTING ITERATIONS.\n");
  while (fscanf(fr, "%d\n", &temp_iter) != EOF) {
      printf("%d\n", temp_iter);
        counter++;
    }
/*
  if(counter < 3) {
    printf("Not enough iterations have been saved.\n");
    exit(BAD_SIMULATION);
  }
  */
  fclose(fr);
  printf("SAVED %d ITERATIONS.\n", counter);
  fr = fopen(filename, "r");
  int* iter = (int*)calloc(counter, sizeof(int));
  for(int i = 0; i < counter; i++) {
    fscanf(fr, "%d\n", &iter[i]);
  }
  fclose(fr);
  return (iter[counter-2]);
}



void BackupSimulation(Conf_type*** sys, MC_type*** mc_step, Clock_type** clock_r, int ind_iter, int N, int Nplaqs, int Ncoppie, double time_MC_sweep) {

  if((ind_iter%NITER_PRINT_CONF==1)&&(ind_iter>0)) {
    SaveSysCriteria(sys,ind_iter);
    SaveMcCriteria(mc_step, ind_iter);
    SaveClockCriteria(clock_r, ind_iter);
    SaveConf(sys, Nplaqs, N,ind_iter);
    SaveMC(mc_step,N,Ncoppie,ind_iter);
    SaveClock(clock_r, ind_iter);
    SaveTime(time_MC_sweep, ind_iter);

    SaveIter(ind_iter);
    //RemoveBackup(ind_iter-SAVEBUFFER);
  }
}
