
/**
   Read input parameters
 */


// for mkdir
#include <sys/stat.h>
#include <sys/types.h>

void trim_whitespace(char* s) {
  const char* d = s;
  do {
    while (*d == ' ')
      ++d;
  } while (*s++ = *d++);
  *s = '\0';
}


void str2array(char *tmps2, double *array){
  char* p;
  int n = 0;
  p = strtok(tmps2,"[,]");
  while (p != NULL){
    array[n] = atof(p);
    p = strtok(NULL, ",");
    n += 1;
  }
}

void read_params(char* path2file)
{
  FILE * fp;
  if ((fp = fopen(path2file, "rt"))) {
    char tempbuff[300];
    while(fgets(tempbuff,300,fp)) {
      trim_whitespace(tempbuff);
      char* tmps1 = strtok(tempbuff, "=");
      char* tmps2 = strtok(NULL, "=");
      // basilisk constants
      if      (strcmp(tmps1,"N")    ==0) { N     = atoi(tmps2); }
//      else if (strcmp(tmps1,"nl")   ==0) { nl    = atoi(tmps2); }
      else if (strcmp(tmps1,"DT")   ==0) { DT    = atof(tmps2); }
      else if (strcmp(tmps1,"CFL")  ==0) { CFL   = atof(tmps2); }
      else if (strcmp(tmps1,"TOLERANCE")==0) { TOLERANCE= atof(tmps2); }
      // qg specific constants
      else if (strcmp(tmps1,"delta_nl")==0){
        delta_nl = atof(tmps2);
        L0 = pow(tau0/(pow(delta_nl*beta,2.)), 1./3.);
        fprintf(stdout, "nonlinear thickness is %g.\n", pow(tau0/(pow(L0, 3.)*pow(beta, 2.)),0.5));
      }
      else if (strcmp(tmps1,"Re")==0){
        Re = atof(tmps2);
        nu = pow(tau0*pow(L0,3.)/pow(Re,4.), 1./2.);
        fprintf(stdout, "Sublayer Re is %g.\n", pow(tau0*pow(L0,3.)/pow(nu,2.),1./4.));
      }
      else if (strcmp(tmps1,"f0")   ==0) { f0    = atof(tmps2); }
      else if (strcmp(tmps1,"beta") ==0) { beta  = atof(tmps2); }
      else if (strcmp(tmps1,"hEkb") ==0) { hEkb  = atof(tmps2); }
      else if (strcmp(tmps1,"tau0") ==0) { tau0  = atof(tmps2); }
      else if (strcmp(tmps1,"sbc")  ==0) { sbc   = atof(tmps2); }
      else if (strcmp(tmps1,"tend") ==0) { tend  = atof(tmps2); }
      else if (strcmp(tmps1,"dtout")==0) { dtout = atof(tmps2); }
      else if (strcmp(tmps1,"dh")   ==0) { str2array(tmps2, dh);}
      else if (strcmp(tmps1,"iend") ==0) { iend  = atof(tmps2); }

//      printf("%s => %s\n", tmps1, tmps2);
    }
    fclose(fp);
  } else {
    fprintf(stdout, "file %s not found\n", path2file);
    exit(0);
  }

  /**
     Viscosity CFL = 0.5
   */
  if (beta  != 0) DT = min(3.1415/(2.*beta*L0),0.5*min(DT,sq(L0/N)/nu/4.));
  else 0.5*min(DT,sq(L0/N)/nu/4.);

  fprintf(stdout, "Config: N = %d, L0 = %g\n", N,  L0);
//  fprintf(stdout, "Config: N = %d, nl = %d, L0 = %g\n", N, nl, L0);
}

/**
   Create output directory and copy input parameter file for backup
*/
void create_outdir()
{
  if (pid() == 0) {
    for (int i=1; i<10000; i++) {
      sprintf(dpath, "outdir_%04d/", i);
      if (mkdir(dpath, 0777) == 0) {
        fprintf(stdout,"Writing output in %s\n",dpath);
        break;
      }
    }
  }
@if _MPI
  MPI_Bcast(&dpath, 80, MPI_CHAR, 0, MPI_COMM_WORLD);
@endif
}

void backup_config()
{
  fprintf(stdout, "Backup config\n");
  char ch;
  char name[120];
  sprintf (name,"%sparams.in", dpath);
  FILE * source = fopen("params.in", "r");
  FILE * target = fopen(name, "w");
  while ((ch = fgetc(source)) != EOF)
    fputc(ch, target);
  fclose(source);
  fclose(target);
}
