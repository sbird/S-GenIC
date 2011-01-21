#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "allvars.h"
#include "proto.h"


void read_parameterfile(char *fname)
{
#define FLOAT 1
#define STRING 2
#define INT 3
#define MAXTAGS 300
#define BUFSIZE 400

  FILE *fd;
  char buf[BUFSIZE], buf1[BUFSIZE], buf2[BUFSIZE], buf3[BUFSIZE];
  int i, j, nt;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int errorFlag = 0;

  /* read parameter file on all processes for simplicty */

  nt = 0;

  strcpy(tag[nt], "Omega");
  addr[nt] = &Omega;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "OmegaLambda");
  addr[nt] = &OmegaLambda;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "OmegaBaryon");
  addr[nt] = &OmegaBaryon;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "OmegaDM_2ndSpecies");
  addr[nt] = &OmegaDM_2ndSpecies;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "HubbleParam");
  addr[nt] = &HubbleParam;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "ShapeGamma");
  addr[nt] = &ShapeGamma;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Sigma8");
  addr[nt] = &Sigma8;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "PrimordialIndex");
  addr[nt] = &PrimordialIndex;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Box");
  addr[nt] = &Box;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Redshift");
  addr[nt] = &Redshift;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Nmesh");
  addr[nt] = &Nmesh;
  id[nt++] = INT;

  strcpy(tag[nt], "Nsample");
  addr[nt] = &Nsample;
  id[nt++] = INT;

  strcpy(tag[nt], "GlassFile");
  addr[nt] = GlassFile;
  id[nt++] = STRING;

  strcpy(tag[nt], "FileWithInputSpectrum");
  addr[nt] = FileWithInputSpectrum;
  id[nt++] = STRING;

  strcpy(tag[nt], "FileWithTransfer");
  addr[nt] = FileWithTransfer;
  id[nt++] = STRING;


  strcpy(tag[nt], "GlassTileFac");
  addr[nt] = &GlassTileFac;
  id[nt++] = INT;

  strcpy(tag[nt], "Seed");
  addr[nt] = &Seed;
  id[nt++] = INT;

  strcpy(tag[nt], "SphereMode");
  addr[nt] = &SphereMode;
  id[nt++] = INT;

  strcpy(tag[nt], "NumFilesWrittenInParallel");
  addr[nt] = &NumFilesWrittenInParallel;
  id[nt++] = INT;

  strcpy(tag[nt], "OutputDir");
  addr[nt] = OutputDir;
  id[nt++] = STRING;

  strcpy(tag[nt], "FileBase");
  addr[nt] = FileBase;
  id[nt++] = STRING;

  strcpy(tag[nt], "WhichSpectrum");
  addr[nt] = &WhichSpectrum;
  id[nt++] = INT;

  strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
  addr[nt] = &UnitVelocity_in_cm_per_s;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "UnitLength_in_cm");
  addr[nt] = &UnitLength_in_cm;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "UnitMass_in_g");
  addr[nt] = &UnitMass_in_g;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "InputSpectrum_UnitLength_in_cm");
  addr[nt] = &InputSpectrum_UnitLength_in_cm;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "ReNormalizeInputSpectrum");
  addr[nt] = &ReNormalizeInputSpectrum;
  id[nt++] = INT;

  strcpy(tag[nt], "WDM_On");
  addr[nt] = &WDM_On;
  id[nt++] = INT;

  strcpy(tag[nt], "WDM_Vtherm_On");
  addr[nt] = &WDM_Vtherm_On;
  id[nt++] = INT;

  strcpy(tag[nt], "WDM_PartMass_in_kev");
  addr[nt] = &WDM_PartMass_in_kev;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "NU_On");
  addr[nt] = &NU_On;
  id[nt++] = INT;
  strcpy(tag[nt], "NU_KSPACE");
  addr[nt] = &neutrinos_ks;
  id[nt++] = INT;


  strcpy(tag[nt], "NU_Vtherm_On");
  addr[nt] = &NU_Vtherm_On;
  id[nt++] = INT;

  strcpy(tag[nt], "NU_PartMass_in_ev");
  addr[nt] = &NU_PartMass_in_ev;
  id[nt++] = FLOAT;
#ifdef SPLINE
  strcpy(tag[nt], "NumKnots");
  addr[nt] = &NumKnots;
  id[nt++] = INT;

  strcpy(tag[nt], "KnotPositions");
  addr[nt] = &KnotPositions;
  id[nt++] = STRING;

  strcpy(tag[nt], "KnotValues");
  addr[nt] = &KnotValues;
  id[nt++] = STRING;
#endif
  if((fd = fopen(fname, "r")))
    {
      while(!feof(fd))
	{
	  buf[0] = 0;
	  fgets(buf, BUFSIZE, fd);
          /*Check to see if the buffer was too small: 
           * If it is large enough, should have a newline at the end*/
          for(i=0;i<BUFSIZE;i++)
          {
             if(buf[i]=='\r' || buf[i]=='\n')
                     break;
          }
          if(buf[BUFSIZE-1]=='\0' && i==BUFSIZE-1){
                  fprintf(stderr, "Error! Param line buffer not large enough.\n"
                                  "Edit read_param.c\n"
                                  "Read was:%s\n",buf);
                  FatalError(43);
          }

	  if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
	    continue;

	  if(buf1[0] == '%')
	    continue;

	  for(i = 0, j = -1; i < nt; i++)
	    if(strcmp(buf1, tag[i]) == 0)
	      {
		j = i;
		tag[i][0] = 0;
		break;
	      }

	  if(j >= 0)
	    {
	      switch (id[j])
		{
		case FLOAT:
		  *((double *) addr[j]) = atof(buf2);
		  break;
		case STRING:
		  strcpy(addr[j], buf2);
		  break;
		case INT:
		  *((int *) addr[j]) = atoi(buf2);
		  break;
		}
	    }
	  else
	    {
		fprintf(stdout, "Error in file %s:   Tag '%s' not allowed or multiple defined.\n", fname,
			buf1);
	      errorFlag = 1;
	    }
	}
      fclose(fd);

    }
  else
    {
	fprintf(stdout, "Parameter file %s not found.\n", fname);
      errorFlag = 1;
    }


  for(i = 0; i < nt; i++)
    {
      if(*tag[i])
	{
	    fprintf(stdout, "Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
	  errorFlag = 1;
	}
    }

  if(errorFlag)
      exit(0);


#undef FLOAT
#undef STRING
#undef INT
#undef MAXTAGS
}
