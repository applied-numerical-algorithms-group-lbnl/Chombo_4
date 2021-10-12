#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <ctype.h>

#include "crunchflow.h"

#ifdef __cplusplus
extern"C" {
#endif

  char *trimwhitespace(char *str);
  int maxval(int m, int n, int a[][n]);

  int initialize(const int ncompchombo, 
		 const int nspecchombo, 
		 const int nkinchombo,
		 double inflowchombo[ncompchombo],
		 double inflowchombosp10[ncompchombo],
		 double initchombo[ncompchombo],
		 double initchombosp[ncompchombo + nspecchombo],
		 double initchombosp10[ncompchombo + nspecchombo],
		 double stoichchombo[ncompchombo][nkinchombo][1], 
		 const int nx, 
		 const int ny, 
		 const int nz)
  {
    double wtaq[ncompchombo + nspecchombo];
    char *namcx[nspecchombo];
    char *namin[nkinchombo];

    nreactmin = (int (*))calloc(nkinchombo, sizeof(int));
    volmol = (double (*))calloc(nkinchombo, sizeof(double));

    ndepend = (int (*)[MAX_PATH])calloc(nkinchombo * MAX_PATH, sizeof(int));
    AffinityDepend1 = (double (*)[MAX_PATH])calloc(nkinchombo * MAX_PATH, sizeof(double));
    rate0 = (double (*)[MAX_PATH])calloc(nkinchombo * MAX_PATH, sizeof(double));

    idepend = (int (*)[MAX_PATH][MAX_DEPEND])calloc(nkinchombo * MAX_PATH * MAX_DEPEND, sizeof(int));
    itot_min = (int (*)[MAX_PATH][MAX_DEPEND])calloc(nkinchombo * MAX_PATH * MAX_DEPEND, sizeof(int));
    depend = (double (*)[MAX_PATH][MAX_DEPEND])calloc(nkinchombo * MAX_PATH * MAX_DEPEND, sizeof(double));

    // ncomp+nspec,nx,ny,nz
    //
    // double sp[1][32][32][15]; // allocate(sp(ncomp+nspec,nx,ny,nz),stat=ierr); sp = 0.0d0
    double *sp = (double *)calloc((ncompchombo + nspecchombo) * nx * ny * nz, sizeof(double));

    // double sp10[1][32][32][15]; // allocate(sp10(ncomp+nspec,nx,ny,nz),stat=ierr); sp10 = 0.0d0
    double *sp10 = (double *)calloc((ncompchombo + nspecchombo) * nx * ny * nz, sizeof(double));

    // double gam[1][32][32][15]; // allocate(gam(ncomp+nspec,nx,ny,nz),stat=ierr); gam = 0.0d0
    double *gam = (double *)calloc((ncompchombo + nspecchombo) * nx * ny * nz, sizeof(double));

    // double s[1][32][32][6]; // allocate(s(neqn,nx,ny,nz),stat=ierr); s = 0.0d0
    double *s = (double *)calloc(ncompchombo * nx * ny * nz, sizeof(double));

    // double sn[1][32][32][6]; // allocate(sn(neqn,nx,ny,nz),stat=ierr); s = 0.0d0
    double *sn = (double *)calloc(ncompchombo * nx * ny * nz, sizeof(double));

    // double sion[1][32][32]; // allocate(sion(nx,ny,nz),stat=ierr); sion = 0.0d0
    double *sion = (double *)calloc(nx * ny * nz, sizeof(double));

    // double t[1][32][32]; // allocate(t(nx,ny,nz),stat=ierr); t = 25.0d0
    double *t = (double *)calloc(nx * ny * nz, sizeof(double));

    double *mumin = (double *)calloc(ncompchombo * nkinchombo * MAX_PATH, sizeof(double));

    double *muaq = (double *)calloc(ncompchombo * nspecchombo, sizeof(double));

    double (* acmp);
    acmp = (double (*))calloc(ncompchombo * nspecchombo, sizeof(double));

    double (* chg);
    chg = (double (*))calloc(ncompchombo * nspecchombo, sizeof(double));

    double (* alnk);
    alnk = (double (*))calloc(nkinchombo * MAX_PATH, sizeof(double));

    double (* eqhom);
    eqhom = (double (*))calloc(nspecchombo, sizeof(double));

    char *(* ulab);
    ulab = (char * (*))calloc(ncompchombo + nspecchombo, sizeof(char *));

    // [NBASIS][nspec + nkin * max_path]
    double *as1 = (double *)calloc(NBASIS * (nspecchombo + nkinchombo * MAX_PATH), sizeof(double));

    // dppt = (double (*)[ny][nx][nkinchombo])calloc(nz * ny * nx * nkinchombo, sizeof(double));
    // u_rate = (double (*)[ny][nx][nkinchombo])calloc(nz * ny * nx * nkinchombo, sizeof(double));
    // area = (double (*)[ny][nx][nkinchombo])calloc(nz * ny * nx * nkinchombo, sizeof(double));
    // volfx = (double (*)[ny][nx][nkinchombo])calloc(nz * ny * nx * nkinchombo, sizeof(double));

    // local variables
    //
    int neqn;
    int iflgint = 0;
    int INT[NTEMP];
    int ntt;
    int nncnt;
    double temptmp[NTEMP];
    int nbasis = NBASIS;
    double temp = 25.0; // C
    int ndim2, ndim3;
    double initgamma[ncompchombo + nspecchombo];
    int np = 0;
    int nd = 0;
    int ndim1;

    read_kinetics = false; // global


    for(size_t i = 0; i < NTEMP; ++i)
      for(size_t j = 0; j < NBASIS; ++j)
      {
        vec[i][j] = 0.0; // clear vec
        vecgam[i][j] = 0.0; // clear vecgam
      }
    
    RunIsothermal = false; 
	
    //     ________________________ DEBYE HUCKEL _____________________
    //
    for(int i = 0; i < TPOINTS; i++)
    {
      tempc[i] = -333.0;
      a_dh[i] = -333.0;
      b_dh[i] = -333.0;
      b_dot[i] = -333.0;
    }

    //       open(unit=112,file='input_test.txt',status='unknown')
    char* filename = "input_test.txt";
    if(DEBUG)
      printf("parsing %s...\n",filename);
    int status = parseInput(filename);
    if(DEBUG)
      printf("  done\n");

    // temp = 25.0; //     default to 25 deg C, non-isothermal when ntemp == 8

    if(ncomp != ncompchombo)
    {
      printf("ncomp in crunchflow is not equal to ncomp in Chombo\n");
      exit(1);
    }
    else if(nspec != nspecchombo)
    {
      printf("nspec in crunchflow is not equal to nspec in Chombo\n");
      exit(1);
    }
    else if(nkin != nkinchombo)
    {
      printf("nkin in crunchflow is not equal to nkin in Chombo\n");
      exit(1);
    }

    ndim1 = nspec + nkin * MAX_PATH;

    for(int i = 0; i < TPOINTS; i++)
      if ( tempc[i] == -333.0 )
        break;
      else
        ntemp++;

    if (ntemp == 0)
    {
      //  backward compatibility - no temp data provided, assume 25 deg C
      ntemp = 1;
      tempc[0] = 25.0;
    }
    else if (ntemp == 1 && tempc[0] != 25.0)
    {
      printf("error: if only one temperature point, temperature must be = 25 deg C\n");
      exit;
    }
    else if (ntemp == 8 )
      ;
    else
    {
      printf("error: number of temp points in input_test.txt must be 1 or 8\n");
      exit;
    }

    for(int i = 0; i < ntemp; i++)
      if(a_dh[i] == -333.0 || b_dh[i]  == -333.0 || b_dot[i] == -333.0)
      {
        printf("error: a_dh, b_dh or b_dot have not been initialized for one or more temp points\n");
        exit;
      }

    if(ntemp == 1)
    {
      adh[1]  = a_dh[1];
      bdh[1]  = b_dh[1];
      bdot[1] = b_dot[1];
    }
    else // ntemp == 8
    {
      for(int i = 0; i < ntemp; i++)
      {
        vec[i][1-1] = log(tempc[i] + TK);
        vec[i][2-1] = 1.0;
        vec[i][3-1] = tempc[i] + TK;
        vec[i][4-1] = 1.0/(tempc[i] + TK);
        vec[i][5-1] = 1.0/((tempc[i] + TK)*(tempc[i] + TK));

        vecgam[i][1-1] = 1.0;
        vecgam[i][2-1] = tempc[i];
        vecgam[i][3-1] = pow(tempc[i],2);
        vecgam[i][4-1] = pow(tempc[i],3);
        vecgam[i][5-1] = pow(tempc[i],4);
      }

      //       fit gamma
      status = fitgamma_cpp(adh);
      for(int i = 0; i < NBASIS; i++)
        adhcoeff[i] = bvec[i];
      status = fitgamma_cpp(bdh);
      for(int i = 0; i < NBASIS; i++)
        bdhcoeff[i] = bvec[i];
      status = fitgamma_cpp(bdot);
      for(int i = 0; i < NBASIS; i++)
        bdtcoeff[i] = bvec[i];
    }

    neqn = ncomp;

    //     _______________________ ALLOCATE _________________________
    //


    //     _______________________ COMPONENTS _______________________
    //
    //     read all components   
    if (ncomp <= 0)
    {
      printf("no components\n");
      exit;
    }
    else
    {
      for(size_t i = 0; i < ncomp; i++)
      {
        ulab[i] = (char *)malloc((strlen(components[i].name) + 1) * sizeof(char));
        ulab[i] = strcpy(ulab[i], components[i].name);
        chg[i]  = components[i].charge;
        // printf("initialize_cpp: component %d, charge %f, name %s\n",i,chg[i],components[i].name);
        wtaq[i] = components[i].molecular_weight;
        acmp[i] = components[i].a_zero;
      }
    }

    //     ______________________ COMPLEXATION _____________________
    if (nspec < 0)
    {
      printf("negative aqueous reactions");
      exit;
    }
    else if (nspec == 0)
      printf("no aqueous reactions\n");
    else
    {
      //     loop over reactions
      for(int i = 0; i < nspec; i++) // i is column index in muaq
      {
        // initialize row j in muaq[ncomp][nspec]
        for(int j = 0; j < ncomp; j++)
          muaq[j*nspec+i] = 0.0;

        // example stoichiometry: -2.0 'H+'  1.0 'CO2(aq)'  1.0 'Ca++'
        char *tokenPtr = strtok(complexes[i].stoichiometry, " "); // begin tokenizing stoichiometry
        int tokenNum = 0;
        double mu = 0.0;
        char *p;
        char *q;

        // continue tokenizing line until tokenPtr becomes NULL
        while (tokenPtr != NULL) 
        {
          if(tokenNum % 2 == 0)
          {
            // printf("stoichiometry: mu %s\n", tokenPtr);
            sscanf(tokenPtr,"%lf",&mu);
          }
          else
          {
            // printf("stoichiometry: component %s\n", tokenPtr);
            if((p = strpbrk(tokenPtr,"'")))
            {
              q = strrchr(p,'\'');
              if(q)
                *q = '\0';
              int j;
              for(j = 0; j < ncomp; j++)
                if(!strcmp(p+1, ulab[j]))
                {
                  muaq[j*nspec+i] = mu;
                  break;
                }
              if(j >= ncomp)
              {
                fprintf(stderr,"could not find component %s for reaction %s\n",p,complexes[i].name);
                exit(1);
              }
            }
            mu = 0.0;
          }
          tokenNum++;
          tokenPtr = strtok(NULL, " "); // get next token
        }

        chg[ncomp+i] = complexes[i].charge;
        wtaq[ncomp+i] = complexes[i].molecular_weight;
        acmp[ncomp+i] = complexes[i].a_zero;
        namcx[i] = (char *)malloc((strlen(complexes[i].name) + 1) * sizeof(char));
        namcx[i] = strcpy(namcx[i], complexes[i].name);

        // fprintf(stderr,"initialize_cpp: species %s (%d) charge: %f, keq: %f, %d %f\n",
        // 	      namcx[i],ncomp+i,chg[ncomp+i],complexes[i].keq,ntemp,temp);

        // equilibrium constant stuff
        if(ntemp > 1)
        {
          fit(&nbasis,&ntemp,complexes[i].keq,bvec,&vec[0][0],&iflgint,INT,&ntt,namcx[i]);

          if(iflgint == 1)
          {
            nncnt = 0;
            for(int j = 0; j < ntemp; ++j)
              if (INT[j] == 1)
              {
                nncnt = nncnt + 1;
                temptmp[nncnt-1] = tempc[j];
              }

            if (ntt == 1)
              if (temp != temptmp[0])
              {
                printf(" Only one logK in database at T(C): %f\n", temptmp[0]);
                printf(" Temperature of condition:         %f\n",temp);
                printf(" Species missing log Ks: %s\n", namcx[i]);
                printf(" Exiting due to input file or parsing error.\n");
                exit(-1);
              }
          }

          // put it in the permanent vector of coefficients
          for(int j = 0; j < NBASIS; ++j)
            // double as1[NBASIS][12]; // transposed from fortran as1(12,5)
            as1[j*nspec+i] = bvec[j];
        }
        else if (ntemp == 1)
          eqhom[i] = complexes[i].keq[0];
      } // for(int i = 0; i < nspec; i++)
    } // COMPLEXATION

      //     ______________________ MINERALS _________________________    
      
      //     pre-dimension
    np = MAX_PATH; // global, needed by reaction_cpp_
    int max_ip;
    int ip;
    int id; // integer(i4b) :: id,nd
    double time_scale;

    //     read all minerals
    if (nkin < 0)
    {
      printf("negative minerals\n");
      exit(-1);
    }
    else if (nkin == 0)
      printf("no minerals\n");
    else 
    {
      //     (pre-)allocate variables
      np = MAX_PATH; // max number of paths for all mineral reactions
      max_ip = 0;
      nd = ncomp+nspec; // max number of dependencies

      for(size_t i = 0; i < nkin; ++i)
        for(size_t j = 0; j < np; ++j)
          AffinityDepend1[i][j] = 1.0; // initialize AffinityDepend1

      int nptot = 0;

      // loop over reactions
      for(int nkinIndex = 0; nkinIndex < nkin; nkinIndex++)
      {
        // initialize mumin
        for(int ncompIndex = 0; ncompIndex < ncomp; ++ncompIndex)
          for(int npIndex = 0; npIndex < np; ++npIndex)
            mumin[ncompIndex*nkin*np+nkinIndex*np+npIndex] = 0.0; // transpose of (np,nkin,ncomp), mumin(:,i,:) = 0.0d0

        // example mineral reaction stoichiometry: 1.0 'Ca++'  -2.0 'H+'  1.0 'CO2(aq)'
        char *tokenPtr = strtok(minerals[nkinIndex].stoichiometry, " "); // begin tokenizing stoichiometry
        int tokenNum = 0;
        char *p;
        char *q;
        double mu;

        // continue tokenizing line until tokenPtr becomes NULL
        while (tokenPtr != NULL) 
        {
          if(tokenNum % 2 == 0)
          {
            // printf("mineral stoichiometry: mu %s\n", tokenPtr);
            sscanf(tokenPtr,"%lf",&mu);
          }
          else
          {
            // printf("mineral stoichiometry: component %s\n", tokenPtr);
            if((p = strpbrk(tokenPtr,"'")))
            {
              q = strrchr(p,'\'');
              if(q)
                *q = '\0';
              if(!strcmp(p+1, "default"))
              {
                break;
              }
              else 
              {
                int j;
                for(j = 0; j < ncomp; ++j)
                  if(!strcmp(p+1, ulab[j]))
                  {
                    for(int k = 0; k < np; ++k)
                    {
                      mumin[j*nkin*np+nkinIndex*np+k] = mu;
                      // fprintf(stderr,"mumin[%d][%d][%d] = %.4e\n",j,nkinIndex,k,mu);
                    }
                    break;
                  }
                if (j == ncomp)
                {
                  printf("mumin: could not find component %s for mineral reaction %s\n",
                         p+1, minerals[nkinIndex].name);
                  exit(-1);
                }
              }
            }
          }
          tokenNum++;
          tokenPtr = strtok(NULL, " "); // get next token
        }

        namin[nkinIndex] = minerals[nkinIndex].name;
        ip = 0;

        // volume of minerals
        volmol[nkinIndex] = minerals[nkinIndex].molar_volume * 1.0e-6;

        // equilibrium constant stuff
        if (ntemp > 1)
        {
          //   void fit_(int *nbasis, int *ntemp, double *alogk0, double *bvecF, double *vecF, int *iflgint, int *inoutint, int *ntt, char *nameTransfer)
          // fprintf(stderr,"fit: cpp %s, %d\n",namin[nkinIndex],nkinIndex);
          fit(&nbasis, &ntemp, minerals[nkinIndex].keq, bvec, &vec[0][0], &iflgint, INT, &ntt, namin[nkinIndex]);

          if (iflgint == 1)
          {
            nncnt = 0;
            for(int j = 0; j < ntemp; ++j)
            {
              if (INT[j] == 1)
              {
                nncnt = nncnt + 1;
                temptmp[nncnt] = tempc[nkinIndex];
              }
            }

            if (ntt == 1)
            {
              if (temp != temptmp[0])
              {
                printf("Only one logK in database at T(C): %f\n",temptmp[0]);
                printf(" Temperature of condition:          %f\n",temp);
                printf(" Species missing log Ks: %s\n", namin[nkinIndex]);
                printf("\n");
                exit(-1);
              }
            }
          } // if (iflgint == 1)
        } // if (ntemp > 1)

        // end of equilibrium constant stuff

      } // for(int nkinIndex = 0; nkinIndex < nkin; nkinIndex++), loop over reactions
    } // read all minerals


    // figure out max dimensions
    nd = maxval(nkin, np, ndepend);
    np = max_ip;
      
    // reallocate things
    ndim1 = np;
    ndim2 = nkin;
    ndim3 = ncomp;

    // fprintf(stderr,
    // 	      "initialize.cpp: reallocate mumin to dimension [%d][%d][%d]\n",
    // 	      ndim1,ndim2,ndim3);

    // printarray("mumin", 6, 1, 10, mumin); // transpose of (np,nkin,ncomp)

    // stoichchombo(1,1:nkin,1:ncomp) = mumin(1,1:nkin,1:ncomp)
    for(int i = 0; i < ncomp; ++i)
      for(int j = 0; j < nkin; ++j)
        for(int k = 0; k < np; ++k)
          size_t idx = i*nkin*np+j*np+k;
          stoichchombo[idx] = mumin[idx];
  
    for(int k = 0; k < nz; ++k)
      for(int j = 0; j < ny; ++j)
        for(int i = 0; i < nx; ++i)
          for(int l = 0; l < neqn; ++l)
          {
            size_t idx = k*ny*nx*neqn+j*nx*neqn*i*neqn+l;
            s[idx] = 0.0; // allocate(s(neqn,nx,ny,nz),stat=ierr); s = 0.0d0
            sn[idx] = 0.0; // allocate(sn(neqn,nx,ny,nz),stat=ierr); sn = 0.0d0
          }

    //
    for(int k = 0; k < nz; ++k)
      for(int j = 0; j < ny; ++j)
        for(int i = 0; i < nx; ++i)
        {
          size_t idx = k*ny*nx+j*nx+i;
          sion[idx] = 0.0;
          t[idx] = 25.0;
        }


    // read all solutions   
    for(unsigned int i=0; i < MAX_SOL; ++i)
    {
      if(solutions[i].name == NULL)
        break; // no more solutions

      Species *species = solutions[i].concentrationList; // head node 

      for(int k = 0; k < ncomp + nspec; ++k)
      {
        // kth dependency in linked list
        //
        if(species == NULL)
          break; // end of list
        else if (!strcmp(species->name, "default"))
          break;  // tail 
        int j;
        for(j = 0; j < ncomp + nspec; ++j)
        {
          if (j < ncomp)
          {
            if (!strcmp(species->name, ulab[j]))
            {
              sp10[j] = species->nu;
              // fprintf(stderr,"sp10[0][0][0][%d] = %.3e (%s)\n",j,species->nu,ulab[j]);
              break;
            }
          }
          else if (j >= ncomp) 
          {   
            if (!strcmp(species->name, namcx[j-ncomp]))
            {
              sp10[j] = species->nu;
              // fprintf(stderr,"sp10[0][0][0][%d] = %.3e (%s)\n",j,species->nu,namcx[j-ncomp]);
              break;
            }
          }
        }

        if (j >= ncomp+nspec)
        {
          fprintf(stderr,"could not find species %s in solution %s\n",
                  species->name, solutions[i].name);
          exit(1);
        }
	      
        if(species->next != NULL)
          species = species->next;
        else
          break; // no more species in list
      }

      t = temp;

      if (!strcmp(solutions[i].name,"inlet"))
      {
        // this is the inlet boundary

        cfgamma(ncomp,
                nspec,
                0,
                0,
                0,
                nx,
                ny,
                sp10,
                t,
                adhcoeff,
                bdhcoeff,
                bdtcoeff,
                sion,
                ulab,
                acmp,
                chg,
                gam);

        // int jx = 0; 
        // int jy = 0; 
        // int jz = 0;


        totconc_cpp_(ncomp, 
                     nspec, 
                     neqn,
                     1, 
                     1, 
                     1, 
                     nx, 
                     ny, 
                     nz, 
                     muaq,
                     sp10,
                     s);



        for(size_t i = 0; i < ncomp; ++i)
        {
          inflowchombo[i] = s[i]/1000.0;
          inflowchombosp10[i] = sp10[i]/1000.0;
        }


        for(int k = 0; k < nz; ++k)
          for(int j = 0; j < ny; ++j)
            for(int i = 0; i < nx; ++i)
              for(int l = 0; l < neqn; ++l)
                s[k*(ny*nx*neqn)+j*(nx*neqn)+i*neqn+l] = 0.0; // allocate(s(neqn,nx,ny,nz),stat=ierr); s = 0.0d0

        for(int k = 0; k < nz; ++k)
          for(int j = 0; j < ny; ++j)
            for(int i = 0; i < nx; ++i)
              for(int l = 0; l < ncomp + nspec; ++l)
	      {
              int nc = ncomp+nspec;
              size_t idx = k*(ny*nx*nc)+j*(nx*nc)+i*nc+l;
              sp10[idx] = 0.0; // allocate(sp10(ncomp+nspec,nx,ny,nz),stat=ierr); sp10 = 0.0d0
              gam[idx] = 0.0; // allocate(gam(ncomp+nspec,nx,ny,nz),stat=ierr); gam = 0.0d0
	      }

      }
      else if (!strcmp(solutions[i].name,"inlet2"))
      {
        // this is a second the inlet boundary

        cfgamma(ncomp,
                nspec,
                0,
                0,
                0,
                nx,
                ny,
                sp10,
                t,
                adhcoeff,
                bdhcoeff,
                bdtcoeff,
                sion,
                ulab,
                acmp,
                chg,
                gam);

        // int jx = 0; 
        // int jy = 0; 
        // int jz = 0;

        totconc_cpp_(ncomp, 
                     nspec, 
                     neqn,
                     1, 
                     1, 
                     1, 
                     nx, 
                     ny, 
                     nz, 
                     muaq,
                     sp10,
                     s);

        for(int k = 0; k < nz; ++k)
          for(int j = 0; j < ny; ++j)
            for(int i = 0; i < nx; ++i)
              for(int l = 0; l < neqn; ++l)
                s[k*(ny*nx*neqn)+j*(nx*neqn)+i*neqn+l] = 0.0; // allocate(s(neqn,nx,ny,nz),stat=ierr); s = 0.0d0

        for(int k = 0; k < nz; ++k)
          for(int j = 0; j < ny; ++j)
            for(int i = 0; i < nx; ++i)
              for(int l = 0; l < ncomp + nspec; ++l)
              {
                int nc = ncomp+nspec;
                size_t idx = k*(ny*nx*nc)+j*(nx*nc)+i*nc+l;
                sp10[idx] = 0.0; // allocate(sp10(ncomp+nspec,nx,ny,nz),stat=ierr); sp10 = 0.0d0
                gam[idx] = 0.0; // allocate(gam(ncomp+nspec,nx,ny,nz),stat=ierr); gam = 0.0d0
              }
	      
      }
      else if (!strcmp(solutions[i].name,"initial"))
      {
        //  this is the initial condition in the domain

        cfgamma(ncomp,
                nspec,
                0,
                0,
                0,
                nx,
                ny,
                sp10,
                t,
                adhcoeff,
                bdhcoeff,
                bdtcoeff,
                sion,
                ulab,
                acmp,
                chg,
                gam);



        totconc_cpp_(ncomp, 
                     nspec, 
                     neqn,
                     1, 
                     1, 
                     1, 
                     nx, 
                     ny, 
                     nz, 
                     muaq,
                     sp10,
                     s);


        for(size_t i = 0; i < ncomp; ++i)
        {
          initchombo[i] = s[i]/1000.0;
        }

        for(size_t i = 0; i < ncomp + nspec; ++i)
        {
          initchombosp10[i] = sp10[i];
          initchombosp[i] = log(sp10[i]);
        }

        for(size_t i = 0; i < ncompchombo + nspecchombo; ++i)
          initgamma[i] = gam[i];

        // reset to zero
        //
        for(size_t i = 0; i < neqn; ++i)
          s[i] = 0.0; // allocate(s(neqn,nx,ny,nz),stat=ierr); s = 0.0d0

        for(size_t i = 0; i < ncompchombo + nspecchombo; ++i)
        {
          sp10[i] = 0.0; // allocate(sp10(ncomp+nspec,nx,ny,nz),stat=ierr); sp10 = 0.0d0
          gam[i] = 0.0; // allocate(gam(ncomp+nspec,nx,ny,nz),stat=ierr); gam = 0.0d0
        }
      }
      else 
      {     
        fprintf(stderr,"name of condition %s not provided or recognized\n",
                solutions[i].name);
        exit(1);
      }
    } // read all solutions

      // retrieve initial concentrations from 2,1,1 to entire domain
    for(size_t k = 0; k < nz; ++k)
      for(size_t j = 0; j < ny; ++j)
        for(size_t i = 0; i < nx; ++i)
        {
          for(size_t l = 0; l < neqn; ++l)
            s[k*ny*nx*neqn+j*nx*neqn+i*neqn+l] = initchombo[l] * 1000.0; // s(1:neqn,2,1,1), allocate(s(neqn,nx,ny,nz),stat=ierr); s = 0.0d0
          for(size_t l = 0; l < ncompchombo + nspecchombo; ++l)
          {
            int nc = ncompchombo+nspecchombo;
            size_t idx = k*(ny*nx*nc)+j*(nx*nc)+i*nc+l;
            sp10[idx] = initchombosp10[l]; // sp10(1:ncomp+nspec,2,1,1)
            gam[idx]  = initgamma[l]; // gam(1:ncomp+nspec,2,1,1)
          }
        }

    // initialize sn, and sp
    int nc = ncompchombo + nspecchombo;
    for(int k = 0; k < nz; ++k)
      for(int j = 0; j < ny; ++j)
        for(int i = 0; i < nx; ++i)
        {
          for(int l = 0; l < neqn; ++l) {
            size_t idx = k*(ny*nx*neqn)+j*(nx*neqn)+i*neqn+l;
            sn[idx] = s[idx];
          }
          for(int l = 0; l < ncompchombo + nspecchombo; ++l) {
            size_t idx = k*(ny*nx*nc)+j*(nx*nc)+i*nc+l;
            sp[idx] = log(sp10[idx]);
          }
        }

    // ______________________ HARDWIRED ________________________
    return(0);
  }


  char *trimwhitespace(char *str)
  {
    char *end;

    // Trim leading space
    while(isspace((unsigned char)*str)) str++;

    if(*str == 0)  // All spaces?
      return str;

    // Trim trailing space
    end = str + strlen(str) - 1;
    while(end > str && isspace((unsigned char)*end)) end--;

    // Write new null terminator character
    end[1] = '\0';

    return str;
  }


#ifdef __cplusplus
}
#endif
