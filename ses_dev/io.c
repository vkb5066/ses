//Implementation of all file input/output functions
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include "def.h"

#define LINESIZE_MAX 128u
#define N_TOKEN_MAX 32u
#define DELIM " \t"
#define BASE 10


//---Input----------------------------------------------------------------------
/*
*  reads string literal INFILE_RUNPARAMS and returns a job object, with all
*  elements set (either by user input or defaults)
*  MEM: none
*  RET: 0 on success, 1 on file not found error
*/
uchar readjob(job*restrict runparams){
    FILE* infile;
    int c;
    char line[LINESIZE_MAX];
    char* next;    
    uint i;
    uchar eq, j;

    if(!(infile = fopen(INFILE_RUNPARAMS, "r"))) return (uchar)1;    
    

    //need to set some defaults as flags to later see if any given
    //tag was initialized or not.  marked members can have sensible defaults
    //set right now, so I do that
    runparams->runmode = (uchar)1; ///okay default (to static bs calc)
    runparams->diagmode = (uchar)2; ///okay default (to davidson) 
    runparams->nthreads = (ushort)1; ///okay default 
    runparams->nbands = (ushort)0; ///needs updated if no estimate is read in
    runparams->itrlim = (ushort)0; ///needs updated if no value read in
    runparams->stab = (uchar)1; ///okay default
    runparams->meff = ONE; ///okay default 
    runparams->nfstargets = (uchar)0; ///okay default
    runparams->entol = -ONE; ///needs updated if no value is read in 
    runparams->encut = -ONE; ///needs updated if no value is read in
    runparams->fftgmul = ONE; ///okay default 
    runparams->timestep = ONE; ///okay default
    runparams->mditrlim = DEF_MD_ITRLIM; ///okay default
    runparams->rcut = (fp)3.75; ///okay default
    runparams->acut = (fp)120.0; ///okay default
    runparams->ftol = DEF_MD_FTOL; ///okay default

    //reading file starts here
    while((c = getc(infile)) != EOF){
        eq = (uchar)0;
        fgets(line+1, LINESIZE_MAX-1u, infile);
        line[0] = (char)c;

        //comment support: want to split the read-in string at the first 
        //instance of a comment: in this case, hard-coded to '#'
        for(i = 0u; i < LINESIZE_MAX; ++i){
            if(line[i] == '=') eq = (uchar)1;
            if(line[i] == '#'){
                ///if there is a comment before an equals sign, we don't need 
                ///any further parsing
                if(!eq) goto endofloop;
                ///otherwise, just ignore anything after the comment and move on
                line[i] = '\0';
                break;
            }
        }       

        next = strtok(line, DELIM"=\n"); ///next = the 'key' to map to job vars
        
        //Map key to job Vars
        ///algo control
        if(!strcmp("mrun", next)){
            next = strtok(NULL, DELIM"=");
            runparams->runmode = (uchar)atoi(next);
            continue;
        }
        if(!strcmp("mdiag", next)){
            next = strtok(NULL, DELIM"=");
            runparams->diagmode = (uchar)atoi(next);
            continue;
        }
        if(!strcmp("nthreads", next)){
            next = strtok(NULL, DELIM"=");
            runparams->nthreads = (ushort)atoi(next);
            continue;
        }
        if(!strcmp("nbands", next)){
            next = strtok(NULL, DELIM"=");
            runparams->nbands = (ushort)atoi(next);
            continue;
        }
        if(!strcmp("ilim", next)){
            next = strtok(NULL, DELIM"=");
            runparams->itrlim = (ushort)atoi(next);
            continue;
        }
        if(!strcmp("lcgstab", next)){
            next = strtok(NULL, DELIM"=");
            runparams->stab = (uchar)atoi(next);
            continue;
        }
        if(!strcmp("meff", next)){
            next = strtok(NULL, DELIM"=");
            runparams->meff = (fp)atof(next);
            continue;
        }
        if(!strcmp("fstarget", next)){
            runparams->fstargets = malloc(N_TOKEN_MAX*sizeof(fp));
            for(j = (uchar)0; (next = strtok(NULL, DELIM"=")); ++j){
                runparams->fstargets[j] = (fp)atof(next);
                (runparams->nfstargets)++;
            }
            runparams->fstargets = realloc(runparams->fstargets, 
                                           runparams->nfstargets*sizeof(fp));
            continue;
        }
        ///accuracy control
        if(!strcmp("entol", next)){
            next = strtok(NULL, DELIM"=");
            runparams->entol = (fp)atof(next);
            continue;
        }
        if(!strcmp("encut", next)){
            next = strtok(NULL, DELIM"=");
            runparams->encut = (fp)atof(next);
            continue;
        }
        if(!strcmp("fftbsmul", next)){
            next = strtok(NULL, DELIM"=");
            runparams->fftgmul = (fp)atof(next);
            continue;
        }
        ///md control
        if(!strcmp("timestep", next)){
            next = strtok(NULL, DELIM"=");
            runparams->timestep = (fp)atof(next);
            continue;
        }
        if(!strcmp("milim", next)){
            next = strtok(NULL, DELIM"=");
            runparams->mditrlim = (ushort)atoi(next);
            continue;
        }
        if(!strcmp("rcut", next)){
            next = strtok(NULL, DELIM"=");
            runparams->rcut = (fp)atof(next);
            continue;
        }
        if(!strcmp("acut", next)){
            next = strtok(NULL, DELIM"=");
            runparams->acut = (fp)atof(next);
            continue;
        }
        if(!strcmp("ftol", next)){
            next = strtok(NULL, DELIM"=");
            runparams->ftol = (fp)atof(next);
            continue;
        }
        ///read and write
/*        if(!strcmp("rpsi", next)){
            next = strtok(NULL, DELIM"=");
            runparams->readpsi = (uchar)atoi(next);
            continue;
        }
        if(!strcmp("wpsi", next)){
            next = strtok(NULL, DELIM"=");
            runparams->writepsi = (uchar)atoi(next);
            continue;
        }
*/

        endofloop:;
    }

    fclose(infile);
    return (uchar)0;
}

/*
*  reads string literal 'INFILE_LATTICE' and returns a lattice object, with
*  A[3][3], nspecs, speccounts set and 'sites' filled
*  MEM: allocation of nsites sites and speccounts 
*  RET: 0 on success, 1 on file not found error
*/
uchar readlattice(lat*restrict lattice){
    //File format is as follows (basically a simple vasp.5 POSCAR file):
    //line 0:     comment (no-op)
    //lines 1-3:  A, read in row-vector order
    //line 4:     list of species
    //line 5:     counts of species, in line with line 4
    //line 6+:    direct positions of atoms a, b, c

    FILE* infile;
    char line[LINESIZE_MAX];
    char* next;    
    uint i, j, acc;

    if(!(infile = fopen(INFILE_LATTICE, "r"))) return (uchar)1;

    //Reading file starts here
    fgets(line, LINESIZE_MAX, infile); ///skip comment line
    
    ///lattice input
    for(i = 0u; i < 3u; ++i){
        fscanf(infile, " "FPFORMAT" "FPFORMAT" "FPFORMAT" ", 
               &lattice->A[i][0], &lattice->A[i][1], &lattice->A[i][2]);
    }

    ///element types (skip for now, just need the count of them)   
    fgets(line, LINESIZE_MAX, infile);
    next = strtok(line, DELIM); lattice->nspecs = (uchar)1;
    while((next = strtok(NULL, DELIM))){
        lattice->nspecs++;
    }

    ///element counts
    lattice->speccounts = malloc((lattice->nspecs)*sizeof(ushort));
    for(i = 0u, acc=0u; i < (uint)(lattice->nspecs); ++i){
        fscanf(infile, " %hu ", &lattice->speccounts[i]);
        acc += (uint)(lattice->speccounts[i]);
    }

    ///lattice sites (only set species and fractional coordinates)
    lattice->sites = malloc(acc*sizeof(site));
    for(i = 0u, acc = 0u; i < (uint)lattice->nspecs; ++i){
        for(j = 0u; j < (uint)lattice->speccounts[i]; ++j, ++acc){
            lattice->sites[acc].spec = malloc(sizeof(uchar));
            *(lattice->sites[acc].spec) = (uchar)i;

            fscanf(infile, " "FPFORMAT" "FPFORMAT" "FPFORMAT" ",
                   &lattice->sites[acc].crdsf[0], &lattice->sites[acc].crdsf[1],
                   &lattice->sites[acc].crdsf[2]);
        }
    }

    
    fclose(infile);
    return (uchar)0;
}

/*
*  reads string literal 'INFILE_POTENTIAL' and returns an array of potentials
*  MEM: allocation of nspecs potentials nspecs potential samples
*  RET: 0 on success, 1 on file not found error
*/
uchar readpotential(const uchar nspecs, pot*restrict*restrict potentials){
     //File format is as follows:
    //line 0:   comment (no-op)
    //line 1:   number of val electrons, mass in amu
    //line 1:   q^2 sample start, q^2 sample stop, num of q^2 steps
    //line 2:   q^2 samples
    //line 3+:  restart from 0 for each 1 <= i <= nspecs

    FILE* infile;
    char line[LINESIZE_MAX];
    uchar i;
    ushort j;

    if(!(infile = fopen(INFILE_POTENTIAL, "r"))) return (uchar)1;

    *potentials = malloc(nspecs*sizeof(pot));
    //Reading file starts here: one block for each species
    for(i = (uchar)0; i < nspecs; ++i){
        fgets(line, LINESIZE_MAX, infile); ///skip comment line
        
        ///num val elecs, atomic mass
        fscanf(infile, " %hhu "FPFORMAT" ", 
               &(*potentials)[i].nvel, &(*potentials)[i].mass);

        ///start, stop, step of pot sample array
        fscanf(infile, " "FPFORMAT" "FPFORMAT" %hu ",
               &(*potentials)[i].q2lo, &(*potentials)[i].q2hi, 
               &(*potentials)[i].nsamples);

        ///potential array
        (*potentials)[i].samples = malloc((*potentials)[i].nsamples*sizeof(fp));
        for(j = (ushort)0; j < (*potentials)[i].nsamples; ++j){
            fscanf(infile, " "FPFORMAT" ", &(*potentials)[i].samples[j]);
        }
    }


    fclose(infile);
    return (uchar)0;
}


/*
*  reads string literal 'INFILE_KPOINTS' and returns an array of kpts
*  MEM: allocation of nkpts k-points
*  RET: 0 on success, 1 on file not found error, 2 on unrecognized mode error
*/
uchar readkpoints(const fp B[restrict 3][3],
                  ushort*restrict nkpoints, kpt*restrict*restrict kpoints){
    //File format is as follows (a bit like a vasp.5 KPOINTS file):
    //line 0:   comment (no-op)
    //line 1:   a character, 'e', 'l', 'm', for 'explicit', 'line', 'mesh'
    //line 2:   numbers ("nkpts" if 'e', "nlines" "nkpts/line" if 'l',
    //          "nkpts along b1" "nkpts along b2" "nkpts along b3" if 'm'
    //line 3+:  depends on line 1 ... see further down in the function

    FILE* infile;
    char line[LINESIZE_MAX];
    char mode;
    ushort i, j, i1, i2, i3, subdivs[3];
    fp jp;
    fp d[3];

    if(!(infile = fopen(INFILE_KPOINTS, "r"))) return (uchar)1;

    //Reading file starts here
    fgets(line, LINESIZE_MAX, infile); ///skip comment line

    ///figure out what mode we're in
    fscanf(infile, " %c ", &mode);
    
    ///here the reading / working with kpoints changes depending on mode
    switch(mode){
        ///'explicit' mode: lines 3+ each have the k point coordinates in 
        ///fractional units along b1, b2, b3, followed by the k-p weight
        case 'e':
            fscanf(infile, " %hu ", nkpoints);
            
            *kpoints = malloc((*nkpoints)*sizeof(kpt));
            for(i = (ushort)0; i < *nkpoints; ++i){
                fscanf(infile, " "FPFORMAT" "FPFORMAT" "FPFORMAT" "FPFORMAT" ",
                       &(*kpoints)[i].crds[0], &(*kpoints)[i].crds[1],
                       &(*kpoints)[i].crds[2], &(*kpoints)[i].wgt);
            }
            break;

        ///'line' mode: lines 3+ each contain triplets of direct kpoint coords
        ///that are interpolated between a-la numpy's 'linspace' function
        case 'l':
            fscanf(infile, " %hu %hu ", &i1, &i2);
            *nkpoints = i1*i2;

            *kpoints = malloc((*nkpoints)*sizeof(kpt));
            ///over all lines
            for(i = (ushort)0, i3 = (ushort)0; i < i1; ++i, i3 += i2){
                fscanf(infile, " "FPFORMAT" "FPFORMAT" "FPFORMAT
                               " "FPFORMAT" "FPFORMAT" "FPFORMAT" ",
                       &(*kpoints)[i3].crds[0],               // |
                       &(*kpoints)[i3].crds[1],               //  > start
                       &(*kpoints)[i3].crds[2],               // |
                       &(*kpoints)[i2+i3-(ushort)1].crds[0],  // |
                       &(*kpoints)[i2+i3-(ushort)1].crds[1],  //  > end
                       &(*kpoints)[i2+i3-(ushort)1].crds[2]); // |             
                
                ////intermediate
                d[0] = ((*kpoints)[i2+i3-(ushort)1].crds[0] - 
                         (*kpoints)[i3].crds[0]) / ((fp)(i2 - (ushort)1));
                d[1] = ((*kpoints)[i2+i3-(ushort)1].crds[1] - 
                         (*kpoints)[i3].crds[1]) / ((fp)(i2 - (ushort)1));
                d[2] = ((*kpoints)[i2+i3-(ushort)1].crds[2] - 
                         (*kpoints)[i3].crds[2]) / ((fp)(i2 - (ushort)1));
                for(j = i3 + (ushort)1, jp = ONE; j < i2 + i3 - (ushort)1; 
                    ++j, jp += ONE){
                    (*kpoints)[j].crds[0] = (*kpoints)[i3].crds[0] + jp*d[0];
                    (*kpoints)[j].crds[1] = (*kpoints)[i3].crds[1] + jp*d[1];
                    (*kpoints)[j].crds[2] = (*kpoints)[i3].crds[2] + jp*d[2];
                }
            }
            ///may as well set appropriate weights, though you'd be a moron for
            ///trying to use a linemode mesh for integrals
            jp = ONE / (fp)(*nkpoints);
            for(i = (ushort)0; i < *nkpoints; ++i){
                (*kpoints)[i].wgt = jp;
            }
            break;

        ///'mesh' mode: line 3 is a triplet of integers that represent the 
        ///number of subdivisions in the (h, k, l) directions a la M.P.
        ///this implementation uses 0-indexing and applies inversion symmetry
        case 'm':
            fscanf(infile, " %hu %hu %hu ", &i1, &i2, &i3);
            *nkpoints = i1*i2*i3/(ushort)2;
            if(i1%(ushort)2 != (ushort)0 && i2%(ushort)2 != (ushort)0 && 
               i3%(ushort)2 != (ushort)0){
                (*nkpoints)++; ///gamma point
            }
            *kpoints = malloc((*nkpoints)*sizeof(kpt));
            subdivs[0] = i1; subdivs[1] = i2; subdivs[2] = i3;
            mpgridfi(B, subdivs, *kpoints);
            break;

        default:
            return (uchar)2;
    }

    fclose(infile);
    return (uchar)0;
}


/*
*  reads string literal 'INFILE_KEATING' and returns an array of keating params
*  WRN: maximum number of species allowed is floor(sqrt(uchar_max))
*  MEM: allocation of nspecs keating objects
*  RET: 0 on success, 1 on file not found error
*/
uchar readkeating(const uchar nspecs, keat*restrict*restrict keatings){
     //File format is as follows:
    //line 0:   comment (no-op)
    //line 1:   r_0(0, 0), r_0(0, 1), ..., r_0(0, nspecs-1)
    //line 2:   alpha(0, 0), alpha(0, 1), ..., alpha(0, nspecs-1)
    //line 3:   theta_0(0, 0, 0), theta_0(0, 0, 1), ..., theta_0(0, 0, nspecs-1)
    //          theta_0(0, 1, 0), theta_0(0, 1, 1), ..., theta_0(0, 1, nspecs-1)
    //          ...
    //          theta_0(0, nspecs-1, 0), ... theta_0(0, nspecs-1, nspecs-1)
    //line 4:   beta(0, 0), beta(0, 1), ..., beta(0, nspecs-1)
    //line 5+:  restart at line 0, increment first index by 1
    //units: angstroms, eV/angstroms^2, and degrees

    FILE* infile;
    char line[LINESIZE_MAX];
    uchar i, j, nspecs2;

    if(!(infile = fopen(INFILE_KEATING, "r"))) return (uchar)1;

    nspecs2 = nspecs*nspecs;
    *keatings = malloc(nspecs*sizeof(keat));
    //reading file starts here: one block of four lines for each species
    for(i = (uchar)0; i < nspecs; ++i){
        fgets(line, LINESIZE_MAX, infile); ///skip comment line

        //ideal interatomic spacing (dist in angstrm to inverse dist)
        (*keatings)[i].r0inv = malloc(nspecs*sizeof(fp));
        for(j = (uchar)0; j < nspecs; ++j){
            fscanf(infile, " "FPFORMAT" ", &(*keatings)[i].r0inv[j]);
            (*keatings)[i].r0inv[j] = ONE / (*keatings)[i].r0inv[j];
        }

        //alpha params
        (*keatings)[i].al = malloc(nspecs*sizeof(fp));
        for(j = (uchar)0; j < nspecs; ++j){
            fscanf(infile, " "FPFORMAT" ", &(*keatings)[i].al[j]);
        }

        //ideal interatomic angles (angle in deg to cos(angle in rad)
        (*keatings)[i].c0 = malloc(nspecs2*sizeof(fp));
        for(j = (uchar)0; j < nspecs2; ++j){
            fscanf(infile, " "FPFORMAT" ", &(*keatings)[i].c0[j]);
            (*keatings)[i].c0[j] = COS(DEG2RAD((*keatings)[i].c0[j]));
        }

        //beta params (val in eV/A^2 to its square root)
        (*keatings)[i].sqrtbe = malloc(nspecs*sizeof(fp));
        for(j = (uchar)0; j < nspecs; ++j){
            fscanf(infile, " "FPFORMAT" ", &(*keatings)[i].sqrtbe[j]);
            (*keatings)[i].sqrtbe[j] = SQRT((*keatings)[i].sqrtbe[j]);
        }
    }

    fclose(infile);
    return (uchar)0;
}

//---Output (to stdout)---------------------------------------------------------
//Write out basic information (version, date, ...)
//a smart person will redirct this to a file so that it can be mapped to results
//later, but even if that does not happen, it looks professional
void writeheader(const uint seed){
    time_t t;
    struct tm tm;

    //version
    puts("ses v0.2");
    //date, time
    t = time(NULL);
    tm = *localtime(&t);
    printf("date: %d-%02d-%02d %02d:%02d:%02d\n", 
           tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, 
           tm.tm_hour, tm.tm_min, tm.tm_sec);
    //precision level
    #if(PREC == 1)
    puts("using SINGLE point precision");
    #elif(PREC == 2)
    puts("using DOUBLE point precision");
    #endif
    //random seed
    printf("rng seed: %u\n", seed);

    puts("");
}



//---Output (to file)-----------------------------------------------------------
#define FFMT "%+.7e"
//Writes out the eigenvalues vs. k for a set of kpoints to 'OUTFILE_EIGENVALS'
void writebands(const uchar id,
                const fp B[restrict 3][3], const ushort nbands, 
                const ushort nkpts, const kpt*restrict kpts){
    //File format is as follows (a bit like a vasp.5 EIGENVAL file):
    //line 0:   b1x b1y b1z
    //line 1:   b2x b2y b2z
    //line 2:   b3x b3y b3z
    //line 3:   nkpts nbands_per_kpt
    //line 4+: the kpoints themselves, written as
    //         'weight k1 k2 k3 \n band 0 energy \n band 1 energy \n ...
    //there are no empty lines anywhere within the file, kx are in whatever 
    //units they were sent in as (prob. fractional)

    FILE* outfile;
    char outname[N_TOKEN_MAX];
    ushort i, j;

    for(i = (ushort)0; OUTFILE_EIGENVAL_BASE[i] != '\0'; ++i){
        outname[i] = OUTFILE_EIGENVAL_BASE[i];
    }
    outname[i] = id;
    outfile = fopen(outname, "w");

    //write the reciprocal space lattice
    fprintf(outfile, FFMT" "FFMT" "FFMT"\n", B[0][0], B[0][1], B[0][2]);
    fprintf(outfile, FFMT" "FFMT" "FFMT"\n", B[1][0], B[1][1], B[1][2]);
    fprintf(outfile, FFMT" "FFMT" "FFMT"\n", B[2][0], B[2][1], B[2][2]);

    //write counts for parsing later
    fprintf(outfile, "%hu %hu\n", nkpts, nbands);

    //write the k-points
    for(i = (ushort)0; i < nkpts; ++i){
        fprintf(outfile, FFMT" "FFMT" "FFMT" "FFMT"\n", 
                kpts[i].wgt, kpts[i].crds[0], kpts[i].crds[1], kpts[i].crds[2]);
        for(j = (ushort)0; j < nbands; ++j){
            fprintf(outfile, FFMT"\n", kpts[i].evals[j]);
        }
    }

    fclose(outfile);
}
#undef FFMT

//writes lattice object to file OUTFILE_LATTICE
//if specmap is not NULL, will replace the atom species line by specmap
void writelattice(const lat lattice, const char*restrict specmap){
    //file format is exactly the same as the infile format except the atomic 
    //species (written as a, b, c, ... instead of species)

    FILE* outfile;
    uchar i; ushort j, acc;

    outfile = fopen(OUTFILE_LATTICE, "w");

    //comment
    fprintf(outfile, "written by writelattice()\n");
    
    //basis
    fprintf(outfile, FPFORMAT" "FPFORMAT" "FPFORMAT"\n", 
            lattice.A[0][0], lattice.A[0][1], lattice.A[0][2]); 
    fprintf(outfile, FPFORMAT" "FPFORMAT" "FPFORMAT"\n", 
            lattice.A[1][0], lattice.A[1][1], lattice.A[1][2]); 
    fprintf(outfile, FPFORMAT" "FPFORMAT" "FPFORMAT"\n", 
            lattice.A[2][0], lattice.A[2][1], lattice.A[2][2]); 

    //atom species, counts
    if(specmap) fprintf(outfile, "%s", specmap);
    else for(i = (uchar)0; i < lattice.nspecs; ++i) 
             fprintf(outfile, "%c ", 'a' + i); 
    fprintf(outfile, "\n");
    for(i = (uchar)0; i < lattice.nspecs; ++i){
        fprintf(outfile, "%hu ", lattice.speccounts[i]);
    } fprintf(outfile, "\n");

    //atoms
    for(i = (uchar)0, acc = 0u, acc = 0u; i < lattice.nspecs; ++i){
    for(j = (ushort)0; j < lattice.speccounts[i]; ++j, ++acc){
        fprintf(outfile, FPFORMAT" "FPFORMAT" "FPFORMAT"\n", 
                lattice.sites[acc].crdsf[0], lattice.sites[acc].crdsf[1], 
                lattice.sites[acc].crdsf[2]); 
    }
    }
}







#undef LINESIZE_MAX
#undef N_TOKEN_MAX
#undef DELIM
#undef BASE
