#ifndef FRAMEWORK_VERSION
//static const char compileDate[] = __DATE__;
//static const char compileTime[] = __TIME__;
#endif /* NOT FRAMEWORK_VERSION */
const char *liquid_CO2_H2O_ver(void) { return "$Id: liquid_CO2_H2O.c,v 1.42 2009/05/14 04:24:00 ghiorso Exp $"; }


/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Routines to compute liquid solution properties
**      (file: LIQUID_CO2_H2O.C)
**
**--
*/

#ifdef DEBUG
#undef DEBUG
#endif

#define USE_KRESS_CARMICHAEL_FO2

#include "melts_gsl.h"

#include "silmin.h"
#include "mthread.h"
#ifdef USESEH
#include <windows.h>
void raise_sigabrt(DWORD dwType);
extern int doInterrupt;
#else
#include <signal.h>
#endif

#include "param_struct_data_CO2_H2O.h"

#define SQUARE(x) ((x)*(x))
#define CUBE(x)   ((x)*(x)*(x))

// from olivine.c
#define USE_SVD 0
#define MAX_ITER 200  /* Maximum number of iterations allowed in order */

/*
 *=============================================================================
 * Private functions and globals:
*/

#define NA 19                    /* Number of liquid components                   */
#define NS  1                    /* Number of ordering parameters for liquid species      */
#define NY  0                    /* Number of ordering parameters for coordination states */

#define NT (NS+NY)               /* Number of ordering parameters                         */
#define NR (NA-1)                /* Number of independent mole fraction variables         */
#define NW ((NA+NS)*(NA+NS-1)/2) /* Number of regular solution interaction parameters     */
#define NV (NR+NS)               /* Number of independent variables in the model          */
#define NP (NA+NS+NW)            /* Number of model parameters                    */

static int nH2O, nCO2;

/* The statics from here to ... */

static int convergedInOrder;

static double *xSpecies;         /* Mole fractions of endmember species                                                  */
static double **dxSpeciesdr;     /* d(xSpecies)/dr                                                                       */
static double **dxSpeciesds;     /* d(xSpecies)/ds                                                                       */
static double ***d2xSpeciesdrds; /* d2(xSpecies)/drds                                                                    */
static double nSpecies;          /* Total moles of all solution species relative to 1 mole of basis species              */
static double *dnSpeciesds;      /* d(nSpecies)/ds                                                                       */
static double **d2nSpeciesds2;   /* d2(nSpecies)/ds2                                                                     */
static double ***d3nSpeciesds3;  /* d3(nSpecies)/ds3                                                                     */

static MTHREAD_MUTEX_T global_data_mutex = MTHREAD_MUTEX_INITIALIZER;

/* ... here are dealt with by creating a mutex lock on the code that requires access to these quantities */

/***********************************************/
/* Statics for class initialization structures */
/***********************************************/

static MTHREAD_ONCE_T initThreadBlock = MTHREAD_ONCE_INIT;

static void initializeLiquid(void);

static void threadInit(void) {
    initializeLiquid();
}

/* The statics from here ... */

static int NE;   /* Number of liquid endmembers (species) */

static double Gconst,       *gr,       *gs,       **grr,        **grs,       **gss;
static double Hconst,       *hr,       *hs,       **hrr,        **hrs,       **hss;
static double Sconst,       *sr,       *ss,       **srr,        **srs,       **sss;
static double Vconst,       *vr,       *vs,       **vrr,        **vrs,       **vss;
static double CPconst,      *cpr,      *cps;
static double DCPDTconst,   *dcpdtr,   *dcpdts;
static double DVDTconst,    *dvdtr,    *dvdts;
static double DVDPconst,    *dvdpr,    *dvdps;
static double D2VDT2const,  *d2vdt2r,  *d2vdt2s;
static double D2VDTDPconst, *d2vdtdpr, *d2vdtdps;
static double D2VDP2const,  *d2vdp2r,  *d2vdp2s;

static double **taylorCoeff;  /* Taylor Expansion coefficients: [endmember species | W(i,j)][g0, gr, gs, grr, grs, gss] */
static double **rsEndmembers; /* r and s coefficients for species endmembers                                            */

/* ... to here are set in initializeLiquid() */

/* This is adapted from olivine.c */
static MTHREAD_ONCE_T initThreadOBlock = MTHREAD_ONCE_INIT;

static MTHREAD_KEY_T tOldKey;
static MTHREAD_KEY_T pOldKey;
static MTHREAD_KEY_T rOldKey;
static MTHREAD_KEY_T sOldKey;
static MTHREAD_KEY_T ptToD2gds2Key;
static MTHREAD_KEY_T d2gds2Key;
static MTHREAD_KEY_T ptToVKey;
static MTHREAD_KEY_T ptToSKey;
static MTHREAD_KEY_T indexD2gds2Key;

static void freeNSvector(void *NSarray) {
    gsl_vector_free((gsl_vector *) NSarray);
}

static void freeNSmatrix(void *NSarray) {
    gsl_matrix_free((gsl_matrix *) NSarray);
}

static void freeIndexD2gds2(void *indexD2gds2) {
    gsl_permutation_free((gsl_permutation *) indexD2gds2);
}

static void threadOInit(void) {
    MTHREAD_KEY_CREATE(&tOldKey,       free);
    MTHREAD_KEY_CREATE(&pOldKey,       free);
    MTHREAD_KEY_CREATE(&rOldKey,       freeNSvector);
    MTHREAD_KEY_CREATE(&sOldKey,       freeNSvector);
    MTHREAD_KEY_CREATE(&ptToD2gds2Key, freeNSmatrix);
    MTHREAD_KEY_CREATE(&d2gds2Key,     free);
    MTHREAD_KEY_CREATE(&ptToVKey,      freeNSmatrix);
    MTHREAD_KEY_CREATE(&ptToSKey,      freeNSvector);
    MTHREAD_KEY_CREATE(&indexD2gds2Key, freeIndexD2gds2);
}

static double getTOld() {
    double *tOldPt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    tOldPt = (double *) MTHREAD_GETSPECIFIC(tOldKey);
    if (tOldPt == NULL) {
        tOldPt  = (double *) malloc(sizeof(double));
        *tOldPt = -9999.0;
        MTHREAD_SETSPECIFIC(tOldKey, (void *) tOldPt);
    }
    return *tOldPt;
}

static void setTOld(double tOld) {
    double *tOldPt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    tOldPt = (double *) MTHREAD_GETSPECIFIC(tOldKey);
    if (tOldPt == NULL) {
        tOldPt  = (double *) malloc(sizeof(double));
        *tOldPt = -9999.0;
        MTHREAD_SETSPECIFIC(tOldKey, (void *) tOldPt);
    }
    *tOldPt = tOld;
}

static double getPOld() {
    double *pOldPt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    pOldPt = (double *) MTHREAD_GETSPECIFIC(pOldKey);
    if (pOldPt == NULL) {
        pOldPt  = (double *) malloc(sizeof(double));
        *pOldPt = -9999.0;
        MTHREAD_SETSPECIFIC(pOldKey, (void *) pOldPt);
    }
    return *pOldPt;
}

static void setPOld(double pOld) {
    double *pOldPt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    pOldPt = (double *) MTHREAD_GETSPECIFIC(pOldKey);
    if (pOldPt == NULL) {
        pOldPt  = (double *) malloc(sizeof(double));
        *pOldPt = -9999.0;
        MTHREAD_SETSPECIFIC(pOldKey, (void *) pOldPt);
    }
    *pOldPt = pOld;
}

static double *getROld() {
    gsl_vector *rOldPt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    rOldPt = (gsl_vector *) MTHREAD_GETSPECIFIC(rOldKey);
    if (rOldPt == NULL) {
        rOldPt = gsl_vector_alloc((size_t) NR);
        gsl_vector_set_all(rOldPt, -9999.0);
        MTHREAD_SETSPECIFIC(rOldKey, (void *) rOldPt);
    }
    return rOldPt->data;
}

static double *getSOld() {
    gsl_vector *sOldPt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    sOldPt = (gsl_vector *) MTHREAD_GETSPECIFIC(sOldKey);
    if (sOldPt == NULL) {
        sOldPt = gsl_vector_alloc((size_t) NT);
        gsl_vector_set_all(sOldPt, 0.0);
        MTHREAD_SETSPECIFIC(sOldKey, (void *) sOldPt);
    }
    return sOldPt->data;
}

static gsl_matrix *getPtToD2gds2() {
    gsl_matrix *ptToD2gds2Pt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    ptToD2gds2Pt = (gsl_matrix *) MTHREAD_GETSPECIFIC(ptToD2gds2Key);
    if (ptToD2gds2Pt == NULL) {
        ptToD2gds2Pt  = gsl_matrix_alloc((size_t) NT, (size_t) NT);
        gsl_matrix_set_zero(ptToD2gds2Pt);
        MTHREAD_SETSPECIFIC(ptToD2gds2Key, (void *) ptToD2gds2Pt);
    }
    return ptToD2gds2Pt;
}

static double **getD2gds2() {
    double **d2gds2Pt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    d2gds2Pt = (double **) MTHREAD_GETSPECIFIC(d2gds2Key);
    if (d2gds2Pt == NULL) {
        int i;
        gsl_matrix *ptToD2gds2Pt = getPtToD2gds2();
        d2gds2Pt  = (double **) malloc((size_t) NT*sizeof(double *));
        for (i=0; i<NT; i++) d2gds2Pt[i] = &(ptToD2gds2Pt->data[i*NT]);
        MTHREAD_SETSPECIFIC(d2gds2Key, (void *) d2gds2Pt);
    }
    return d2gds2Pt;
}

static gsl_matrix *getPtToV() {
    gsl_matrix *ptToVPt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    ptToVPt = (gsl_matrix *) MTHREAD_GETSPECIFIC(ptToVKey);
    if (ptToVPt == NULL) {
        ptToVPt  = gsl_matrix_alloc((size_t) NT, (size_t) NT);
        gsl_matrix_set_zero(ptToVPt);
        MTHREAD_SETSPECIFIC(ptToVKey, (void *) ptToVPt);
    }
    return ptToVPt;
}

static gsl_vector *getPtToS() {
    gsl_vector *ptToSPt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    ptToSPt = (gsl_vector *) MTHREAD_GETSPECIFIC(ptToSKey);
    if (ptToSPt == NULL) {
        ptToSPt  = gsl_vector_alloc((size_t) NT);
        gsl_vector_set_zero(ptToSPt);
        MTHREAD_SETSPECIFIC(ptToSKey, (void *) ptToSPt);
    }
    return ptToSPt;
}

static gsl_permutation *getIndexD2gds2() {
    gsl_permutation *indexD2gds2Pt;
    MTHREAD_ONCE(&initThreadOBlock, threadOInit);

    indexD2gds2Pt = (gsl_permutation *) MTHREAD_GETSPECIFIC(indexD2gds2Key);
    if (indexD2gds2Pt == NULL) {
        indexD2gds2Pt = gsl_permutation_alloc((size_t) NT);
        gsl_permutation_init(indexD2gds2Pt);
        MTHREAD_SETSPECIFIC(indexD2gds2Key, (void *) indexD2gds2Pt);
    }
    return indexD2gds2Pt;
}

static double **matrix_alloc(int n1, int n2) {
    int i;
    double *m0 = (double *) malloc((size_t) n1*n2*sizeof(double));
    double **m = (double **) malloc((size_t) n1*sizeof(double *));
    for (i=0; i<n1; i++) m[i] = &m0[i*n2];
    return m;
}

static void matrix_free(double **m, int n1, int n2) {
    free(m[0]);
    free(m);
}

static double *vector_alloc(int n) {
    double *v = (double *) malloc((size_t) n*sizeof(double));
    return v;
}

static void vector_free(double *v, int n) {
    free(v);
}

static void initializeTaylor(gsl_vector *array) {
    gsl_vector_set_zero(array);
}

static void loadTaylor(gsl_vector *vTaylor, double coeff, double *p) {
    int i, j, n;
    double array[NP];
    gsl_vector_view temp;

    array[0] = coeff;
    for (i=0, n=1; i<NV; i++, n++) array[n] = (p[i] != 0.0) ? coeff*p[i] : 0.0;
    for (i=0; i<NV; i++) for (j=i; j<NV; j++, n++)
        array[n] =  (p[i] != 0.0 && p[j] != 0.0) ? coeff*p[i]*p[j] : 0.0;

    temp = gsl_vector_view_array(array, (size_t) NP);
    gsl_vector_add(vTaylor, &temp.vector);

}

static void initializeLiquid(void) {
    int i, j, k, m, n;
    static int contentsOK = FALSE;
    double *temp, *coeff;
    gsl_matrix *mTaylor, *bTaylor;

#ifdef DEBUG
    printf("Entering function initializeLiquid...\n");
#endif

    /* Initialize global constants */
    NE = nls; /* Number of liquid endmembers (species) */
    if ( (NE == 0) || (NA != nlc) || (NS != (nls-nlc)) ) {
        if (NE == 0)         printf("Error in initializeLiquid() [liquid.c].  Number of species is zero.\n");
        if (NA != nlc)       printf("Error in initializeLiquid() [liquid.c].  Constant NA is not equal to nlc (%d).\n", nlc);
        if (NS != (nls-nlc)) printf("Error in initializeLiquid() [liquid.c].  Constant NS is not equal to nls-nlc (%d).\n", nls-nlc);
        exit(0);
    }

    nH2O = -1;
    for (i=0; i<NE; i++) if ((strcmp(liquid[i].label, "H2O") == 0) || (strcmp(liquid[i].label, "h2o") == 0)) { nH2O = i; break; }
    nCO2 = -1;
    for (i=0; i<NE; i++) if ((strcmp(liquid[i].label, "CO2") == 0) || (strcmp(liquid[i].label, "co2") == 0)) { nCO2 = i; break; }

    /* Static global storage for endmember species mole fractions */
    xSpecies       = vector_alloc(NE);
    dxSpeciesdr    = matrix_alloc(NE, NR);
    dxSpeciesds    = matrix_alloc(NE, NS);
    d2xSpeciesdrds = (double ***) malloc((size_t) NE*sizeof(double **));
                   for (i=0; i<NE; i++) d2xSpeciesdrds[i] = matrix_alloc(NR, NS);
    dnSpeciesds    = vector_alloc(NS);
    d2nSpeciesds2  = matrix_alloc(NS, NS);
    d3nSpeciesds3  = (double ***) malloc((size_t) NS*sizeof(double **));
                   for (i=0; i<NS; i++) d3nSpeciesds3[i] = matrix_alloc(NS, NS);

    /* Allocate static storage for the derivatives of g */

    gr = vector_alloc(NR);
    gs = vector_alloc(NS);
    grr = matrix_alloc(NR, NR);
    grs = matrix_alloc(NR, NS);
    gss = matrix_alloc(NS, NS);
    hr = vector_alloc(NR);
    hs = vector_alloc(NS);
    hrr = matrix_alloc(NR, NR);
    hrs = matrix_alloc(NR, NS);
    hss = matrix_alloc(NS, NS);
    sr = vector_alloc(NR);
    ss = vector_alloc(NS);
    srr = matrix_alloc(NR, NR);
    srs = matrix_alloc(NR, NS);
    sss = matrix_alloc(NS, NS);
    vr = vector_alloc(NR);
    vs = vector_alloc(NS);
    vrr = matrix_alloc(NR, NR);
    vrs = matrix_alloc(NR, NS);
    vss = matrix_alloc(NS, NS);
    cpr = vector_alloc(NR);
    cps = vector_alloc(NS);
    dcpdtr = vector_alloc(NR);
    dcpdts = vector_alloc(NS);
    dvdtr = vector_alloc(NR);
    dvdts = vector_alloc(NS);
    dvdpr = vector_alloc(NR);
    dvdps = vector_alloc(NS);
    d2vdt2r = vector_alloc(NR);
    d2vdt2s = vector_alloc(NS);
    d2vdtdpr = vector_alloc(NR);
    d2vdtdps = vector_alloc(NS);
    d2vdp2r = vector_alloc(NR);
    d2vdp2s = vector_alloc(NS);

#ifdef DEBUG
    printf("...Constructing r,s --> Endmember species matrix.\n");
#endif
    /* r and s coefficients for endmember species */
    rsEndmembers = matrix_alloc (NE, NV);
    coeff        = vector_alloc (NE);
    for (i=0; i<NE; i++) for (j=0; j<NV; j++) rsEndmembers[i][j] = 0.0;
    for (i=0; i<NA; i++) coeff[i] = 1.0;
    /* --> basis species     */
    for (i=1; i<NA; i++) rsEndmembers[i][i-1] = 1.0;
    /* --> dependent species */
    for (i=NA; i<NE; i++) {
        for (k=0, coeff[i]=0.0; k<nc; k++) coeff[i] += (bulkSystem[k].oxToLiq)[0]*(liquid[i].liqToOx)[k];
        for (j=1; j<NA; j++) for (k=0, rsEndmembers[i][j-1]=0.0; k<nc; k++)
            rsEndmembers[i][j-1] += (bulkSystem[k].oxToLiq)[j]*(liquid[i].liqToOx)[k];
        for (j=0; j<NR; j++) coeff[i] += rsEndmembers[i][j];
        if (coeff[i] != 0.0) for (j=0; j<NR; j++) rsEndmembers[i][j] /= coeff[i];
        rsEndmembers[i][i-NA+NR] = 1.0;  /* ordering parameter */
    }
#ifdef DEBUG
    printf("...rsEndmembers matrix:\n");
    for (i=0; i<NE; i++) {
        printf("coeff[%-15.15s] = %10.4f\n", liquid[i].label, coeff[i]);
        for (j=0; j<NV; j++) printf("%4.1f ", rsEndmembers[i][j]);
        printf("\n");
    }
#endif

#ifdef DEBUG
    printf("...Zeroing dxdr, dxds, and d2xdrds matrices.\n");
#endif
    for (i=0; i<NE; i++) {
        for (j=0; j<NR; j++) {
            dxSpeciesdr[i][j] = 0.0;
            for (k=0; k<NS; k++) d2xSpeciesdrds[i][j][k] = 0.0;
        }
        for (j=0; j<NS; j++) dxSpeciesds[i][j] = 0.0;
    }

    /********************************************
    * Create the liquid model Taylor expansion *
    ********************************************/

    if (!contentsOK) {

        gsl_vector *S;
        gsl_matrix *V;
        gsl_matrix_view A;
        gsl_vector_view x, b;


#ifdef DEBUG
        printf("...Constructing Taylor expansion coefficients of liquid model.\n");
#endif
        mTaylor = gsl_matrix_alloc (NP, NP);

        /* First-order Taylor terms */
        for (i=0; i<NE; i++) {
            x = gsl_matrix_row(mTaylor, i);
            initializeTaylor(&x.vector);
            loadTaylor (&x.vector, coeff[i], &rsEndmembers[i][0]);
        }

        /* Second-order Taylor terms: binary interaction parameters
        join A-B ==>  4 ( G(A/2+B/2) - G(A)/2 - G(B)/2 )        */
        m = NE; n = NP; //1 + NV + NV*NV;
        temp = vector_alloc(NE);
        for (i=0; i<NE; i++) {
            for (j=i+1; j<NE; j++) {
                x = gsl_matrix_row(mTaylor, m);
                initializeTaylor(&x.vector);
                for (k=0; k<NV; k++) temp[k] = (rsEndmembers[i][k] + rsEndmembers[j][k])/2.0;
                loadTaylor (&x.vector, (double)  4.0, temp);
                loadTaylor (&x.vector, (double) -2.0, &rsEndmembers[i][0]);
                loadTaylor (&x.vector, (double) -2.0, &rsEndmembers[j][0]);
                m++;
            }
        }
        vector_free(temp,  NE);
        vector_free(coeff, NE);

        // SVD not implemented for M < N in GSL
        if (m < n) gsl_matrix_transpose(mTaylor);

        /* Fill in uninitialized entries */
        bTaylor    = gsl_matrix_alloc((size_t) NP, (size_t) NP);

        taylorCoeff = matrix_alloc(NP, NP);
        S = gsl_vector_alloc((size_t) NP);
        V = gsl_matrix_alloc((size_t) NP, (size_t) NP);
        A = gsl_matrix_submatrix(mTaylor, (size_t) 0, (size_t) 0, (size_t) NP, (size_t) NP);
        x = gsl_vector_view_array(taylorCoeff[0], (size_t) NP);

#ifdef DEBUG
        printf("...Performing singular value analysis of coefficient matrix...\n");
#endif
        gsl_linalg_SV_decomp(&A.matrix, V, S, &x.vector); // first row of taylorCoeff used as working space
        for (i=0, j=0; i<NP; i++) if (fabs(gsl_vector_get(S, i)) < DBL_EPSILON) { gsl_vector_set(S, i, 0.0); j++; }
#ifdef DEBUG
        if (NP-j < n) printf("...Problem is rank deficient! rank = %d\n", NP-j);
        printf("...Performing back-substitution phase...\n");
#endif

        gsl_matrix_set_identity(bTaylor);
        for (i=0; i<NP; i++) {
            b = gsl_matrix_row(bTaylor, i);
            x = gsl_vector_view_array(taylorCoeff[i], (size_t) NP);
            // SVD not implemented for M < N in GSL
            if (m >= n) gsl_linalg_SV_solve(&A.matrix, V, S, &b.vector, &x.vector);
            else gsl_linalg_SV_solve(V, &A.matrix, S, &b.vector, &x.vector);
        }

        for (i=0; i<NP; i++) for (j=0; j<NP; j++)
            if (fabs(taylorCoeff[i][j]) < sqrt(DBL_EPSILON)) taylorCoeff[i][j] = 0.0;

        gsl_matrix_free(V); gsl_vector_free(S);
        gsl_matrix_free(bTaylor); gsl_matrix_free(mTaylor);
        contentsOK = TRUE;

    }

#ifdef DEBUG
    printf("...Exiting function initializeLiquid.\n");
#endif
}

static int rANDsTOx (double r[NR], double s[NT]) {
    static double tolerance;
    double rSum, coeff, dcoeffds[NS], denom, dDenomds[NS];
    int i, j, k, okay = TRUE;
    /* const double y = 0.3; */ /* Fe2SiO(4+2*0.3) or Fe2AlO(3.5+2*0.3) */

    for (i=0, rSum=0.0; i<NR; i++) rSum += r[i];

    coeff = 1.0;
        for (i=0; i<NS; i++) dcoeffds[i] = 0.0;

    /* xSpecies */
    xSpecies[ 0] = 1.0 - rSum*coeff;                   /* SiO2  */
    for (i=0; i<NR; i++) xSpecies[ i+1] = r[i]*coeff;  /* basis */
    for (i=0; i<NS; i++) xSpecies[NA+i] = s[i];        /* depen */

    xSpecies[ 0] += s[0]; // special case SiO2
    xSpecies[10] -= s[0]; // special case CaSiO3
    xSpecies[14] -= s[0]; // special case CO2

    /* Catch bad input data */
    for (i=0;  i<NE; i++) okay &= (xSpecies[i] >= 0.0);
    if (!okay) return okay;

    for (i=NS; i<NT; i++) okay &= ((s[i] > (10.0*DBL_EPSILON)) && (s[i] < (1.0-10.0*DBL_EPSILON)));
    if (okay && (NT > NS)) {
        double yIV = 1.0;
        for (i=NS; i<NT; i++) yIV -= s[i];
        okay &= ((yIV > (10.0*DBL_EPSILON)) && (yIV < (1.0-10.0*DBL_EPSILON)));
    }
    if (!okay) return okay;

    /* Correct roundoff problems - removed check on 4/10/02 when MgO species was included */
    if (tolerance == 0.0) tolerance = pow(DBL_EPSILON, (double) (2.0/3.0));
    /*  for (i=0; i<(NA+NS); i++) if (fabs(xSpecies[i]) < tolerance) xSpecies[i] = 0.0; */

    /* d xSpecies / dr */
    for (i=0; i<NR; i++) {
        dxSpeciesdr[  0][i] = - coeff;  /* SiO2  */
        dxSpeciesdr[i+1][i] =   coeff;  /* other */
    }

    /* d xSpecies / ds */
    for (i=0; i<NS; i++) {
        dxSpeciesds[   0][i] = - rSum*dcoeffds[i]; /* SiO2  */
        for (j=0; j<NR; j++) dxSpeciesds[ j+1][i] =   r[j]*dcoeffds[i]; /* basis */
        dxSpeciesds[NA+i][i] = 1.0;                /* depen */
    }

    dxSpeciesds[ 0][0] += 1.0; // special case SiO2
    dxSpeciesds[10][0] -= 1.0; // special case CaCO3
    dxSpeciesds[14][0] -= 1.0; // special case CO2

    /* d2 xSpecies / dr ds */
    for (i=0; i<NR; i++) {
        for (j=0; j<NS; j++) {
            d2xSpeciesdrds[  0][i][j] = -dcoeffds[j];  /* SiO2  */
            d2xSpeciesdrds[i+1][i][j] =  dcoeffds[j];  /* other */
        }
    }

    /* Total moles of species relative to 1 mole of basis components */
    denom = 12.0;                                    /* Special case */
    dDenomds[0] = 0.0;  // CaCO3

    nSpecies = 12.0/denom;
    for (i=0; i<NS; i++) {
        dnSpeciesds[i] = -12.0*dDenomds[i]/(denom*denom);
        for (j=0; j<NS; j++) {
            d2nSpeciesds2[i][j] = 24.0*dDenomds[i]*dDenomds[j]/(denom*denom*denom);
            for (k=0; k<NS; k++) d3nSpeciesds3[i][j][k] = -72.0*dDenomds[i]*dDenomds[j]*dDenomds[k]/(denom*denom*denom*denom);
        }
    }

    return okay;
}

#define WH(k) (meltsAndCO2_H2OModelParameters[k].enthalpy)
#define WS(k) (meltsAndCO2_H2OModelParameters[k].entropy)
#define WV(k) (meltsAndCO2_H2OModelParameters[k].volume)
#define W(k)  (meltsAndCO2_H2OModelParameters[k].enthalpy - t*meltsAndCO2_H2OModelParameters[k].entropy  + (p-1.0)*meltsAndCO2_H2OModelParameters[k].volume)

#define G(i)       ((liquid[i].cur).g + meltsAndCO2_H2OModelParameters[NW+i].enthalpy - t*meltsAndCO2_H2OModelParameters[NW+i].entropy + (p-1.0)*meltsAndCO2_H2OModelParameters[NW+i].volume)
#define H(i)       ((liquid[i].cur).h + meltsAndCO2_H2OModelParameters[NW+i].enthalpy)
#define S(i)       ((liquid[i].cur).s + meltsAndCO2_H2OModelParameters[NW+i].entropy)
#define V(i)       ((liquid[i].cur).v + meltsAndCO2_H2OModelParameters[NW+i].volume)
#define CP(i)      ((liquid[i].cur).cp)
#define DCPDT(i)   ((liquid[i].cur).dcpdt)
#define DVDT(i)    ((liquid[i].cur).dvdt)
#define DVDP(i)    ((liquid[i].cur).dvdp)
#define D2VDT2(i)  ((liquid[i].cur).d2vdt2)
#define D2VDTDP(i) ((liquid[i].cur).d2vdtdp)
#define D2VDP2(i)  ((liquid[i].cur).d2vdp2)

static void loadTaylorCoefficients(double t, double p)
{
    int i, j, k, l, m, n;

#ifdef DEBUG
    printf("Call to loadTaylorCoefficients ...\n");
#endif

    for (i=0; i<NE; i++)  gibbs(t,   p, (char *) liquid[i].label, &(liquid[i].ref), &(liquid[i].liq), &(liquid[i].fus), &(liquid[i].cur));

    Gconst       = 0.0;  Hconst       = 0.0;  Sconst       = 0.0;  Vconst       = 0.0;
    CPconst      = 0.0;  DCPDTconst   = 0.0;  DVDTconst    = 0.0;  DVDPconst    = 0.0;
    D2VDT2const  = 0.0;  D2VDTDPconst = 0.0;  D2VDP2const  = 0.0;

    memset(gr      , '\0', (size_t) NR*sizeof(double));
    memset(hr      , '\0', (size_t) NR*sizeof(double));
    memset(sr      , '\0', (size_t) NR*sizeof(double));
    memset(vr      , '\0', (size_t) NR*sizeof(double));
    memset(cpr     , '\0', (size_t) NR*sizeof(double));
    memset(dcpdtr  , '\0', (size_t) NR*sizeof(double));
    memset(dvdtr   , '\0', (size_t) NR*sizeof(double));
    memset(dvdpr   , '\0', (size_t) NR*sizeof(double));
    memset(d2vdt2r , '\0', (size_t) NR*sizeof(double));
    memset(d2vdtdpr, '\0', (size_t) NR*sizeof(double));
    memset(d2vdp2r , '\0', (size_t) NR*sizeof(double));
    for (j=0; j<NR; j++) {
        memset(grr[j], '\0', (size_t) NR*sizeof(double));
        memset(hrr[j], '\0', (size_t) NR*sizeof(double));
        memset(srr[j], '\0', (size_t) NR*sizeof(double));
        memset(vrr[j], '\0', (size_t) NR*sizeof(double));
        memset(grs[j], '\0', (size_t) NS*sizeof(double));
        memset(hrs[j], '\0', (size_t) NS*sizeof(double));
        memset(srs[j], '\0', (size_t) NS*sizeof(double));
        memset(vrs[j], '\0', (size_t) NS*sizeof(double));
    }
    memset(gs      , '\0', (size_t) NS*sizeof(double));
    memset(hs      , '\0', (size_t) NS*sizeof(double));
    memset(ss      , '\0', (size_t) NS*sizeof(double));
    memset(vs      , '\0', (size_t) NS*sizeof(double));
    memset(cps     , '\0', (size_t) NS*sizeof(double));
    memset(dcpdts  , '\0', (size_t) NS*sizeof(double));
    memset(dvdts   , '\0', (size_t) NS*sizeof(double));
    memset(dvdps   , '\0', (size_t) NS*sizeof(double));
    memset(d2vdt2s , '\0', (size_t) NS*sizeof(double));
    memset(d2vdtdps, '\0', (size_t) NS*sizeof(double));
    memset(d2vdp2s , '\0', (size_t) NS*sizeof(double));
    for (j=0; j<NS; j++) {
        memset(gss[j], '\0', (size_t) NS*sizeof(double));
        memset(hss[j], '\0', (size_t) NS*sizeof(double));
        memset(sss[j], '\0', (size_t) NS*sizeof(double));
        memset(vss[j], '\0', (size_t) NS*sizeof(double));
    }

    for (i=0; i<NE; i++) {
        Gconst       += taylorCoeff[i][0]*G(i);
        Hconst       += taylorCoeff[i][0]*H(i);
        Sconst       += taylorCoeff[i][0]*S(i);
        Vconst       += taylorCoeff[i][0]*V(i);
        CPconst      += taylorCoeff[i][0]*CP(i);
        DCPDTconst   += taylorCoeff[i][0]*DCPDT(i);
        DVDTconst    += taylorCoeff[i][0]*DVDT(i);
        DVDPconst    += taylorCoeff[i][0]*DVDP(i);
        D2VDT2const  += taylorCoeff[i][0]*D2VDT2(i);
        D2VDTDPconst += taylorCoeff[i][0]*D2VDTDP(i);
        D2VDP2const  += taylorCoeff[i][0]*D2VDP2(i);
        for (j=0; j<NR; j++) {
            gr[j]       += taylorCoeff[i][1+j]*G(i);
            hr[j]       += taylorCoeff[i][1+j]*H(i);
            sr[j]       += taylorCoeff[i][1+j]*S(i);
            vr[j]       += taylorCoeff[i][1+j]*V(i);
            cpr[j]      += taylorCoeff[i][1+j]*CP(i);
            dcpdtr[j]   += taylorCoeff[i][1+j]*DCPDT(i);
            dvdtr[j]    += taylorCoeff[i][1+j]*DVDT(i);
            dvdpr[j]    += taylorCoeff[i][1+j]*DVDP(i);
            d2vdt2r[j]  += taylorCoeff[i][1+j]*D2VDT2(i);
            d2vdtdpr[j] += taylorCoeff[i][1+j]*D2VDTDP(i);
            d2vdp2r[j]  += taylorCoeff[i][1+j]*D2VDP2(i);
        }
        for (j=0; j<NS; j++) {
            gs[j]       += taylorCoeff[i][1+NR+j]*G(i);
            hs[j]       += taylorCoeff[i][1+NR+j]*H(i);
            ss[j]       += taylorCoeff[i][1+NR+j]*S(i);
            vs[j]       += taylorCoeff[i][1+NR+j]*V(i);
            cps[j]      += taylorCoeff[i][1+NR+j]*CP(i);
            dcpdts[j]   += taylorCoeff[i][1+NR+j]*DCPDT(i);
            dvdts[j]    += taylorCoeff[i][1+NR+j]*DVDT(i);
            dvdps[j]    += taylorCoeff[i][1+NR+j]*DVDP(i);
            d2vdt2s[j]  += taylorCoeff[i][1+NR+j]*D2VDT2(i);
            d2vdtdps[j] += taylorCoeff[i][1+NR+j]*D2VDTDP(i);
            d2vdp2s[j]  += taylorCoeff[i][1+NR+j]*D2VDP2(i);
        }
    }
    /* Code below is optimized for speed of execution ... */
    for (i=0, n=0; i<NE; i++) for (l=i+1; l<NE; l++, n++) {
        register double w  = W(n);
        if (w != 0.0) {
            Gconst += taylorCoeff[n+NE][0]*w;
            for (j=0, m=0; j<NR; j++) {
                gr[j] += taylorCoeff[n+NE][1+j]*w;
                for (k=j; k<NR; k++, m++) grr[j][k] += taylorCoeff[n+NE][1+NR+NS+m]*w;
                for (k=0; k<NS; k++, m++) grs[j][k] += taylorCoeff[n+NE][1+NR+NS+m]*w;
            }
            for (j=0; j<NS; j++) {
                gs[j] += taylorCoeff[n+NE][1+NR+j]*w;
                for (k=j; k<NS; k++, m++) if (taylorCoeff[n+NE][1+NR+NS+m] != 0.0) gss[j][k] += taylorCoeff[n+NE][1+NR+NS+m]*w;
            }
        }
    }
    for (i=0, n=0; i<NE; i++) for (l=i+1; l<NE; l++, n++) {
        register double wh = WH(n);
        if (wh != 0.0) {
        Hconst += taylorCoeff[n+NE][0]*wh;
        for (j=0, m=0; j<NR; j++) {
            hr[j] += taylorCoeff[n+NE][1+j]*wh;
            for (k=j; k<NR; k++, m++) hrr[j][k] += taylorCoeff[n+NE][1+NR+NS+m]*wh;
            for (k=0; k<NS; k++, m++) hrs[j][k] += taylorCoeff[n+NE][1+NR+NS+m]*wh;
        }
        for (j=0; j<NS; j++) {
            hs[j] += taylorCoeff[n+NE][1+NR+j]*wh;
            for (k=j; k<NS; k++, m++) if (taylorCoeff[n+NE][1+NR+NS+m] != 0.0) hss[j][k] += taylorCoeff[n+NE][1+NR+NS+m]*wh;
        }
        }
    }
    for (i=0, n=0; i<NE; i++) for (l=i+1; l<NE; l++, n++) {
        register double ws = WS(n);
        if (ws != 0.0) {
        Sconst += taylorCoeff[n+NE][0]*ws;
        for (j=0, m=0; j<NR; j++) {
            sr[j] += taylorCoeff[n+NE][1+j]*ws;
            for (k=j; k<NR; k++, m++) srr[j][k] += taylorCoeff[n+NE][1+NR+NS+m]*ws;
            for (k=0; k<NS; k++, m++) srs[j][k] += taylorCoeff[n+NE][1+NR+NS+m]*ws;
        }
        for (j=0; j<NS; j++) {
            ss[j] += taylorCoeff[n+NE][1+NR+j]*ws;
            for (k=j; k<NS; k++, m++) if (taylorCoeff[n+NE][1+NR+NS+m] != 0.0) sss[j][k] += taylorCoeff[n+NE][1+NR+NS+m]*ws;
        }
        }
    }
    for (i=0, n=0; i<NE; i++) for (l=i+1; l<NE; l++, n++) {
        register double wv = WV(n);
        if (wv != 0.0) {
        Vconst += taylorCoeff[n+NE][0]*wv;
        for (j=0, m=0; j<NR; j++) {
            vr[j] += taylorCoeff[n+NE][1+j]*wv;
            for (k=j; k<NR; k++, m++) vrr[j][k] += taylorCoeff[n+NE][1+NR+NS+m]*wv;
            for (k=0; k<NS; k++, m++) vrs[j][k] += taylorCoeff[n+NE][1+NR+NS+m]*wv;
        }
        for (j=0; j<NS; j++) {
            vs[j] += taylorCoeff[n+NE][1+NR+j]*wv;
            for (k=j; k<NS; k++, m++) if (taylorCoeff[n+NE][1+NR+NS+m] != 0.0) vss[j][k] += taylorCoeff[n+NE][1+NR+NS+m]*wv;
        }
        }
    }
    /* ... end speed optimization. */

    for (j=0; j<NR; j++) for (k=j+1; k<NR; k++) {
        grr[k][j]       = grr[j][k];
        hrr[k][j]       = hrr[j][k];
        srr[k][j]       = srr[j][k];
        vrr[k][j]       = vrr[j][k];
    }

    for (j=0; j<NS; j++) for (k=j+1; k<NS; k++) {
        gss[k][j]       = gss[j][k];
        hss[k][j]       = hss[j][k];
        sss[k][j]       = sss[j][k];
        vss[k][j]       = vss[j][k];
    }

}

static double fillG (double r[NR], double s[NT], double t, double p) {
    double result, config;
    int i, j;

    /* Taylor expansion and standard state terms */
    result = Gconst;
    for (i=0; i<NR; i++) {
        result += gr[i]*r[i];
        for (j=i; j<NR; j++) result += grr[i][j]*r[i]*r[j];
        for (j=0; j<NS; j++) result += grs[i][j]*r[i]*s[j];
    }
    for (i=0; i<NS; i++) {
        result += gs[i]*s[i];
        for (j=i; j<NS; j++) result += gss[i][j]*s[i]*s[j];
    }

    /* Configurational entropy terms */
    for (i=0, config=0.0; i<NE; i++) if (xSpecies[i] > 0.0) config += xSpecies[i]*log(xSpecies[i]);
    if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0) config += xSpecies[nH2O]*log(xSpecies[nH2O]) + (1.0-xSpecies[nH2O])*log(1.0-xSpecies[nH2O]);
    result += R*t*nSpecies*config;

    return result;
}

static double fillS (double r[NR], double s[NT], double t, double p) {
    double result, config;
    int i, j;

    /* Taylor expansion and standard state terms */
    result = Sconst;
    for (i=0; i<NR; i++) {
        result += sr[i]*r[i];
        for (j=i; j<NR; j++) result += srr[i][j]*r[i]*r[j];
        for (j=0; j<NS; j++) result += srs[i][j]*r[i]*s[j];
    }
    for (i=0; i<NS; i++) {
        result += ss[i]*s[i];
        for (j=i; j<NS; j++) result += sss[i][j]*s[i]*s[j];
    }

    /* Configurational entropy terms */
    for (i=0, config=0.0; i<NE; i++) if (xSpecies[i] > 0.0) config += xSpecies[i]*log(xSpecies[i]);
    if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0) config += xSpecies[nH2O]*log(xSpecies[nH2O]) + (1.0-xSpecies[nH2O])*log(1.0-xSpecies[nH2O]);
    result -= R*nSpecies*config;

    return result;
}

static double fillV (double r[NR], double s[NT], double t, double p) {
    double result;
    int i, j;

    /* Taylor expansion and standard state terms */
    result = Vconst;
    for (i=0; i<NR; i++) {
        result += vr[i]*r[i];
        for (j=i; j<NR; j++) result += vrr[i][j]*r[i]*r[j];
        for (j=0; j<NS; j++) result += vrs[i][j]*r[i]*s[j];
    }
    for (i=0; i<NS; i++) {
        result += vs[i]*s[i];
        for (j=i; j<NS; j++) result += vss[i][j]*s[i]*s[j];
    }

    return result;
}

static void fillDGDR (double r[NR], double s[NT], double t, double p, double *result) {
    int i, j;

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) {
        result[i] = gr[i] + grr[i][i]*r[i];
        for (j=0; j<NR; j++) result[i] += grr[i][j]*r[j];
        for (j=0; j<NS; j++) result[i] += grs[i][j]*s[j];
    }

    /* Configurational entropy terms */
    for (j=0; j<NR; j++) {
        double config = 0.0;
        for (i=0; i<NE; i++) if (xSpecies[i] > 0.0 && dxSpeciesdr[i][j] != 0.0) config += dxSpeciesdr[i][j]*(1.0 + log(xSpecies[i]));
        if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0 && dxSpeciesdr[nH2O][j] != 0.0)
            config += dxSpeciesdr[nH2O][j]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]));
        result[j] += R*t*nSpecies*config;
    }
}

static void fillDGDS (double r[NR], double s[NT], double t, double p, double *result) {
    int i, j;

    memset(result, '\0', (size_t) NT*sizeof(double));

    /* Taylor expansion and standard state terms */
    for (i=0; i<NS; i++) {
        result[i] = gs[i] + gss[i][i]*s[i];
        for (j=0; j<NS; j++) result[i] += gss[i][j]*s[j];
        for (j=0; j<NR; j++) result[i] += grs[j][i]*r[j];
    }

    /* Configurational entropy terms */
    for (j=0; j<NS; j++) {
        double config = 0.0;
        for (i=0; i<NE; i++) if (xSpecies[i] > 0.0)
            config += dnSpeciesds[j]*xSpecies[i]*log(xSpecies[i]) + nSpecies*dxSpeciesds[i][j]*(1.0 + log(xSpecies[i]));
        if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0)
            config += dnSpeciesds[j]*(xSpecies[nH2O]*log(xSpecies[nH2O]) + (1.0-xSpecies[nH2O])*log(1.0-xSpecies[nH2O]))
                            + nSpecies*dxSpeciesds[nH2O][j]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]));
        result[j] += R*t*config;
    }

}

static void fillDGDW (double r[NR], double s[NT], double t, double p, double *result) {
    int i, j, k, l, m, n;

    /*******************************
    * Parameters: NW WH(), NE H() *
    *             NW WS(), NE S() *
    *             NW WV(), NE V() *
    *******************************/
    memset(result, '\0', (size_t) 3*NP*sizeof(double));

    /**************************************
    * NW W parameters solution are first *
    **************************************/
    for (i=0, n=0; i<NE; i++) for (l=i+1; l<NE; l++, n++) {
        result[     n] +=         taylorCoeff[n+NE][0];
        result[  NP+n] +=      -t*taylorCoeff[n+NE][0];
        result[2*NP+n] += (p-1.0)*taylorCoeff[n+NE][0];
        m = 0;
        for (j=0, m=0; j<NR; j++) {
            result[     n] +=         taylorCoeff[n+NE][1+j]*r[j];
            result[  NP+n] +=      -t*taylorCoeff[n+NE][1+j]*r[j];
            result[2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+j]*r[j];
            for (k=j; k<NR; k++, m++) {
                result[     n] +=         taylorCoeff[n+NE][1+NR+NS+m]*r[j]*r[k];
                result[  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+NS+m]*r[j]*r[k];
                result[2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+NS+m]*r[j]*r[k];
            }
            for (k=0; k<NS; k++, m++) {
                result[     n] +=         taylorCoeff[n+NE][1+NR+NS+m]*r[j]*s[k];
                result[  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+NS+m]*r[j]*s[k];
                result[2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+NS+m]*r[j]*s[k];
            }
        }
        for (j=0; j<NS; j++) {
            result[     n] +=         taylorCoeff[n+NE][1+NR+j]*s[j];
            result[  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+j]*s[j];
            result[2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+j]*s[j];
            for (k=j; k<NS; k++, m++) {
                result[     n] +=         taylorCoeff[n+NE][1+NR+NS+m]*s[j]*s[k];
                result[  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+NS+m]*s[j]*s[k];
                result[2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+NS+m]*s[j]*s[k];
            }
        }
    }

    /**************************************
    * NE standard state terms are second *
    **************************************/
    for (i=0; i<NE; i++, n++) {
        result[        n] +=         taylorCoeff[i][0];
        result[  NP+n] +=      -t*taylorCoeff[i][0];
        result[2*NP+n] += (p-1.0)*taylorCoeff[i][0];
        for (j=0; j<NR; j++) {
            result[        n] +=           taylorCoeff[i][1+j]*r[j];
            result[  NP+n] +=      -t*taylorCoeff[i][1+j]*r[j];
            result[2*NP+n] += (p-1.0)*taylorCoeff[i][1+j]*r[j];
        }
        for (j=0; j<NS; j++) {
            result[        n] +=         taylorCoeff[i][1+NR+j]*s[j];
            result[  NP+n] +=      -t*taylorCoeff[i][1+NR+j]*s[j];
            result[2*NP+n] += (p-1.0)*taylorCoeff[i][1+NR+j]*s[j];
        }
    }
}

static void fillD2GDR2 (double r[NR], double s[NT], double t, double p, double result[NR][NR]) {
    int i, j, k;

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) {
        result[i][i] = 2.0*grr[i][i];
        for (j=i+1; j<NR; j++) {
            result[i][j] = grr[i][j];
            result[j][i] = grr[i][j];
        }
    }

    /* Configurational entropy terms */
    for (j=0; j<NR; j++) {
        for (k=j; k<NR; k++) {
            double config = 0.0;
            for (i=0; i<NE; i++) if (xSpecies[i] > 0.0 && dxSpeciesdr[i][j] != 0.0 && dxSpeciesdr[i][k] != 0.0)
                config += dxSpeciesdr[i][j]*dxSpeciesdr[i][k]/xSpecies[i];
            if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0 && dxSpeciesdr[nH2O][j] != 0.0 && dxSpeciesdr[nH2O][k] != 0.0)
                config += dxSpeciesdr[nH2O][j]*dxSpeciesdr[nH2O][k]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O]));
            result[j][k] += R*t*nSpecies*config;
            result[k][j]  = result[j][k];
        }
    }
}

static void fillD2GDRDS (double r[NR], double s[NT], double t, double p, double result[NR][NT]) {
    int i, j, k;

    for (j=0; j<NR; j++) memset(result[j], '\0', (size_t) NT*sizeof(double));

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) for (j=0; j<NS; j++) result[i][j] = grs[i][j];

    /* Configurational entopy terms */
    for (j=0; j<NR; j++) {
        for (k=0; k<NS; k++) {
            double config = 0.0;
            for (i=0; i<NE; i++) if (xSpecies[i] > 0.0)
                config += nSpecies*(d2xSpeciesdrds[i][j][k]*log(xSpecies[i]) + dxSpeciesdr[i][j]*dxSpeciesds[i][k]/xSpecies[i])
                    + dnSpeciesds[k]*(dxSpeciesdr[i][j]*(1.0 + log(xSpecies[i])));
            if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0)
                config += dnSpeciesds[k]*dxSpeciesdr[nH2O][j]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                    + nSpecies*(d2xSpeciesdrds[nH2O][j][k]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                    + dxSpeciesds[nH2O][k]*dxSpeciesdr[nH2O][j]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O])));
            result[j][k] += R*t*config;
        }
    }
}

static void fillD2GDRDT (double r[NR], double s[NT], double t, double p, double *result) {
    int i, j;

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) {
        result[i] = -sr[i] + -srr[i][i]*r[i];
        for (j=0; j<NR; j++) result[i] += -srr[i][j]*r[j];
        for (j=0; j<NS; j++) result[i] += -srs[i][j]*s[j];
    }

    /* Configurational entropy terms */
    for (j=0; j<NR; j++) {
        double config = 0.0;
        for (i=0; i<NE; i++) if (xSpecies[i] > 0.0 && dxSpeciesdr[i][j] != 0.0) config += dxSpeciesdr[i][j]*(1.0 + log(xSpecies[i]));
        if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0 && dxSpeciesdr[nH2O][j] != 0.0)
            config += dxSpeciesdr[nH2O][j]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]));
        result[j] += R*nSpecies*config;
    }
}

static void fillD2GDRDP (double r[NR], double s[NT], double t, double p, double *result) {
    int i, j;

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) {
        result[i] = vr[i] + vrr[i][i]*r[i];
        for (j=0; j<NR; j++) result[i] += vrr[i][j]*r[j];
        for (j=0; j<NS; j++) result[i] += vrs[i][j]*s[j];
    }
}

static void fillD2GDRDW (double r[NR], double s[NT], double t, double p, double result[NR][3*NP]) {
    int i, j, k, l, m, n, ii;

    /*******************************
    * Parameters: NW WH(), NE H() *
    *             NW WS(), NE S() *
    *             NW WV(), NE V() *
    *******************************/
    for (ii=0; ii<NR; ii++) memset(result[ii], '\0', (size_t) 3*NP*sizeof(double));

    /**************************************
    * NW W parameters solution are first *
    **************************************/
    for (i=0, n=0; i<NE; i++) for (l=i+1; l<NE; l++, n++) {
        for (ii=0; ii<NR; ii++) {
            result[ii][     n] +=        taylorCoeff[n+NE][1+ii];
            result[ii][  NP+n] +=     -t*taylorCoeff[n+NE][1+ii];
            result[ii][2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+ii];
            m = 0;
            for (j=0, m=0; j<NR; j++) {
                   for (k=j; k<NR; k++, m++) {
                    if (j == ii) {
                        result[ii][     n] +=         taylorCoeff[n+NE][1+NR+NS+m]*r[k];
                        result[ii][  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+NS+m]*r[k];
                        result[ii][2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+NS+m]*r[k];
                    }
                    if (k == ii) {
                        result[ii][     n] +=         taylorCoeff[n+NE][1+NR+NS+m]*r[j];
                        result[ii][  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+NS+m]*r[j];
                        result[ii][2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+NS+m]*r[j];
                    }
                }
                for (k=0; k<NS; k++, m++) {
                    if (j == ii) {
                        result[ii][     n] +=         taylorCoeff[n+NE][1+NR+NS+m]*s[k];
                        result[ii][  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+NS+m]*s[k];
                        result[ii][2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+NS+m]*s[k];
                    }
                }
            }
        }
    }

    /**************************************
   * NE standard state terms are second *
   **************************************/
    for (i=0; i<NE; i++, n++) {
        for (ii=0; ii<NR; ii++) {
            result[ii][     n] +=         taylorCoeff[i][1+ii];
            result[ii][  NP+n] +=      -t*taylorCoeff[i][1+ii];
            result[ii][2*NP+n] += (p-1.0)*taylorCoeff[i][1+ii];
        }
    }
}

static void fillD2GDS2New (double r[NR], double s[NT], double t, double p, double **result) {
    int i, j, k;

    for (j=0; j<NT; j++) memset(result[j], '\0', (size_t) NT*sizeof(double));

    /* Taylor expansion and standard state terms */
    for (i=0; i<NS; i++) {
        result[i][i] = 2.0*gss[i][i];
        for (j=i+1; j<NS; j++) {
            result[i][j] = gss[i][j];
            result[j][i] = gss[i][j];
        }
    }

    /* Configurational entropy terms */
    for (j=0; j<NS; j++) {
        for (k=j; k<NS; k++) {
            double config = 0.0;
            for (i=0; i<NE; i++) if (xSpecies[i] > 0.0)
                config += nSpecies*dxSpeciesds[i][j]*dxSpeciesds[i][k]/xSpecies[i]
                        + dnSpeciesds[k]*dxSpeciesds[i][j]*(1.0 + log(xSpecies[i]))
                        + dnSpeciesds[j]*dxSpeciesds[i][k]*(1.0 + log(xSpecies[i]))
                        + d2nSpeciesds2[j][k]*xSpecies[i]*log(xSpecies[i]);
            if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0)
                config += d2nSpeciesds2[j][k]*(xSpecies[nH2O]*log(xSpecies[nH2O]) + (1.0-xSpecies[nH2O])*log(1.0-xSpecies[nH2O]))
                        + dnSpeciesds[j]*dxSpeciesds[nH2O][k]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                        + dnSpeciesds[k]*dxSpeciesds[nH2O][j]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                        + nSpecies*dxSpeciesds[nH2O][j]*dxSpeciesds[nH2O][k]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O]));
            result[j][k] += R*t*config;
            result[k][j]  = result[j][k];
        }
    }
}

static void fillD2GDS2 (double r[NR], double s[NT], double t, double p, double result[NT][NT]) {
    int i, j, k;

    for (j=0; j<NT; j++) memset(result[j], '\0', (size_t) NT*sizeof(double));

    /* Taylor expansion and standard state terms */
    for (i=0; i<NS; i++) {
        result[i][i] = 2.0*gss[i][i];
        for (j=i+1; j<NS; j++) {
            result[i][j] = gss[i][j];
            result[j][i] = gss[i][j];
        }
    }

    /* Configurational entropy terms */
    for (j=0; j<NS; j++) {
        for (k=j; k<NS; k++) {
            double config = 0.0;
            for (i=0; i<NE; i++) if (xSpecies[i] > 0.0)
                config += nSpecies*dxSpeciesds[i][j]*dxSpeciesds[i][k]/xSpecies[i]
                        + dnSpeciesds[k]*dxSpeciesds[i][j]*(1.0 + log(xSpecies[i]))
                        + dnSpeciesds[j]*dxSpeciesds[i][k]*(1.0 + log(xSpecies[i]))
                        + d2nSpeciesds2[j][k]*xSpecies[i]*log(xSpecies[i]);
            if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0)
                config += d2nSpeciesds2[j][k]*(xSpecies[nH2O]*log(xSpecies[nH2O]) + (1.0-xSpecies[nH2O])*log(1.0-xSpecies[nH2O]))
                        + dnSpeciesds[j]*dxSpeciesds[nH2O][k]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                        + dnSpeciesds[k]*dxSpeciesds[nH2O][j]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                        + nSpecies*dxSpeciesds[nH2O][j]*dxSpeciesds[nH2O][k]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O]));
            result[j][k] += R*t*config;
            result[k][j]  = result[j][k];
        }
    }
}

static void fillD2GDSDT (double r[NR], double s[NT], double t, double p, double *result) {
    int i, j;

    memset(result, '\0', (size_t) NT*sizeof(double));

    /* Taylor expansion and standard state terms */
    for (i=0; i<NS; i++) {
        result[i] = -ss[i] - sss[i][i]*s[i];
        for (j=0; j<NS; j++) result[i] += -sss[i][j]*s[j];
        for (j=0; j<NR; j++) result[i] += -srs[j][i]*r[j];
    }

    /* Configurational entropy terms */
    for (j=0; j<NS; j++) {
        double config = 0.0;
        for (i=0; i<NE; i++) if (xSpecies[i] > 0.0)
            config += dnSpeciesds[j]*xSpecies[i]*log(xSpecies[i]) + nSpecies*dxSpeciesds[i][j]*(1.0 + log(xSpecies[i]));
        if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0)
            config += dnSpeciesds[j]*(xSpecies[nH2O]*log(xSpecies[nH2O]) + (1.0-xSpecies[nH2O])*log(1.0-xSpecies[nH2O]))
                            + nSpecies*dxSpeciesds[nH2O][j]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]));
        result[j] += R*config;
    }
}

static void fillD2GDSDP (double r[NR], double s[NT], double t, double p, double *result) {
    int i, j;

    memset(result, '\0', (size_t) NT*sizeof(double));

    /* Taylor expansion and standard state terms */
    for (i=0; i<NS; i++) {
        result[i] = vs[i] + vss[i][i]*s[i];
        for (j=0; j<NS; j++) result[i] += vss[i][j]*s[j];
        for (j=0; j<NR; j++) result[i] += vrs[j][i]*r[j];
    }
}

static void fillD2GDSDW (double r[NR], double s[NT], double t, double p, double result[NT][3*NP]) {
    int i, j, k, l, m, n, ii;

    /*******************************
   * Parameters: NW WH(), NE H() *
   *             NW WS(), NE S() *
   *             NW WV(), NE V() *
   *******************************/
    for (ii=0; ii<NT; ii++) memset(result[ii], '\0', (size_t) 3*NP*sizeof(double));

    /**************************************
   * NW W parameters solution are first *
   **************************************/
    for (i=0, n=0; i<NE; i++) for (l=i+1; l<NE; l++, n++) {
        for (ii=0; ii<NS; ii++) {
            for (j=0, m=0; j<NR; j++) {
                m += NR - j;
                m += ii;
                result[ii][     n] +=         taylorCoeff[n+NE][1+NR+NS+m]*r[j];
                result[ii][  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+NS+m]*r[j];
                result[ii][2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+NS+m]*r[j];
                m += NS - ii;
            }
            result[ii][     n] +=         taylorCoeff[n+NE][1+NR+ii];
            result[ii][  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+ii];
            result[ii][2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+ii];

            for (k=ii; k<NS; k++) {
                m = ii*NS+(k+1)-(ii+1)*(ii+2)/2+(ii+1)-1+NR*(NR-1)/2+NR+NR*NS;
                if (taylorCoeff[n+NE][1+NR+NS+m] != 0.0) {
                    result[ii][     n] +=         taylorCoeff[n+NE][1+NR+NS+m]*s[k];
                    result[ii][  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+NS+m]*s[k];
                    result[ii][2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+NS+m]*s[k];
                }
            }

            for (j=0; j<=ii; j++) {
                m = j*NS+(ii+1)-(j+1)*(j+2)/2+(j+1)-1+NR*(NR-1)/2+NR+NR*NS;
                if (taylorCoeff[n+NE][1+NR+NS+m] != 0.0) {
                    result[ii][     n] +=         taylorCoeff[n+NE][1+NR+NS+m]*s[j];
                    result[ii][  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+NS+m]*s[j];
                    result[ii][2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+NS+m]*s[j];
                }
            }

        }
    }

    /**************************************
    * NE standard state terms are second *
    **************************************/
    for (i=0; i<NE; i++, n++) {
        for (ii=0; ii<NS; ii++) {
            result[ii][     n] +=         taylorCoeff[i][1+NR+ii];
            result[ii][  NP+n] +=      -t*taylorCoeff[i][1+NR+ii];
            result[ii][2*NP+n] += (p-1.0)*taylorCoeff[i][1+NR+ii];
        }
    }
}

static double fillD2GDT2 (double r[NR], double s[NT], double t, double p) {
    double result;
    int i;

    /* Taylor expansion and standard state terms */
    result = CPconst;
    for (i=0; i<NR; i++) result += cpr[i]*r[i];
    for (i=0; i<NS; i++) result += cps[i]*s[i];
    result /= -t;

    return result;
}

static double fillD2GDTDP (double r[NR], double s[NT], double t, double p) {
    double result;
    int i;

    /* Taylor expansion and standard state terms */
    result = DVDTconst;
    for (i=0; i<NR; i++) result += dvdtr[i]*r[i];
    for (i=0; i<NS; i++) result += dvdts[i]*s[i];

    return result;
}

static void fillD2GDTDW (double r[NR], double s[NT], double t, double p, double result[3*NP]) {
    int i, j, k, l, m, n;

    /*******************************
    * Parameters: NW WH(), NE H() *
    *             NW WS(), NE S() *
    *             NW WV(), NE V() *
    *******************************/
    memset(result, '\0', (size_t) 3*NP*sizeof(double));

    /**************************************
    * NW W parameters solution are first *
    **************************************/
    for (i=0, n=0; i<NE; i++) for (l=i+1; l<NE; l++, n++) {
        result[NP+n] += -taylorCoeff[n+NE][0];
        m = 0;
        for (j=0, m=0; j<NR; j++) {
            result[NP+n] += -taylorCoeff[n+NE][1+j]*r[j];
            for (k=j; k<NR; k++, m++) result[NP+n] += -taylorCoeff[n+NE][1+NR+NS+m]*r[j]*r[k];
            for (k=0; k<NS; k++, m++) result[NP+n] += -taylorCoeff[n+NE][1+NR+NS+m]*r[j]*s[k];
        }
        for (j=0; j<NS; j++) {
            result[NP+n] += -taylorCoeff[n+NE][1+NR+j]*s[j];
            for (k=j; k<NS; k++, m++) result[NP+n] += -taylorCoeff[n+NE][1+NR+NS+m]*s[j]*s[k];
        }
    }

    /**************************************
    * NE standard state terms are second *
    **************************************/
    for (i=0; i<NE; i++, n++) {
        result[NP+n] += -taylorCoeff[i][0];
        for (j=0; j<NR; j++) result[NP+n] += -taylorCoeff[i][1+j]*r[j];
        for (j=0; j<NS; j++) result[NP+n] += -taylorCoeff[i][1+NR+j]*s[j];
    }
}

static void fillD2GDPDW (double r[NR], double s[NT], double t, double p, double *result) {

    /*******************************
    * Parameters: NW WH(), NE H() *
    *             NW WS(), NE S() *
    *             NW WV(), NE V() *
    *******************************/
    memset(result, '\0', (size_t) 3*NP*sizeof(double));

    /**************************************
    * NW W parameters solution are first *
    **************************************/
{
    int i, j, k, l, m, n;
    for (i=0, n=0; i<NE; i++) for (l=i+1; l<NE; l++, n++) {
        result[2*NP+n] += taylorCoeff[n+NE][0];
        m = 0;
        for (j=0, m=0; j<NR; j++) {
            result[2*NP+n] += taylorCoeff[n+NE][1+j]*r[j];
            for (k=j; k<NR; k++, m++) result[2*NP+n] += taylorCoeff[n+NE][1+NR+NS+m]*r[j]*r[k];
            for (k=0; k<NS; k++, m++) result[2*NP+n] += taylorCoeff[n+NE][1+NR+NS+m]*r[j]*s[k];
        }
        for (j=0; j<NS; j++) {
            result[2*NP+n] += taylorCoeff[n+NE][1+NR+j]*s[j];
            for (k=j; k<NS; k++, m++) result[2*NP+n] += taylorCoeff[n+NE][1+NR+NS+m]*s[j]*s[k];
        }
    }

    /**************************************
    * NE standard state terms are second *
    **************************************/
    for (i=0; i<NE; i++, n++) {
        result[2*NP+n] += taylorCoeff[i][0];
        for (j=0; j<NR; j++) result[2*NP+n] += taylorCoeff[i][1+j]*r[j];
        for (j=0; j<NS; j++) result[2*NP+n] += taylorCoeff[i][1+NR+j]*s[j];
    }
}
}

static double fillD2GDP2 (double r[NR], double s[NT], double t, double p) {
    double result;
    int i;

    /* Taylor expansion and standard state terms */
    result = DVDPconst;
    for (i=0; i<NR; i++) result += dvdpr[i]*r[i];
    for (i=0; i<NS; i++) result += dvdps[i]*s[i];

    return result;
}

static void fillD3GDR3 (double r[NR], double s[NT], double t, double p, double result[NR][NR][NR]) {
    int i, j, k, l;

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) for (j=0; j<NR; j++) memset(result[i][j], '\0', (size_t) NR*sizeof(double));

    /* Configurational entropy terms */
    for (j=0; j<NR; j++) {
        for (k=j; k<NR; k++) {
            for (l=k; l<NR; l++) {
                double config = 0.0;
                for (i=0; i<NE; i++) if (xSpecies[i] > 0.0 && dxSpeciesdr[i][j] != 0.0 && dxSpeciesdr[i][k] != 0.0 && dxSpeciesdr[i][l] != 0.0)
                    config += -dxSpeciesdr[i][j]*dxSpeciesdr[i][k]*dxSpeciesdr[i][l]/(xSpecies[i]*xSpecies[i]);
                if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0 && dxSpeciesdr[nH2O][j] != 0.0 && dxSpeciesdr[nH2O][k] != 0.0 && dxSpeciesdr[i][l] != 0.0)
                    config += -dxSpeciesdr[nH2O][j]*dxSpeciesdr[nH2O][k]*dxSpeciesdr[nH2O][l]
                *(1.0/(xSpecies[nH2O]*xSpecies[nH2O]) - 1.0/((1.0-xSpecies[nH2O])*(1.0-xSpecies[nH2O])));
                result[j][k][l] += R*t*nSpecies*config;
                result[k][j][l]  = result[j][k][l];
                result[l][j][k]  = result[j][k][l];
                result[l][k][j]  = result[j][k][l];
                result[j][l][k]  = result[j][k][l];
                result[k][l][j]  = result[j][k][l];
            }
        }
    }
}

static void fillD3GDR2DS (double r[NR], double s[NT], double t, double p, double result[NR][NR][NT]) {
    int i, j, k, l;

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) for (j=0; j<NR; j++) memset(result[i][j], '\0', (size_t) NT*sizeof(double));

    /* Configurational entropy terms */
    for (j=0; j<NR; j++) {
        for (k=j; k<NR; k++) {
            for (l=0; l<NS; l++) {
                double config = 0.0;
                for (i=0; i<NE; i++) if (xSpecies[i] > 0.0)
                    config += nSpecies*(
                        (d2xSpeciesdrds[i][j][l]*dxSpeciesdr[i][k]+dxSpeciesdr[i][j]*d2xSpeciesdrds[i][k][l])/xSpecies[i]
                        - dxSpeciesdr[i][j]*dxSpeciesdr[i][k]*dxSpeciesds[i][l]/(xSpecies[i]*xSpecies[i]) )
                        + dnSpeciesds[l]*dxSpeciesdr[i][j]*dxSpeciesdr[i][k]/xSpecies[i];
                if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0)
                    config += nSpecies*( (d2xSpeciesdrds[nH2O][j][l]*dxSpeciesdr[nH2O][k] + dxSpeciesdr[nH2O][j]*d2xSpeciesdrds[nH2O][k][l])
                            *(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O]))
                        - dxSpeciesdr[nH2O][j]*dxSpeciesdr[nH2O][k]*dxSpeciesds[nH2O][l]
                            *(1.0/(xSpecies[nH2O]*xSpecies[nH2O]) - 1.0/((1.0-xSpecies[nH2O])*(1.0-xSpecies[nH2O]))) )
                        + dnSpeciesds[l]*dxSpeciesdr[nH2O][j]*dxSpeciesdr[nH2O][k]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O]));

                result[j][k][l] += R*t*config;
                result[k][j][l]  = result[j][k][l];
            }
        }
    }
}

static void fillD3GDR2DT (double r[NR], double s[NT], double t, double p, double result[NR][NR]) {
    int i, j, k;

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) for (j=i; j<NR; j++) {
        result[i][j] = (i == j) ? -2.0*srr[i][i] : -srr[i][j];
        result[j][i] = result[i][j];
    }

    /* Configurational entropy terms */
    for (j=0; j<NR; j++) {
        for (k=j; k<NR; k++) {
            double config = 0.0;
            for (i=0; i<NE; i++) if (xSpecies[i] > 0.0 && dxSpeciesdr[i][j] != 0.0 && dxSpeciesdr[i][k] != 0.0)
                config += dxSpeciesdr[i][j]*dxSpeciesdr[i][k]/xSpecies[i];
            if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0 && dxSpeciesdr[nH2O][j] != 0.0 && dxSpeciesdr[nH2O][k] != 0.0)
                config += dxSpeciesdr[nH2O][j]*dxSpeciesdr[nH2O][k]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O]));
            result[j][k] += R*nSpecies*config;
            result[k][j]  = result[j][k];
        }
    }
}

static void fillD3GDR2DP (double r[NR], double s[NT], double t, double p, double result[NR][NR]) {
    int i, j;

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) for (j=i; j<NR; j++) {
        result[i][j] = (i == j) ? 2.0*vrr[i][i] : vrr[i][j];
        result[j][i] = result[i][j];
    }

}

static void fillD3GDRDS2 (double r[NR], double s[NT], double t, double p, double result[NR][NT][NT]) {
    int i, j, k, l;

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) for (j=0; j<NT; j++) memset(result[i][j], '\0', (size_t) NT*sizeof(double));

    /* Configurational entropy terms */
    for (j=0; j<NR; j++) {
        for (k=0; k<NS; k++) {
            for (l=k; l<NS; l++) {
                double config = 0.0;
                for (i=0; i<NE; i++) if (xSpecies[i] > 0.0)
                    config += nSpecies*(
                        (d2xSpeciesdrds[i][j][k]*dxSpeciesds[i][l] + dxSpeciesds[i][k]*d2xSpeciesdrds[i][j][l])/xSpecies[i]
                        - dxSpeciesdr[i][j]*dxSpeciesds[i][k]*dxSpeciesds[i][l]/(xSpecies[i]*xSpecies[i]) )
                        + dnSpeciesds[l]*(d2xSpeciesdrds[i][j][k]*log(xSpecies[i]) + dxSpeciesdr[i][j]*dxSpeciesds[i][k]/xSpecies[i])
                        + dnSpeciesds[k]*(d2xSpeciesdrds[i][j][l]*log(xSpecies[i]) + dxSpeciesdr[i][j]*dxSpeciesds[i][l]/xSpecies[i])
                        + d2nSpeciesds2[k][l]*(dxSpeciesdr[i][j]*(1.0 + log(xSpecies[i])));
                if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0)
                    config += d2nSpeciesds2[k][l]*(dxSpeciesdr[nH2O][j]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O])))
                        + dnSpeciesds[k]*(d2xSpeciesdrds[nH2O][j][l]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                        + dxSpeciesdr[nH2O][j]*dxSpeciesds[nH2O][l]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O])))
                        + dnSpeciesds[l]*(d2xSpeciesdrds[nH2O][j][k]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                        + dxSpeciesds[nH2O][k]*dxSpeciesdr[nH2O][j]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O])))
                        + nSpecies*(
                            d2xSpeciesdrds[nH2O][j][k]*dxSpeciesds[nH2O][l]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O]))
                            + dxSpeciesds[nH2O][k]*d2xSpeciesdrds[nH2O][j][l]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O]))
                            - dxSpeciesds[nH2O][k]*dxSpeciesdr[nH2O][j]*dxSpeciesds[nH2O][l]
                                *(1.0/(xSpecies[nH2O]*xSpecies[nH2O]) - 1.0/((1.0-xSpecies[nH2O])*(1.0-xSpecies[nH2O]))) );
                result[j][k][l] += R*t*config;
                result[j][l][k]  = result[j][k][l];
            }
        }
    }
}

static void fillD3GDRDSDT (double r[NR], double s[NT], double t, double p, double result[NR][NT]) {
    int i, j, k;

    for (i=0; i<NR; i++) memset(result[i], '\0', (size_t) NT*sizeof(double));

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) for (j=0; j<NS; j++) result[i][j] = -srs[i][j];

    /* Configurational entropy terms */
    for (j=0; j<NR; j++) {
        for (k=0; k<NS; k++) {
            double config = 0.0;
            for (i=0; i<NE; i++) if (xSpecies[i] > 0.0)
                config += nSpecies*(d2xSpeciesdrds[i][j][k]*log(xSpecies[i]) + dxSpeciesdr[i][j]*dxSpeciesds[i][k]/xSpecies[i])
                    + dnSpeciesds[k]*(dxSpeciesdr[i][j]*(1.0 + log(xSpecies[i])));
            if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0)
                config += dnSpeciesds[k]*dxSpeciesdr[nH2O][j]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                    + nSpecies*(d2xSpeciesdrds[nH2O][j][k]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                    + dxSpeciesds[nH2O][k]*dxSpeciesdr[nH2O][j]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O])));
            result[j][k] += R*config;
        }
    }
}

static void fillD3GDRDSDP (double r[NR], double s[NT], double t, double p, double result[NR][NT]) {
    int i, j;

    for (i=0; i<NR; i++) memset(result[i], '\0', (size_t) NT*sizeof(double));

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) for (j=0; j<NS; j++) result[i][j] = vrs[i][j];
}

static void fillD3GDS3 (double r[NR], double s[NT], double t, double p, double result[NT][NT][NT]) {
    int i, j, k, l;

    /* Taylor expansion and standard state terms */
    for (i=0; i<NT; i++) for (j=0; j<NT; j++) for (k=0; k<NT; k++) result[i][j][k] = 0.0;

    /* Configurational entropy terms */
    for (j=0; j<NS; j++) {
        for (k=j; k<NS; k++) {
            for (l=k; l<NS; l++) {
                double config = 0.0;
                for (i=0; i<NE; i++) if (xSpecies[i] > 0.0)
                    config += -nSpecies*dxSpeciesds[i][j]*dxSpeciesds[i][k]*dxSpeciesds[i][l]/(xSpecies[i]*xSpecies[i])
                        + dnSpeciesds[l]*dxSpeciesds[i][j]*dxSpeciesds[i][k]/xSpecies[i]
                        + dnSpeciesds[k]*dxSpeciesds[i][j]*dxSpeciesds[i][l]/xSpecies[i]
                        + d2nSpeciesds2[k][l]*dxSpeciesds[i][j]*(1.0 + log(xSpecies[i]))
                        + dnSpeciesds[j]*dxSpeciesds[i][k]*dxSpeciesds[i][l]/xSpecies[i]
                        + d2nSpeciesds2[j][l]*dxSpeciesds[i][k]*(1.0 + log(xSpecies[i]))
                        + d2nSpeciesds2[j][k]*dxSpeciesds[i][l]*(1.0 + log(xSpecies[i]))
                        + d3nSpeciesds3[j][k][l]*xSpecies[i]*log(xSpecies[i]);
                if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0)
                    config += d3nSpeciesds3[j][k][l]*(xSpecies[nH2O]*log(xSpecies[nH2O]) + (1.0-xSpecies[nH2O])*log(1.0-xSpecies[nH2O]))
                        + d2nSpeciesds2[j][k]*dxSpeciesds[nH2O][l]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                        + d2nSpeciesds2[j][l]*dxSpeciesds[nH2O][k]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                        + dnSpeciesds[j]*dxSpeciesds[nH2O][k]*dxSpeciesds[nH2O][l]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O]))
                        + d2nSpeciesds2[k][l]*dxSpeciesds[nH2O][j]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                        + dnSpeciesds[k]*dxSpeciesds[nH2O][j]*dxSpeciesds[nH2O][l]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O]))
                        + dnSpeciesds[l]*dxSpeciesds[nH2O][j]*dxSpeciesds[nH2O][k]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O]))
                        - nSpecies*dxSpeciesds[nH2O][j]*dxSpeciesds[nH2O][k]*dxSpeciesds[nH2O][l]
                        *(1.0/(xSpecies[nH2O]*xSpecies[nH2O]) - 1.0/((1.0-xSpecies[nH2O])*(1.0-xSpecies[nH2O])));
                result[j][k][l] += R*t*config;
                result[k][j][l]  = result[j][k][l];
                result[l][j][k]  = result[j][k][l];
                result[l][k][j]  = result[j][k][l];
                result[j][l][k]  = result[j][k][l];
                result[k][l][j]  = result[j][k][l];
            }
        }
    }
}

static void fillD3GDS2DT (double r[NR], double s[NT], double t, double p, double result[NT][NT]) {
    int i, j, k;

    for (i=0; i<NS; i++) memset(result[i], '\0', (size_t) NS*sizeof(double));

    /* Taylor expansion and standard state terms */
    for (i=0; i<NS; i++) for (j=i; j<NS; j++) {
            result[i][j] = (i == j) ? -2.0*sss[i][i] : -sss[i][j];
            result[j][i] = result[i][j];
    }

    /* Configurational entropy terms */
    for (j=0; j<NS; j++) {
        for (k=j; k<NS; k++) {
            double config = 0.0;
            for (i=0; i<NE; i++) if (xSpecies[i] > 0.0)
                config += nSpecies*dxSpeciesds[i][j]*dxSpeciesds[i][k]/xSpecies[i]
                    + dnSpeciesds[k]*dxSpeciesds[i][j]*(1.0 + log(xSpecies[i]))
                    + dnSpeciesds[j]*dxSpeciesds[i][k]*(1.0 + log(xSpecies[i]))
                    + d2nSpeciesds2[j][k]*xSpecies[i]*log(xSpecies[i]);
            if (nH2O != -1 && xSpecies[nH2O] > 0.0 && xSpecies[nH2O] < 1.0)
                config += d2nSpeciesds2[j][k]*(xSpecies[nH2O]*log(xSpecies[nH2O]) + (1.0-xSpecies[nH2O])*log(1.0-xSpecies[nH2O]))
                    + dnSpeciesds[j]*dxSpeciesds[nH2O][k]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                    + dnSpeciesds[k]*dxSpeciesds[nH2O][j]*(log(xSpecies[nH2O]) - log(1.0-xSpecies[nH2O]))
                    + nSpecies*dxSpeciesds[nH2O][j]*dxSpeciesds[nH2O][k]*(1.0/xSpecies[nH2O] + 1.0/(1.0-xSpecies[nH2O]));
            result[j][k] += R*config;
            result[k][j]  = result[j][k];
        }
    }
}

static void fillD3GDS2DP (double r[NR], double s[NT], double t, double p, double result[NT][NT]) {
    int i, j;

    for (i=0; i<NS; i++) memset(result[i], '\0', (size_t) NS*sizeof(double));

    /* Taylor expansion and standard state terms */
    for (i=0; i<NS; i++) for (j=i; j<NS; j++) {
            result[i][j] = (i == j) ? 2.0*vss[i][i] : vss[i][j];
            result[j][i] = result[i][j];
    }
}

static void fillD3GDSDT2 (double r[NR], double s[NT], double t, double p, double *result) {
    int i;

    memset(result, '\0', (size_t) NT*sizeof(double));

    /* Taylor expansion and standard state terms */
    for (i=0; i<NS; i++) result[i] = -cps[i]/t;
}

static void fillD3GDSDTDP (double r[NR], double s[NT], double t, double p, double *result) {
    int i;

    memset(result, '\0', (size_t) NT*sizeof(double));

    /* Taylor expansion and standard state terms */
    for (i=0; i<NS; i++) result[i] = dvdts[i];
}

static void fillD3GDSDP2 (double r[NR], double s[NT], double t, double p, double *result) {
    int i;

    memset(result, '\0', (size_t) NT*sizeof(double));

    /* Taylor expansion and standard state terms */
    for (i=0; i<NS; i++) result[i] = dvdps[i];
}

static void fillD3GDRDT2 (double r[NR], double s[NT], double t, double p, double *result) {
    int i;

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) result[i] = -cpr[i]/t;
}

static void fillD3GDRDTDP (double r[NR], double s[NT], double t, double p, double *result) {
    int i;

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) result[i] = dvdtr[i];
}

static void fillD3GDRDP2 (double r[NR], double s[NT], double t, double p, double *result) {
    int i;

    /* Taylor expansion and standard state terms */
    for (i=0; i<NR; i++) result[i] = dvdpr[i];
}

static double fillD3GDT3 (double r[NR], double s[NT], double t, double p) {
    double result;
    int i;

    /* Taylor expansion and standard state terms */
    result = CPconst/(t*t) - DCPDTconst/t;
    for (i=0; i<NR; i++) result += cpr[i]*r[i]/(t*t) - dcpdtr[i]*r[i]/t;
    for (i=0; i<NS; i++) result += cps[i]*s[i]/(t*t) - dcpdts[i]*s[i]/t;

    return result;
}

static double fillD3GDT2DP (double r[NR], double s[NT], double t, double p) {
    double result;
    int i;

    /* Taylor expansion and standard state terms */
    result = D2VDT2const;
    for (i=0; i<NR; i++) result += d2vdt2r[i]*r[i];
    for (i=0; i<NS; i++) result += d2vdt2s[i]*s[i];

    return result;
}

static double fillD3GDTDP2 (double r[NR], double s[NT], double t, double p) {
    double result;
    int i;

    /* Taylor expansion and standard state terms */
    result = D2VDTDPconst;
    for (i=0; i<NR; i++) result += d2vdtdpr[i]*r[i];
    for (i=0; i<NS; i++) result += d2vdtdps[i]*s[i];

    return result;
}

static double fillD3GDP3 (double r[NR], double s[NT], double t, double p) {
    double result;
    int i;

    /* Taylor expansion and standard state terms */
    result = D2VDP2const;
    for (i=0; i<NR; i++) result += d2vdp2r[i]*r[i];
    for (i=0; i<NS; i++) result += d2vdp2s[i]*s[i];

    return result;
}

static void fillD3GDS2DW (double r[NR], double s[NT], double t, double p, double result[NT][NT][3*NP]) {
    int i, j, k, l, m, n, ii, iii;

    /*******************************
    * Parameters: NW WH(), NE H() *
    *             NW WS(), NE S() *
    *             NW WV(), NE V() *
    *******************************/
    for (k=0; k<NT; k++) for (j=0; j<NT; j++) memset(result[k][j], '\0', (size_t) 3*NP*sizeof(double));

    /**************************************
    * NW W parameters solution are first *
    **************************************/
    for (i=0, n=0; i<NE; i++) for (l=i+1; l<NE; l++, n++) {
        for (ii=0; ii<NS; ii++) {
            for (iii=ii; iii<NS; iii++) {
                m = ii*NS+(iii+1)-(ii+1)*(ii+2)/2+(ii+1)-1+NR*(NR-1)/2+NR+NR*NS;
                if (taylorCoeff[n+NE][1+NR+NS+m] != 0.0) {
                    result[ii][iii][     n] +=         taylorCoeff[n+NE][1+NR+NS+m];
                    result[ii][iii][  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+NS+m];
                    result[ii][iii][2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+NS+m];
                    result[iii][ii][     n] +=         taylorCoeff[n+NE][1+NR+NS+m];
                    result[iii][ii][  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+NS+m];
                    result[iii][ii][2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+NS+m];
                }
            }
        }
    }
}

static void fillD3GDSDTDW (double r[NR], double s[NT], double t, double p, double result[NT][3*NP]) {
    int i, j, k, l, m, n, ii;

    /*******************************
    * Parameters: NW WH(), NE H() *
    *             NW WS(), NE S() *
    *             NW WV(), NE V() *
    *******************************/
    for (ii=0; ii<NT; ii++) memset(result[ii], '\0', (size_t) 3*NP*sizeof(double));

    /**************************************
    * NW W parameters solution are first *
    **************************************/
    for (i=0, n=0; i<NE; i++) for (l=i+1; l<NE; l++, n++) {
        for (ii=0; ii<NS; ii++) {
            for (j=0, m=0; j<NR; j++) {
                m += NR - j;
                m += ii;
                result[ii][  NP+n] += -taylorCoeff[n+NE][1+NR+NS+m]*r[j];
                m += NS - ii;
            }
            result[ii][  NP+n] += -taylorCoeff[n+NE][1+NR+ii];
            for (k=ii; k<NS; k++) {
                m = ii*NS+(k+1)-(ii+1)*(ii+2)/2+(ii+1)-1+NR*(NR-1)/2+NR+NR*NS;
                result[ii][  NP+n] += -taylorCoeff[n+NE][1+NR+NS+m]*s[k];
            }
            for (j=0; j<=ii; j++) {
                m = j*NS+(ii+1)-(j+1)*(j+2)/2+(j+1)-1+NR*(NR-1)/2+NR+NR*NS;
                result[ii][  NP+n] += -taylorCoeff[n+NE][1+NR+NS+m]*s[j];
            }
        }
    }

    /**************************************
    * NE standard state terms are second *
    **************************************/
    for (i=0; i<NE; i++, n++) {
        for (ii=0; ii<NS; ii++) result[ii][  NP+n] += -taylorCoeff[i][1+NR+ii];
    }
}

static void fillD3GDRDTDW (double r[NR], double s[NT], double t, double p, double result[NR][3*NP]) {
    int i, j, k, l, m, n, ii;

    /*******************************
    * Parameters: NW WH(), NE H() *
    *             NW WS(), NE S() *
    *             NW WV(), NE V() *
    *******************************/
    for (ii=0; ii<NR; ii++) memset(result[ii], '\0', (size_t) 3*NP*sizeof(double));

    /**************************************
    * NW W parameters solution are first *
    **************************************/
    for (i=0, n=0; i<NE; i++) for (l=i+1; l<NE; l++, n++) {
        for (ii=0; ii<NR; ii++) {
            result[ii][  NP+n] += -taylorCoeff[n+NE][1+ii];
            m = 0;
            for (j=0, m=0; j<NR; j++) {
                for (k=j; k<NR; k++, m++) {
                    if (j == ii) result[ii][  NP+n] += -taylorCoeff[n+NE][1+NR+NS+m]*r[k];
                    if (k == ii) result[ii][  NP+n] += -taylorCoeff[n+NE][1+NR+NS+m]*r[j];
                }
                for (k=0; k<NS; k++, m++) {
                    if (j == ii) result[ii][  NP+n] += -taylorCoeff[n+NE][1+NR+NS+m]*s[k];
                }
            }
        }
    }

    /**************************************
    * NE standard state terms are second *
    **************************************/
    for (i=0; i<NE; i++, n++) {
        for (ii=0; ii<NR; ii++) result[ii][  NP+n] += -taylorCoeff[i][1+ii];
    }
}

static void fillD3GDRDSDW (double r[NR], double s[NT], double t, double p, double result[NR][NT][3*NP]) {
    int i, l, m, n, ii, iii;

    /*******************************
    * Parameters: NW WH(), NE H() *
    *             NW WS(), NE S() *
    *             NW WV(), NE V() *
    *******************************/
    for (iii=0; iii<NR; iii++) for (ii=0; ii<NT; ii++) memset(result[iii][ii], '\0', (size_t) 3*NP*sizeof(double));

    /**************************************
    * NW W parameters solution are first *
    **************************************/
    for (i=0, n=0; i<NE; i++) for (l=i+1; l<NE; l++, n++) {
        for (ii=0; ii<NR; ii++) {
            for (iii=0; iii<NS; iii++) {
                m = (ii+1)*NR-ii*(ii+1)/2+iii*(ii+1)+(NS-iii)*ii;
                result[ii][iii][     n] +=         taylorCoeff[n+NE][1+NR+NS+m];
                result[ii][iii][  NP+n] +=      -t*taylorCoeff[n+NE][1+NR+NS+m];
                result[ii][iii][2*NP+n] += (p-1.0)*taylorCoeff[n+NE][1+NR+NS+m];
            }
        }
    }
}

static void initialGuessOrdering(double r[NR], double s[NT]) {
    int i;
    static double *sCorr;
    double factor = 1.0;
#ifndef TESTDYNAMICLIB
    double rSum, s0;
#endif

    /* Update to match Thermoengine */
    double denom = 1.0 + (((r[9] > 0.0) && (r[13] > 0.0)) ? 1.0 : 0.0); /* [ 0] SiO2 [10] CaO [14] CO2 */
    for (i=0; i<NS; i++)  s[i] = r[13]/denom; // if CaO  is present 1/denom of all CO2 is CaCO3

    if (NT == 0) return;

#ifdef DEBUG
    printf("Call to initialGuessOrdering in liquid.c\n");
#endif
    if (sCorr == NULL) sCorr = (double *) malloc((size_t) NS*sizeof(double));

    /*for (i=0; i<NS; i++)  s[i] = 0.0;
    for (i=NS; i<NT; i++) s[i] = 1.0/(((double) NY)+1.0);*/

#ifndef TESTDYNAMICLIB
    /* A tweak for alkali-rich liquids where SiO2 tends to have negative mole fraction
        in the initial guess. Used for dynamically loaded library to help avoid doAbort. */
    for (i=0, rSum=0.0; i<NR; i++) rSum += r[i];
    /*coeff = 1.0;
    xSpecies[ 0] = 1.0 - rSum*coeff;*/                   /* SiO2  */
    if (rSum > 1.0) s[0] = rSum - 1.0 + sqrt(DBL_EPSILON);
#endif

    if (!rANDsTOx (r, s)) {
        /* [ 0] SiO2 [10] CaO [14] CO2 */
        static double tolerance = 0.0;

        if (tolerance == 0.0) tolerance = pow(DBL_EPSILON, (double) (2.0/3.0));

        /***********************************/
        /* Find minimum length Soln vector */
        /***********************************/

        /* objective fucntion row for finding minimum length solution vector that satisfies constraints */

        /***********************************/
        /* Find maximum length Soln vector */
        /***********************************/

        /* objective fucntion row for finding maximum length solution vector that satisfies constraints */

        /***************************************************************/
        /* Return average of two solutions - should always be feasible */
        /***************************************************************/

        for (i=0; i<NS; i++) {
            //s[i] = (sMin[i] + sMax[i])/2.0;
            if (fabs(s[i]) < tolerance) s[i] = 0.0;
        }

        if(!rANDsTOx (r, s)) {
#ifdef DEBUG
            fprintf(stderr, "Failed to find feasible solution in initialGuessOrdering.\n");
#endif
        }

        //for (i=NS; i<NT; i++) s[i] = 1.0/(((double) NY)+1.0);

#ifdef DEBUG
        printf("Results of call to initialGuessOrdering:\n");
        printf("   %20.20s %13.13s %13.13s %13.13s\n", "Species", "Mole frac", "r", "s");
        printf("   %20.20s %13.6g\n", liquid[0].label, xSpecies[0]);
        for (i=0;  i<NR; i++) printf("   %20.20s %13.6g %13.6g\n", liquid[i+1].label, xSpecies[i+1], r[i]);
        for (i=0;  i<NS; i++) printf("   %20.20s %13.6g %13.13s %13.6g\n", liquid[i+NA].label, xSpecies[i+NA], "", s[i]);
        for (i=NS; i<NT; i++) printf("   %20.20s %13.13s %13.13s %13.6g\n", "order CN[*]", "", "", s[i]);
#endif

        //return;
    } /* end block on simplex method */

    for (i=0; i<NS; i++) {
#ifndef TESTDYNAMICLIB
        s0 = MAX(rSum - 1.0, 0.0);
        sCorr[i] = s0; s[i] = sCorr[i] + sqrt(DBL_EPSILON);
        if (!rANDsTOx (r, s)) s[i] = sCorr[i];
        else {
            sCorr[i] = 0.5; s[i] = 1.0;
            while (sCorr[i] > sqrt(DBL_EPSILON)) {
                if(!rANDsTOx (r, s)) s[i] -= sCorr[i]; else s[i] += sCorr[i];
                sCorr[i] /= 2.0;
            }
        }
        sCorr[i] = (s[i] + s0)/2.0;
        s[i]     = 0.0;
    }

    for (i=0; i<NS; i++) s[i] = sCorr[i];
    for (i=0; i<NS; i++) sCorr[i] -= s0;
    while (factor > sqrt(DBL_EPSILON) && !rANDsTOx (r, s)) {
        factor /= 2.0;
        for (i=0; i<NS; i++) s[i] = s0 + sCorr[i]*factor;
    }

#else
        sCorr[i] = 0.0; s[i] = sqrt(DBL_EPSILON);
        if (!rANDsTOx (r, s)) s[i] = 0.0;
        else {
            sCorr[i] = 0.5; s[i] = 1.0;
            while (sCorr[i] > sqrt(DBL_EPSILON)) {
                if(!rANDsTOx (r, s)) s[i] -= sCorr[i]; else s[i] += sCorr[i];
                sCorr[i] /= 2.0;
            }
        }
        sCorr[i] = s[i]/2.0;
        s[i]     = 0.0;
    }

    for (i=0; i<NS; i++) s[i] = sCorr[i];
    while (factor > sqrt(DBL_EPSILON) && !rANDsTOx (r, s)) {
        factor /= 2.0;
        for (i=0; i<NS; i++) s[i] = sCorr[i]*factor;
    }

#endif
    //for (i=NS; i<NT; i++) s[i] = 1.0/(((double) NY)+1.0);

#ifdef DEBUG
    printf("Results of call to initialGuessOrdering:\n");
    printf("   %20.20s %13.13s %13.13s %13.13s\n", "Species", "Mole frac", "r", "s");
    printf("   %20.20s %13.6g\n", liquid[0].label, xSpecies[0]);
    for (i=0;  i<NR; i++) printf("   %20.20s %13.6g %13.6g\n", liquid[i+1].label, xSpecies[i+1], r[i]);
    for (i=0;  i<NS; i++) printf("   %20.20s %13.6g %13.13s %13.6g\n", liquid[i+NA].label, xSpecies[i+NA], "", s[i]);
    for (i=NS; i<NT; i++) printf("   %20.20s %13.13s %13.13s %13.6g\n", "order CN[*]", "", "", s[i]);
#endif
}

static void
order(int mask, double t, double p, double r[NR],
            double s[NT],            /* s[NT]                  BINARY MASK: 000000000001 */
            double dr[NT][NR] ,      /* ds[NT]/dr[NR]          BINARY MASK: 000000000010 */
            double dt[NT],           /* ds[NT]/dt              BINARY MASK: 000000000100 */
            double dp[NT],           /* ds[NT]/dp              BINARY MASK: 000000001000 */
            double dr2[NT][NR][NR],  /* d2s[NT]/dr[NR]dr[NR]   BINARY MASK: 000000010000 */
            double drt[NT][NR],      /* d2s[NT]/dr[NR]dt       BINARY MASK: 000000100000 */
            double drp[NT][NR],      /* d2s[NT]/dr[NR]dp       BINARY MASK: 000001000000 */
            double dt2[NT],          /* d2s[NT]/dt2            BINARY MASK: 000010000000 */
            double dtp[NT],          /* d2s[NT]/dtp            BINARY MASK: 000100000000 */
            double dp2[NT],          /* d2s[NT]/dp2            BINARY MASK: 001000000000 */
            double dw[NT][3*NP],     /* ds[NT]/dw[3*NP]        BINARY MASK: 010000000000 */
            double dtw[NT][3*NP]     /* ds[NT]/dtdw[3*NP]      BINARY MASK: 100000000000 */
     )
{
    double tOld         = getTOld();
    double pOld         = getPOld();
    double *rOld        = getROld();
    double *sOld        = getSOld();
    double **d2gds2     = getD2gds2();
    gsl_matrix       *ptToD2gds2 = getPtToD2gds2();
    gsl_matrix       *ptToV = getPtToV();
    gsl_vector       *ptToS = getPtToS();
    gsl_permutation  *indexLU = getIndexD2gds2();
    int i, j, iter = 0, doAbort=FALSE, update=FALSE, loop, dLU;

    update |= (t != tOld);
    update |= (p != pOld);
    if (update) loadTaylorCoefficients(t, p);         /* if T, P, or meltsAndCO2_H2OModelParameters change                    */
    for (i=0; i<NR; i++) update |= (r[i] != rOld[i]);

    /* look-up or compute the current ordering state */
    if (update) {                                     /* if T, P, meltsAndCO2_H2OModelParameters or liquid composition change */
        double dgds[NT], sNew[NT], dgdsNORM=0.0, residual[NT];
        gsl_matrix *ptToD2gds2Copy = gsl_matrix_alloc((size_t) NT, (size_t) NT);
        gsl_vector_view vvToDgds = gsl_vector_view_array(dgds, (size_t) NT);
        for (i=0; i<NT; i++) { sNew[i] = sOld[i]; sOld[i] = 2.0; }
        convergedInOrder = TRUE;

        initialGuessOrdering(r, sNew);
        if (!rANDsTOx (r, sNew)) {
            printf("Initial guess to ordering state iteration is infeasible\n");
            doAbort = TRUE;
        }
        if (doAbort) {
            printf("Results of call to ordering with bad initial guess:\n");
            printf("   %20.20s %13.13s %13.13s %13.13s\n", "Species", "Mole frac", "r", "s");
            printf("   %20.20s %13.6g\n", liquid[0].label, xSpecies[0]);
            for (i=0;  i<NR; i++) printf("   %20.20s %13.6g %13.6g\n",         liquid[i+1].label,  xSpecies[i+1],  r[i]);
            for (i=0;  i<NS; i++) printf("   %20.20s %13.6g %13.13s %13.6g\n", liquid[i+NA].label, xSpecies[i+NA], "", sNew[i]);
            for (i=NS; i<NT; i++) printf("   %20.20s %13.13s %13.13s %13.6g\n", "order CN[*]", "", "", sNew[i]);
#ifdef USESEH
            doInterrupt = TRUE;
            raise_sigabrt(EXCEPTION_FLT_INVALID_OPERATION);
#else
            (void) raise(SIGABRT);
#endif
        }

        loop = TRUE;
        while (loop) {
            double s[NT], deltaS[NT], lambda=1.0, dgdsCopy[NT]; // dgdsCopy used as work space
            gsl_vector_view vvToDeltaS = gsl_vector_view_array(deltaS, (size_t) NT);
            gsl_vector_view vvToDgdsCopy = gsl_vector_view_array(dgdsCopy, (size_t) NT);
            int cycle = TRUE /*, errGaussj */;

            for (i=0; i<NT; i++) s[i] = sNew[i];

            fillDGDS   (r, s, t, p, dgds);
            fillD2GDS2New (r, s, t, p, d2gds2);

            //dgdsNORM=0.0
            for (i=0; i<NT; i++) {
                sOld[i] = s[i];
                for (j=0; j<NT; j++) if (s[i] == 0.0 || s[j] == 0.0) d2gds2[i][j] = 0.0;
                if (s[i] == 0.0) {
                    d2gds2[i][i] = 1.0;
                    dgds[i]      = 0.0;
                }
                //dgdsNORM += pow(dgds[i]/MAX(1.0, fabs(eosIntDGDS[i])), (double) 2.0);
            }
            //dgdsNORM = sqrt(dgdsNORM);

            gsl_matrix_scale(ptToD2gds2, -1.0);
            gsl_matrix_memcpy(ptToD2gds2Copy, ptToD2gds2);

            if (USE_SVD) {
                gsl_linalg_SV_decomp(ptToD2gds2, ptToV, ptToS, &vvToDgdsCopy.vector);
                gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToDgds.vector, &vvToDeltaS.vector);
            }
            else {
                melts_LU_decomp(ptToD2gds2, indexLU, &dLU);
                melts_LU_solve(ptToD2gds2, indexLU, &vvToDgds.vector, &vvToDeltaS.vector);
                melts_LU_refine(ptToD2gds2Copy, ptToD2gds2, indexLU,  &vvToDgds.vector, &vvToDeltaS.vector, &vvToDgdsCopy.vector);
            }
            for (i=0; i<NS; i++) s[i] += deltaS[i];


#ifdef DEBUG
            printf("--->dgds:   ");
            for (i=0; i<NT; i++) if (fabs(s[i]) > 10.0*DBL_EPSILON) printf("%20.13g", dgds[i]);
            printf("\n");
            for (j=0; j<NT; j++) if (fabs(s[j]) > 10.0*DBL_EPSILON) {
                printf("--->d2gds2[%d][]: ", j);
                for (i=0; i<NT; i++) if (fabs(s[i]) > 10.0*DBL_EPSILON) printf("%20.13g", d2gds2Copy[j][i]);
                printf("\n");
            }
            printf("--->s:      ");
            for (i=0; i<NT; i++) if (fabs(s[i]) > 10.0*DBL_EPSILON) printf("%20.13g", s[i]);
            printf("\n");
            printf("--->dels:   ");
            for (i=0; i<NT; i++) if (fabs(s[i]) > 10.0*DBL_EPSILON) printf("%20.13g", deltaS[i]);
            printf("\n");
#endif

            while (cycle && !rANDsTOx (r, s)) {
                lambda /= 2.0;
                for (j=0; j<NT; j++) s[j] = sOld[j] + lambda*deltaS[j];
                /* if (lambda < DBL_EPSILON) { */
                if (lambda < DBL_MIN) {
                    cycle = FALSE;
                    s[0] = (double) iter;
                    iter = MAX_ITER - 1;
                    fprintf(stderr, "\n*****lambda -> zero in ORDER. Terminating search loop.\n");
                  }
            }
#ifdef DEBUG
            printf("steplength correction:  = %20.13g\n", lambda);
            printf("--->s(adj): ");
            for (i=0; i<NT; i++) if (fabs(s[i]) > 10.0*DBL_EPSILON) printf("%20.13g", s[i]);
            printf("\n----- end -----\n");
#endif

            for (i=0; i<NT; i++) sNew[i] = s[i];
            iter++;
            loop = FALSE;
            if (iter < MAX_ITER) for (i=0; i<NT; i++) loop |= (fabs(sNew[i]-sOld[i]) > 10.0*DBL_EPSILON);
        }
        tOld = t;
        pOld = p;
        for (i=0; i<NR; i++) rOld[i] = r[i];

        (void) rANDsTOx (rOld, sOld);

        if (iter == MAX_ITER) {
            double sNorm;
            for (i=0, sNorm=0.0; i<NT; i++) sNorm += pow(sNew[i]-sOld[i], (double) 2.0);
            sNorm = sqrt(sNorm);
            if (sNorm > sqrt(DBL_EPSILON)) {
                /* convergedInOrder = FALSE; */
                fprintf(stderr, "ERROR in LIQUID.C (function ORDER). Failed to converge!\n");
                if (iter >= MAX_ITER) fprintf(stderr, " Iteration limit (%4d) exceeded.\n", iter);
                  fprintf(stderr, "   T (C) = %8.2f, P (GPa) = %10.4f\n", t-273.15, p/10000.0);
                fprintf(stderr, "   %20.20s %13.13s %13.13s %13.13s %13.13s %13.13s\n", "Species", "Mole frac", "r", "s", "dgds", "deltaS");
                fprintf(stderr, "   %20.20s %13.6g\n", liquid[0].label, xSpecies[0]);
                for (i=0;  i<NR; i++) fprintf(stderr, "   %20.20s %13.6g %13.6g\n",         liquid[i+1].label,  xSpecies[i+1],  r[i]);
                for (i=0;  i<NS; i++) fprintf(stderr, "   %20.20s %13.6g %13.13s %13.6g %13.6g %13.6g\n", liquid[i+NA].label, xSpecies[i+NA], "", sOld[i], dgds[i], sNew[i]-sOld[i]);
                for (i=NS; i<NT; i++) fprintf(stderr, "   %20.20s %13.13s %13.13s %13.6g %13.6g %13.6g\n", "order CN[*]", "", "", sOld[i], dgds[i], sNew[i]-sOld[i]);
                fprintf(stderr, " sNorm             = %20.13g\n", sNorm);
                fprintf(stderr, " dgdsNorm          = %20.13g\n", dgdsNORM);
                fprintf(stderr, " 10*DBL_EPSILON    = %20.13g\n", 10.0*DBL_EPSILON);
                fprintf(stderr, " DBL_EPSILON^(2/3) = %20.13g\n", pow(DBL_EPSILON, 2.0/3.0));
                fprintf(stderr, " DBL_EPSILON^(1/2) = %20.13g\n", sqrt(DBL_EPSILON));
            } else if (sNorm > pow(DBL_EPSILON, 2.0/3.0)) {
                fprintf(stderr, "WARNING in LIQUID.C (function ORDER). sNorm = %g, dgdsNorm = %g [eps = %g, sqrt(eps) = %g]\n", sNorm, dgdsNORM, DBL_EPSILON, sqrt(DBL_EPSILON));
            }
        }

#ifdef DEBUG
        printf("Results of ordering state calculation:\n");
        printf("   T (C) = %8.2f, P (GPa) = %10.4f\n", t-273.15, p/10000.0);
        printf("   %20.20s %13.13s %13.13s %13.13s %13.13s %13.13s %13.13s\n", "Species", "Mole frac", "r", "s", "dgds", "deltaS", "eosIntDGDS");
        printf("   %20.20s %13.6g\n", liquid[0].label, xSpecies[0]);
        for (i=0;  i<NR; i++) printf("   %20.20s %13.6g %13.6g\n", liquid[i+1].label, xSpecies[i+1], r[i]);
        for (i=0;  i<NS; i++) printf("   %20.20s %13.6g %13.13s %13.6g %13.6g %13.6g %13.6g\n", liquid[i+NA].label, xSpecies[i+NA], "", sOld[i], dgds[i], sNew[i]-sOld[i], eosIntDGDS[i]);
        for (i=NS; i<NT; i++) printf("   %20.20s %13.13s %13.13s %13.6g %13.6g %13.6g\n", "order CN[*]", "", "", sOld[i], dgds[i], sNew[i]-sOld[i]);
        printf(" 10*DBL_EPSILON    = %20.13g\n", 10.0*DBL_EPSILON);
        printf(" DBL_EPSILON^(2/3) = %20.13g\n", pow(DBL_EPSILON, 2.0/3.0));
        printf(" DBL_EPSILON^(1/2) = %20.13g\n", sqrt(DBL_EPSILON));
        printf(" eosIntegralBranch = %s\n", (eosIntegralBranch == GMAPeosBRANCH) ? "GMAP" : "LMAP");
        for (i=0; i<NS; i++) {
            double s[NT];
            for (j=0; j<NT; j++) s[j] = sOld[j];
            s[i] += 30.0*DBL_EPSILON;
            fillDGDS(r, s, t, p, dgds);
            printf(" s[%d]+30eps gives dgds = %20.13g\n", i, dgds[i]);
            s[i] -= 10.0*DBL_EPSILON;
            fillDGDS(r, s, t, p, dgds);
            printf(" s[%d]+20eps gives dgds = %20.13g\n", i, dgds[i]);
            s[i] -= 10.0*DBL_EPSILON;
            fillDGDS(r, s, t, p, dgds);
            printf(" s[%d]+10eps gives dgds = %20.13g\n", i, dgds[i]);
            printf(" ----------------------\n");
            s[i] = sOld[i];
            fillDGDS(r, s, t, p, dgds);
            printf(" s[%d]       gives dgds = %20.13g\n", i, dgds[i]);
            printf(" ----------------------\n");
            s[i] -= 10.0*DBL_EPSILON;
            fillDGDS(r, s, t, p, dgds);
            printf(" s[%d]-10eps gives dgds = %20.13g\n", i, dgds[i]);
            s[i] -= 10.0*DBL_EPSILON;
            fillDGDS(r, s, t, p, dgds);
            printf(" s[%d]-20eps gives dgds = %20.13g\n", i, dgds[i]);
            s[i] -= 10.0*DBL_EPSILON;
            fillDGDS(r, s, t, p, dgds);
            printf(" s[%d]-30eps gives dgds = %20.13g\n", i, dgds[i]);
            printf(" ++++++++++++++++++++++\n");
        }
#endif

        gsl_matrix_free(ptToD2gds2Copy);

        setTOld(tOld);
        setPOld(pOld);
        /* arrays (rOld, sOld, d2gds2, etc) should be preserved automatically */

    }

    if (mask & FIRST  ) {   /* return s        */
        for (i=0; i<NT; i++) s[i] = sOld[i];
    }
    if (mask & SECOND ) {   /* compute ds/dr:  */
        double *s = sOld;
        double d2gdrds[NR][NT];
        gsl_matrix_view mvToDr = gsl_matrix_view_array((double *) dr, (size_t) NT, (size_t) NR),
            mvToD2gdrds = gsl_matrix_view_array((double *) d2gdrds, (size_t) NR, (size_t) NT);

        fillD2GDRDS (r, s, t, p,  d2gdrds);

        for (j=0; j<NR; j++) {
            gsl_vector_view vvToDr = gsl_matrix_column(&mvToDr.matrix, j),
                vvToD2gdrds = gsl_matrix_row(&mvToD2gdrds.matrix, j);
            if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdrds.vector, &vvToDr.vector);
            else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdrds.vector, &vvToDr.vector);
            for (i=0; i<NT; i++) dr[i][j] = (s[i] > 0.0) ? dr[i][j] : 0.0;
        }
    }

    if (mask & THIRD  ) {   /* compute ds/dt:  */
        double *s = sOld;
        double d2gdsdt[NT];
        gsl_vector_view vvToDt = gsl_vector_view_array(dt, (size_t) NT),
            vvToD2gdsdt = gsl_vector_view_array(d2gdsdt, (size_t) NT);

        fillD2GDSDT (r, s, t, p, d2gdsdt);

        if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdt.vector, &vvToDt.vector);
        else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdt.vector, &vvToDt.vector);
        for (i=0; i<NT; i++) dt[i] = (s[i] > 0.0) ? dt[i] : 0.0;

    }
    if (mask & FOURTH ) {   /* compute ds/dp:  */
        double *s = sOld;
        double d2gdsdp[NT];
        gsl_vector_view vvToDp = gsl_vector_view_array(dp, (size_t) NT),
            vvToD2gdsdp = gsl_vector_view_array(d2gdsdp, (size_t) NT);

        fillD2GDSDP (r, s, t, p, d2gdsdp);

        if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdp.vector, &vvToDp.vector);
        else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdp.vector, &vvToDp.vector);
        for (i=0; i<NT; i++) dp[i] = (s[i] > 0.0) ? dp[i] : 0.0;
    }
    if (mask & FIFTH  ) {   /* compute d2s/dr2 */
        double *s = sOld;
        double d2gdrds[NR][NT];
        double d3gdr2ds[NR][NR][NT];
        double d3gdrds2[NR][NT][NT];
        double d3gds3[NT][NT][NT];
        double dsdr[NT][NR], temp[NT], temp2[NT];
        gsl_matrix_view mvToDsdr = gsl_matrix_view_array((double *) dsdr, (size_t) NT, (size_t) NR),
            mvToD2gdrds = gsl_matrix_view_array((double *) d2gdrds, (size_t) NR, (size_t) NT);
        gsl_vector_view vvToTemp = gsl_vector_view_array(temp, (size_t) NT), vvToTemp2 = gsl_vector_view_array(temp2, (size_t) NT);
        int k, l, m, n;

        fillD2GDRDS  (r, s, t, p,  d2gdrds);
        fillD3GDR2DS (r, s, t, p,  d3gdr2ds);
        fillD3GDRDS2 (r, s, t, p,  d3gdrds2);
        fillD3GDS3   (r, s, t, p,  d3gds3);

        /* compute dsdr matrix */
        for (j=0; j<NR; j++) {
            gsl_vector_view vvToDsdr = gsl_matrix_column(&mvToDsdr.matrix, j),
              vvToD2gdrds = gsl_matrix_row(&mvToD2gdrds.matrix, j);
            if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdrds.vector, &vvToDsdr.vector);
            else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdrds.vector, &vvToDsdr.vector);
            for (i=0; i<NT; i++) dsdr[i][j] = (s[i] > 0.0) ? dsdr[i][j] : 0.0;
        }

        /* compute dsdr2 cube */
        for (j=0; j<NR; j++) {
            for (k=0; k<NR; k++) {
                for (l=0; l<NT; l++) {
                    temp[l] = d3gdr2ds[j][k][l];
                    for (m=0; m<NT; m++) {
                        temp[l] += d3gdrds2[j][l][m]*dsdr[m][k]
                            + d3gdrds2[k][l][m]*dsdr[m][j];
                        for (n=0; n<NT; n++)
                            temp[l] += d3gds3[l][m][n]*dsdr[m][j]*dsdr[n][k];
                    }
                }
                if (USE_SVD) {
                    gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToTemp.vector, &vvToTemp2.vector);
                    for (i=0; i<NT; i++) dr2[i][j][k] = (s[i] > 0.0) ? temp2[i] : 0.0;
                }
                else {
                    melts_LU_svx(ptToD2gds2, indexLU, &vvToTemp.vector);
                    for (i=0; i<NT; i++) dr2[i][j][k] = (s[i] > 0.0) ? temp[i] : 0.0;
                }
            }
        }
    }

    if (mask & SIXTH  ) {   /* compute d2s/drt */
        double *s = sOld;
        double d2gdrds[NR][NT];
        double d2gdsdt[NT];
        double d3gdrds2[NR][NT][NT];
        double d3gdrdsdt[NR][NT];
        double d3gds2dt[NT][NT];
        double d3gds3[NT][NT][NT];
        double dsdr[NT][NR], dsdt[NT], temp[NT], temp2[NT];
        gsl_matrix_view mvToDsdr = gsl_matrix_view_array((double *) dsdr, (size_t) NT, (size_t) NR),
            mvToD2gdrds = gsl_matrix_view_array((double *) d2gdrds, (size_t) NR, (size_t) NT);
        gsl_vector_view vvToDsdt = gsl_vector_view_array(dsdt, (size_t) NT),
            vvToD2gdsdt = gsl_vector_view_array(d2gdsdt, (size_t) NT),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NT), vvToTemp2 = gsl_vector_view_array(temp2, (size_t) NT);
        int k, l, m;

        fillD2GDRDS   (r, s, t, p,  d2gdrds);
        fillD2GDSDT   (r, s, t, p,  d2gdsdt);
        fillD3GDRDS2  (r, s, t, p,  d3gdrds2);
        fillD3GDRDSDT (r, s, t, p,  d3gdrdsdt);
        fillD3GDS3    (r, s, t, p,  d3gds3);
        fillD3GDS2DT  (r, s, t, p,  d3gds2dt);

        /* compute dsdr matrix */
        for (j=0; j<NR; j++) {
            gsl_vector_view vvToDsdr = gsl_matrix_column(&mvToDsdr.matrix, j),
              vvToD2gdrds = gsl_matrix_row(&mvToD2gdrds.matrix, j);
            if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdrds.vector, &vvToDsdr.vector);
            else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdrds.vector, &vvToDsdr.vector);
            for (i=0; i<NT; i++) dsdr[i][j] = (s[i] > 0.0) ? dsdr[i][j] : 0.0;
        }

        /* compute dsdt vector */
        if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdt.vector, &vvToDsdt.vector);
        else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdt.vector, &vvToDsdt.vector);
        for (i=0; i<NT; i++) dsdt[i] = (s[i] > 0.0) ? dsdt[i] : 0.0;

        /* compute dsdrdt matrix */
        for (j=0; j<NR; j++) {
            for (k=0; k<NT; k++) {
                temp[k] = d3gdrdsdt[j][k];
                for (l=0; l<NT; l++) {
                    temp[k] += d3gdrds2[j][k][l]*dsdt[l] + d3gds2dt[k][l]*dsdr[l][j];
                    for (m=0; m<NT; m++) temp[k] += d3gds3[k][l][m]*dsdr[l][j]*dsdt[m];
                }
            }
            if (USE_SVD) {
                gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToTemp.vector, &vvToTemp2.vector);
                for (i=0; i<NT; i++) drt[i][j] = (s[i] > 0.0) ? temp2[i] : 0.0;
            }
            else {
                melts_LU_svx(ptToD2gds2, indexLU, &vvToTemp.vector);
                for (i=0; i<NT; i++) drt[i][j] = (s[i] > 0.0) ? temp[i] : 0.0;
            }
        }

    }
    if (mask & SEVENTH) {   /* compute d2s/drp */
        double *s = sOld;
        double d2gdrds[NR][NT];
        double d2gdsdp[NT];
        double d3gdrds2[NR][NT][NT];
        double d3gdrdsdp[NR][NT];
        double d3gds3[NT][NT][NT];
        double d3gds2dp[NT][NT];
        double dsdr[NT][NR], dsdp[NT], temp[NT], temp2[NT];
        gsl_matrix_view mvToDsdr = gsl_matrix_view_array((double *) dsdr, (size_t) NT, (size_t) NR),
            mvToD2gdrds = gsl_matrix_view_array((double *) d2gdrds, (size_t) NR, (size_t) NT);
        gsl_vector_view vvToDsdp = gsl_vector_view_array(dsdp, (size_t) NT),
            vvToD2gdsdp = gsl_vector_view_array(d2gdsdp, (size_t) NT),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NT), vvToTemp2 = gsl_vector_view_array(temp2, (size_t) NT);

        int k, l, m;

        fillD2GDRDS   (r, s, t, p,  d2gdrds);
        fillD2GDSDP   (r, s, t, p,  d2gdsdp);
        fillD3GDRDS2  (r, s, t, p,  d3gdrds2);
        fillD3GDRDSDP (r, s, t, p,  d3gdrdsdp);
        fillD3GDS3    (r, s, t, p,  d3gds3);
        fillD3GDS2DP  (r, s, t, p,  d3gds2dp);

        /* compute dsdr matrix */
        for (j=0; j<NR; j++) {
            gsl_vector_view vvToDsdr = gsl_matrix_column(&mvToDsdr.matrix, j),
                vvToD2gdrds = gsl_matrix_row(&mvToD2gdrds.matrix, j);
            if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdrds.vector, &vvToDsdr.vector);
            else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdrds.vector, &vvToDsdr.vector);
            for (i=0; i<NT; i++) dsdr[i][j] = (s[i] > 0.0) ? dsdr[i][j] : 0.0;
        }

        /* compute dsdp vector */
        if (USE_SVD)  gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdp.vector, &vvToDsdp.vector);
        else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdp.vector, &vvToDsdp.vector);
        for (i=0; i<NT; i++) dsdp[i] = (s[i] > 0.0) ? dsdp[i] : 0.0;

        /* compute dsdrdp matrix */
        for (j=0; j<NR; j++) {
            for (k=0; k<NT; k++) {
                temp[k] = d3gdrdsdp[j][k];
                for (l=0; l<NT; l++) {
                        temp[k] += d3gdrds2[j][k][l]*dsdp[l] + d3gds2dp[k][l]*dsdr[l][j];
                        for (m=0; m<NT; m++) temp[k] += d3gds3[k][l][m]*dsdr[l][j]*dsdp[m];
                }
            }
            if (USE_SVD) {
                gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToTemp.vector, &vvToTemp2.vector);
                for (i=0; i<NT; i++) drp[i][j] = (s[i] > 0.0) ? temp2[i] : 0.0;
            }
            else {
                melts_LU_svx(ptToD2gds2, indexLU, &vvToTemp.vector);
                for (i=0; i<NT; i++) drp[i][j] = (s[i] > 0.0) ? temp[i] : 0.0;
            }
        }
    }

    if (mask & EIGHTH ) {   /* compute d2s/dt2 */
        double *s = sOld;
        double d2gdsdt[NT];
        double d3gds3[NT][NT][NT];
        double d3gds2dt[NT][NT];
        double d3gdsdt2[NT];
        double dsdt[NT], temp[NT], temp2[NT];
        gsl_vector_view vvToDsdt = gsl_vector_view_array(dsdt, (size_t) NT),
            vvToD2gdsdt = gsl_vector_view_array(d2gdsdt, (size_t) NT),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NT), vvToTemp2 = gsl_vector_view_array(temp2, (size_t) NT);
        int k, l;

        fillD2GDSDT  (r, s, t, p,  d2gdsdt);
        fillD3GDS3   (r, s, t, p,  d3gds3);
        fillD3GDS2DT (r, s, t, p,  d3gds2dt);
        fillD3GDSDT2 (r, s, t, p,  d3gdsdt2);

        /* compute dsdt vector */
        if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdt.vector, &vvToDsdt.vector);
        else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdt.vector, &vvToDsdt.vector);
        for (i=0; i<NT; i++) dsdt[i] = (s[i] > 0.0) ? dsdt[i] : 0.0;

        /* compute dsdt2 vector */
        for (j=0; j<NT; j++) {
            temp[j] = d3gdsdt2[j];
            for (k=0; k<NT; k++) {
                temp[j] +=  2.0*d3gds2dt[j][k]*dsdt[k];
                for (l=0; l<NT; l++) temp[j] += d3gds3[j][k][l]*dsdt[k]*dsdt[l];
            }
        }
        if (USE_SVD) {
            gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToTemp.vector, &vvToTemp2.vector);
            for (i=0; i<NT; i++) dt2[i] = (s[i] > 0.0) ? temp2[i] : 0.0;
        }
        else {
            melts_LU_svx(ptToD2gds2, indexLU, &vvToTemp.vector);
            for (i=0; i<NT; i++) dt2[i] = (s[i] > 0.0) ? temp[i] : 0.0;
        }
    }

    if (mask & NINTH  ) {   /* compute d2s/dtp */
        double *s = sOld;
        double d2gdsdt[NT];
        double d2gdsdp[NT];
        double d3gds3[NT][NT][NT];
        double d3gds2dt[NT][NT];
        double d3gds2dp[NT][NT];
        double d3gdsdtdp[NT];
        double dsdt[NT], dsdp[NT], temp[NT], temp2[NT];
        gsl_vector_view vvToDsdt = gsl_vector_view_array(dsdt, (size_t) NT),
            vvToD2gdsdt = gsl_vector_view_array(d2gdsdt, (size_t) NT),
            vvToDsdp = gsl_vector_view_array(dsdp, (size_t) NT),
            vvToD2gdsdp = gsl_vector_view_array(d2gdsdp, (size_t) NT),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NT), vvToTemp2 = gsl_vector_view_array(temp2, (size_t) NT);
        int k, l;

        fillD2GDSDT   (r, s, t, p,  d2gdsdt);
        fillD2GDSDP   (r, s, t, p,  d2gdsdp);
        fillD3GDS3    (r, s, t, p,  d3gds3);
        fillD3GDS2DT  (r, s, t, p,  d3gds2dt);
        fillD3GDS2DP  (r, s, t, p,  d3gds2dp);
        fillD3GDSDTDP (r, s, t, p,  d3gdsdtdp);

        /* compute dsdt vector */
        if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdt.vector, &vvToDsdt.vector);
        else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdt.vector, &vvToDsdt.vector);
        for (i=0; i<NT; i++) dsdt[i] = (s[i] > 0.0) ? dsdt[i] : 0.0;

        /* compute dsdp vector */
        if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdp.vector, &vvToDsdp.vector);
        else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdp.vector, &vvToDsdp.vector);
        for (i=0; i<NT; i++) dsdp[i] = (s[i] > 0.0) ? dsdp[i] : 0.0;

        /* compute dsdtp vector */
        for (j=0; j<NT; j++) {
            temp[j] = d3gdsdtdp[j];
            for (k=0; k<NT; k++) {
                temp[j] += d3gds2dt[j][k]*dsdp[k] + d3gds2dp[j][k]*dsdt[k];
                for (l=0; l<NT; l++) temp[j] += d3gds3[j][k][l]*dsdt[k]*dsdp[l];
            }
        }
        if (USE_SVD) {
            gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToTemp.vector, &vvToTemp2.vector);
            for (i=0; i<NT; i++) dtp[i] = (s[i] > 0.0) ? temp2[i] : 0.0;
        }
        else {
            melts_LU_svx(ptToD2gds2, indexLU, &vvToTemp.vector);
            for (i=0; i<NT; i++) dtp[i] = (s[i] > 0.0) ? temp[i] : 0.0;
        }
    }

    if (mask & TENTH  ) {   /* compute d2s/dp2 */
        double *s = sOld;
        double d2gdsdp[NT];
        double d3gds3[NT][NT][NT];
        double d3gds2dp[NT][NT];
        double d3gdsdp2[NT];
        double dsdp[NT], temp[NT], temp2[NT];
        gsl_vector_view vvToDsdp = gsl_vector_view_array(dsdp, (size_t) NT),
            vvToD2gdsdp = gsl_vector_view_array(d2gdsdp, (size_t) NT),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NT), vvToTemp2 = gsl_vector_view_array(temp2, (size_t) NT);
        int k, l;

        fillD2GDSDP  (r, s, t, p,  d2gdsdp);
        fillD3GDS3   (r, s, t, p,  d3gds3);
        fillD3GDS2DP (r, s, t, p,  d3gds2dp);
        fillD3GDSDP2 (r, s, t, p,  d3gdsdp2);

        /* compute dsdp vector */
        if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdp.vector, &vvToDsdp.vector);
        else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdp.vector, &vvToDsdp.vector);
        for (i=0; i<NT; i++) dsdp[i] = (s[i] > 0.0) ? dsdp[i] : 0.0;

        /* compute dsdp2 vector */
        for (j=0; j<NT; j++) {
            temp[j] = d3gdsdp2[j];
            for (k=0; k<NT; k++) {
                temp[j] +=  2.0*d3gds2dp[j][k]*dsdp[k];
                for (l=0; l<NT; l++) temp[j] += d3gds3[j][k][l]*dsdp[k]*dsdp[l];
            }
        }
        if (USE_SVD) {
            gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToTemp.vector, &vvToTemp2.vector);
            for (i=0; i<NT; i++) dp2[i] = (s[i] > 0.0) ? temp2[i] : 0.0;
        }
        else {
            melts_LU_svx(ptToD2gds2, indexLU, &vvToTemp.vector);
            for (i=0; i<NT; i++) dp2[i] = (s[i] > 0.0) ? temp[i] : 0.0;
        }
    }

    if (mask & ELEVENTH ) {   /* compute ds/dw:  */
        double *s = sOld;
        double d2gdsdw[NT][3*NP];
        gsl_matrix_view mvToDw = gsl_matrix_view_array((double *) dw, (size_t) NT, (size_t) 3*NP),
            mvToD2gdsdw = gsl_matrix_view_array((double *) d2gdsdw, (size_t) NT, (size_t) 3*NP);

        fillD2GDSDW (r, s, t, p,  d2gdsdw);
        for (j=0; j<NP; j++) {
            gsl_vector_view vvToDw, vvToD2gdsdw;
            vvToDw = gsl_matrix_column(&mvToDw.matrix,      j);
            vvToD2gdsdw = gsl_matrix_column(&mvToD2gdsdw.matrix,      j);
            if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdw.vector, &vvToDw.vector);
            else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdw.vector, &vvToDw.vector);
            for (i=0; i<NT; i++) dw[i][     j] = (s[i] > 0.0) ? dw[i][     j] : 0.0;
            vvToDw = gsl_matrix_column(&mvToDw.matrix,   NP+j);
            vvToD2gdsdw = gsl_matrix_column(&mvToD2gdsdw.matrix,   NP+j);
            if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdw.vector, &vvToDw.vector);
            else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdw.vector, &vvToDw.vector);
            for (i=0; i<NT; i++) dw[i][  NP+j] = (s[i] > 0.0) ? dw[i][  NP+j] : 0.0;
            vvToDw = gsl_matrix_column(&mvToDw.matrix, 2*NP+j);
            vvToD2gdsdw = gsl_matrix_column(&mvToD2gdsdw.matrix, 2*NP+j);
            if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdw.vector, &vvToDw.vector);
            else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdw.vector, &vvToDw.vector);
            for (i=0; i<NT; i++) dw[i][2*NP+j] = (s[i] > 0.0) ? dw[i][2*NP+j] : 0.0;
        }
    }

    if (mask & TWELFTH ) {   /* compute ds/dtdw:  */
        double *s = sOld;
        double d2gdsdt[NT];
        double d2gdsdw[NT][3*NP];
        double d3gds3[NT][NT][NT];
        double d3gds2dt[NT][NT];
        double d3gds2dw[NT][NT][3*NP];
        double d3gdsdtdw[NT][3*NP];
        double dsdt[NT], dsdw[NT][3*NP], temp[NT], temp2[NT];
        gsl_matrix_view mvToDsdw = gsl_matrix_view_array((double *) dsdw, (size_t) NT, (size_t) 3*NP),
            mvToD2gdsdw = gsl_matrix_view_array((double *) d2gdsdw, (size_t) NT, (size_t) 3*NP);
        gsl_vector_view vvToDsdt = gsl_vector_view_array(dsdt, (size_t) NT),
            vvToD2gdsdt = gsl_vector_view_array(d2gdsdt, (size_t) NT),
            vvToTemp = gsl_vector_view_array(temp, (size_t) NT), vvToTemp2 = gsl_vector_view_array(temp2, (size_t) NT);

        fillD2GDSDW (r, s, t, p,  d2gdsdw);
        for (j=0; j<NP; j++) {
            gsl_vector_view vvToDsdw, vvToD2gdsdw;
            vvToDsdw = gsl_matrix_column(&mvToDsdw.matrix,      j);
            vvToD2gdsdw = gsl_matrix_column(&mvToD2gdsdw.matrix,      j);
            if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdw.vector, &vvToDsdw.vector);
            else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdw.vector, &vvToDsdw.vector);
            for (i=0; i<NT; i++) dsdw[i][     j] = (s[i] > 0.0) ? dsdw[i][     j] : 0.0;
            vvToDsdw = gsl_matrix_column(&mvToDsdw.matrix,   NP+j);
            vvToD2gdsdw = gsl_matrix_column(&mvToD2gdsdw.matrix,   NP+j);
            if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdw.vector, &vvToDsdw.vector);
            else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdw.vector, &vvToDsdw.vector);
            for (i=0; i<NT; i++) dsdw[i][  NP+j] = (s[i] > 0.0) ? dsdw[i][  NP+j] : 0.0;
            vvToDsdw = gsl_matrix_column(&mvToDsdw.matrix, 2*NP+j);
            vvToD2gdsdw = gsl_matrix_column(&mvToD2gdsdw.matrix, 2*NP+j);
            if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdw.vector, &vvToDsdw.vector);
            else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdw.vector, &vvToDsdw.vector);
            for (i=0; i<NT; i++) dsdw[i][2*NP+j] = (s[i] > 0.0) ? dsdw[i][2*NP+j] : 0.0;
        }

        fillD2GDSDT (r, s, t, p, d2gdsdt);
        if (USE_SVD) gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToD2gdsdt.vector, &vvToDsdt.vector);
        else melts_LU_solve(ptToD2gds2, indexLU, &vvToD2gdsdt.vector, &vvToDsdt.vector);
        for (i=0; i<NT; i++) dsdt[i] = (s[i] > 0.0) ? dsdt[i] : 0.0;

        fillD3GDSDTDW (r, s, t, p,  d3gdsdtdw);
        fillD3GDS2DW  (r, s, t, p,  d3gds2dw);
        fillD3GDS2DT  (r, s, t, p,  d3gds2dt);
        fillD3GDS3    (r, s, t, p,  d3gds3);
        for (j=0; j<NP; j++) {
            for (i=0; i<NT; i++) {
                int k, l;
                temp[i] = d3gdsdtdw[i][j];
                for (k=0; k<NT; k++) {
                    if (dsdw[k][j] != 0.0) temp[i] += d3gds2dt[i][k]*dsdw[k][j];
                    if (dsdt[k] != 0.0) {
                        temp[i] += d3gds2dw[i][k][j]*dsdt[k];
                        for (l=0; l<NT; l++) if (dsdw[l][j] != 0.0) temp[i] += d3gds3[i][k][l]*dsdt[k]*dsdw[l][j];
                    }
                }
            }
            if (USE_SVD) {
                gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToTemp.vector, &vvToTemp2.vector);
                for (i=0; i<NT; i++) dtw[i][j] = (s[i] > 0.0) ? temp2[i] : 0.0;
            }
            else {
                melts_LU_svx(ptToD2gds2, indexLU, &vvToTemp.vector);
                for (i=0; i<NT; i++) dtw[i][j] = (s[i] > 0.0) ? temp[i] : 0.0;
            }
            for (i=0; i<NT; i++) {
                int k, l;
                temp[i] = d3gdsdtdw[i][j+NP];
                for (k=0; k<NT; k++) {
                    if (dsdw[k][j+NP] != 0.0) temp[i] += d3gds2dt[i][k]*dsdw[k][j+NP];
                    if (dsdt[k] != 0.0) {
                        temp[i] += d3gds2dw[i][k][j+NP]*dsdt[k];
                        for (l=0; l<NT; l++) if (dsdw[l][j+NP] != 0.0) temp[i] += d3gds3[i][k][l]*dsdt[k]*dsdw[l][j+NP];
                    }
                }
            }
            if (USE_SVD) {
                gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToTemp.vector, &vvToTemp2.vector);
                for (i=0; i<NT; i++) dtw[i][j+NP] = (s[i] > 0.0) ? temp2[i] : 0.0;
            }
            else {
                melts_LU_svx(ptToD2gds2, indexLU, &vvToTemp.vector);
                for (i=0; i<NT; i++) dtw[i][j+NP] = (s[i] > 0.0) ? temp[i] : 0.0;
            }
            for (i=0; i<NT; i++) {
                int k, l;
                temp[i] = d3gdsdtdw[i][j+2*NP];
                for (k=0; k<NT; k++) {
                    if (dsdw[k][j+2*NP] != 0.0) temp[i] += d3gds2dt[i][k]*dsdw[k][j+2*NP];
                    if (dsdt[k] != 0.0) {
                        temp[i] += d3gds2dw[i][k][j+2*NP]*dsdt[k];
                        for (l=0; l<NT; l++) if (dsdw[l][j+2*NP] != 0.0) temp[i] += d3gds3[i][k][l]*dsdt[k]*dsdw[l][j+2*NP];
                    }
                }
            }
            if (USE_SVD) {
                gsl_linalg_SV_solve(ptToD2gds2, ptToV, ptToS, &vvToTemp.vector, &vvToTemp2.vector);
                for (i=0; i<NT; i++) dtw[i][j+2*NP] = (s[i] > 0.0) ? temp2[i] : 0.0;
            }
            else {
                melts_LU_svx(ptToD2gds2, indexLU, &vvToTemp.vector);
                for (i=0; i<NT; i++) dtw[i][j+2*NP] = (s[i] > 0.0) ? temp[i] : 0.0;
            }
        }
    }

}

/*
 *=============================================================================
 * Public functions:
 *    inpMask -  bitwise mask for specifying input parameters
 *    outMask -  bitwise mask for selecting output
 *    mask    -  bitwise mask for selecting output
 *    t       -  Temperature (K)
 *    p       -  Pressure (bars)
 *    *r      -  (pointer to x[]) Array of independent compositional variables
 */

void
conLiq_CO2_H2O(int inpMask, int outMask, double t, double p,
    double *o,      /* comp of liquid in moles of oxides                        */
    double *m,      /* comp of liquid in moles of endmember components          */
    double *r,      /* comp of liquid in terms of the independent comp var      */
    double *x,      /* comp of liquid in mole fractions of endmember components */
    double **dm,    /* Jacobian matrix: dm[i][j] = dr[i]/dm[j]                  */
    double ***d2m,  /* vector of matrices: d2m[i][j][k] = d2r[i]/dm[j]dm[k]     */
    double *logfo2) /* base 10 logarithm of the oxygen fugacity                 */
{
    /*---------------------------------------------------------------------------
    Not all combinations of inpMask and outMask are feasible. Valid
        combinations are:

       inpMask          outMask
    (1)  FIRST | SEVENTH  FIRST
    (2)  FIRST            SEVENTH
    (3)  SECOND           THIRD | FOURTH | FIFTH | SIXTH
    (4)  THIRD            FOURTH
    (5)  THIRD            FOURTH | EIGHTH
    (6)  THIRD            FOURTH | NINTH

    (1) converts a vector of moles of oxides into a vector of moles of oxides
            with the correct redox state for the given t, p, and logfo2. Note that
            the original vector is used as output.
    (2) calculates from a vector of moles of oxides and the given t and p, the
            appropriate logfo2
    (3) calculates from a vector of moles of endmember components, one or
            all of: r[], x[], dr[]/dm[], or d2r[]/dm[]dm[]
    (4) calculates from a vector of independent compositional variables
            mole fractions of endmember components
    (5) calculates from a vector of independent compositional variables
            mole fractions of endmember species
    (6) calculates from a vector of independent compositional variables
            fraction of system in CN state (vector 0 ... NT-1)
    ----------------------------------------------------------------------------*/

    int i, j, k;

    MTHREAD_ONCE(&initThreadBlock, threadInit);

    if (inpMask & FIRST) {
        /*-------------------------------------------------------------------------
            Oxide+logfo2 -> oxide or oxide -> logfo2. The algorithm used is that
            given by:
            Kress, VC and Carmichael, ISE (1991) The compressibility of silicate
                liquids containing Fe2O3 and the effect of composition, temperature,
                oxygen fugacity and pressure on their redox states.
                Contributions to Mineralogy and Petrology (in press)
            Coefficients for the oxides are initialized in LIQ_STRUCT_DATA.H
        --------------------------------------------------------------------------*/
        static int indexFeO = 0, indexFe2O3 = 0;
#ifdef USE_KRESS_CARMICHAEL_FO2
        static const double t0 = 1673.15,                      /* K       */
                         a =    0.196,
                         b =    1.1492e4,                  /* K       */
                         c =   -6.675,
                         e =   -3.364,
                         f =   -7.01e-7  * 1.0e5,          /* K/bar   */
                         g =   -1.54e-10 * 1.0e5,          /* 1/bar   */
                         h =    3.85e-17 * 1.0e5 * 1.0e5;  /* K/bar^2 */
#endif /* USE_KRESS_CARMICHAEL_FO2 */
        double sum = 0.0, temp;

        if (indexFeO == 0 || indexFe2O3 == 0) {
            for (i=0; i<nc; i++) {
                if (bulkSystem[i].type == FEO)   indexFeO   = i;
                if (bulkSystem[i].type == FE2O3) indexFe2O3 = i;
            }
            if (indexFeO == 0 || indexFe2O3 == 0) {
         printf("Fatal error in conLiq_CO2_H2O (LIQUID.C)\n");
         printf("The oxides FeO and Fe2O3 cannot be identified.\n");
         return;
            }
        }

        if (inpMask == (FIRST | SEVENTH)  && outMask == FIRST) {
            /*----------------------------------------------------------------------
                Converts a vector of moles of oxides (as defined in LIQ_STRUCT_DATA.H
                for the structure bulkSystem) into a vector of moles of oxides with
                the correct redox state (ferric/ferrous ratio) for the given bulk
                composition, t and p.
            ------------------------------------------------------------------------*/

            o[indexFeO]   += 2.0*o[indexFe2O3];
            o[indexFe2O3]  = 0.0;
            if (o[indexFeO] == 0.0) return;

            for (i=0; i<nc; i++) sum += o[i];
            if (sum == 0.0) return;

#ifdef USE_KRESS_CARMICHAEL_FO2
            if (p < 50000.0) temp = a*log(10.0)*(*logfo2) + b/t + c + e*(1.0 - t0/t - log(t/t0))
                                                        + f*p/t + g*(t-t0)*p/t + h*SQUARE(p)/t;
            else             temp = a*log(10.0)*(*logfo2) + b/t + c + e*(1.0 - t0/t - log(t/t0))
                                                        + f*50000.0/t + g*(t-t0)*50000.0/t + h*SQUARE(50000.0)/t
                - a*log(10.0)*(608.966*p/10000.0-608.966*5.0)/t;

            for (i=0; i<nc; i++) temp += bulkSystem[i].coeff*o[i]/sum;
            temp = exp(temp);

            o[indexFe2O3]  = temp*o[indexFeO]/(1.0 + 2.0*temp);
            o[indexFeO]   -= 2.0*o[indexFe2O3];
#else
            {
                double y = 0.3;
                double intVFe2O3  = integralV_GKsp(iOxFe2O3,  t, p);
                double intVFeO    = integralV_GKsp(iOxFeO,    t, p);
                double intVFeO1_3 = integralV_GKsp(iOxFeO1_3, t, p);
                double deltaG = -106200.0 - t*(-55.1) + 31.86*(t - 1673.15 - t*log(t/1673.15)) + intVFe2O3/2.0 - intVFeO;
                double KD1 = exp(-deltaG/(R*t)
                            -(39860.0*o[iOxAl2O3] - 62520.0*o[iOxCaO] - 102000.0*o[iOxNa2O] - 119000.0*o[iOxK2O])/(R*t*sum));
                double K2 = 0.4*exp(-(intVFeO1_3 - (1.0-2.0*y)*intVFeO - y*intVFe2O3)/(R*t));
                double fo2 = exp((*logfo2)*log(10.0));

                temp =  (KD1*pow(fo2, (double) 0.25) + 2.0*y*K2*pow(KD1, 2.0*y)*pow(fo2, y/2.0))
                /(1.0 + (1.0-2.0*y)*K2*pow(KD1, 2.0*y)*pow(fo2, y/2.0));

                o[indexFe2O3]  = o[indexFeO]*(1.0-1.0/(1.0+temp))/2.0;
                o[indexFeO]   -= 2.0*o[indexFe2O3];
            }
#endif /* USE_KRESS_CARMICHAEL_FO2 */

        } else if (inpMask == FIRST && outMask == SEVENTH) {
            /*----------------------------------------------------------------------
                Calculates from the given t and p and a vector of moles of oxides
                (as defined in LIQ_STRUCT_DATA.H for the structure bulkSystem) the
                appropriate log10fo2 for the given t and p.
            ------------------------------------------------------------------------*/

            if (o[indexFeO] == 0.0 || o[indexFe2O3] == 0.0) { *logfo2 = 0.0; return; }
            for (i=0; i<nc; i++) sum += o[i];
            sum += o[indexFe2O3];
            if (sum == 0.0) { *logfo2 = 0.0; return; }

#ifdef USE_KRESS_CARMICHAEL_FO2
            if (p< 50000.0) temp = b/t + c + e*(1.0 - t0/t - log(t/t0)) + f*p/t + g*(t-t0)*p/t + h*SQUARE(p)/t;
            else            temp = b/t + c + e*(1.0 - t0/t - log(t/t0)) + f*50000.0/t + g*(t-t0)*50000.0/t + h*SQUARE(50000.0)/t
                             - a*log(10.0)*(608.966*p/10000.0-608.966*5.0)/t;
            for (i=0; i<nc; i++) temp += bulkSystem[i].coeff*o[i]/sum;
            temp += 2.0*bulkSystem[indexFeO].coeff*o[indexFe2O3]/sum
                                - bulkSystem[indexFe2O3].coeff*o[indexFe2O3]/sum;
            *logfo2 = (log(o[indexFe2O3]/o[indexFeO]) - temp)/(a*log(10.0));
#else
            {
                double y = 0.3;
                double intVFe2O3  = integralV_GKsp(iOxFe2O3,  t, p);
                double intVFeO    = integralV_GKsp(iOxFeO,    t, p);
                double intVFeO1_3 = integralV_GKsp(iOxFeO1_3, t, p);
                double deltaG = -106200.0 - t*(-55.1) + 31.86*(t - 1673.15 - t*log(t/1673.15)) + intVFe2O3/2.0 - intVFeO;
                double KD1 = exp(-deltaG/(R*t)
                            -(39860.0*o[iOxAl2O3] - 62520.0*o[iOxCaO] - 102000.0*o[iOxNa2O] - 119000.0*o[iOxK2O])/(R*t*sum));
                double K2 = 0.4*exp(-(intVFeO1_3 - (1.0-2.0*y)*intVFeO - y*intVFe2O3)/(R*t));
                int converged = FALSE;
                int iter = 0;

                *logfo2 = -10.0;

                while (!converged && (iter < 200)) {
                double fo2    = exp((*logfo2)*log(10.0));
                double dfo2   = exp((*logfo2)*log(10.0))*log(10.0);
                double numer  = KD1*pow(fo2, (double) 0.25) + 2.0*y*K2*pow(KD1, 2.0*y)*pow(fo2, y/2.0);
                double dnumer = 0.25*KD1*dfo2/pow(fo2, (double) 0.75) + 2.0*y*K2*pow(KD1, 2.0*y)*(y/2.0)*dfo2/pow(fo2, 1.0-y/2.0);
                double denom  = 1.0 + (1.0-2.0*y)*K2*pow(KD1, 2.0*y)*pow(fo2, y/2.0);
                double ddenom = (1.0-2.0*y)*K2*pow(KD1, 2.0*y)*(y/2.0)*dfo2/pow(fo2, 1.0-y/2.0);

                double f  =  numer/denom - 2.0*o[indexFe2O3]/o[indexFeO];
                double df = dnumer/denom - numer*ddenom/(denom*denom);
                double corr = -f/df;

                if (fabs(corr) > sqrt(DBL_EPSILON)) *logfo2 += corr; else converged = TRUE;
                if (*logfo2 >  10.0) *logfo2 =  10.0;
                if (*logfo2 < -50.0) *logfo2 = -50.0;
                iter++;
            }

    if (!converged) printf("Convergence failure in conLiq_CO2_H2O(FIRST,SEVENTH)\n");
            }
#endif /* USE_KRESS_CARMICHAEL_FO2 */

        } else
            printf("Illegal call to conLiq_CO2_H2O with inpMask = %o and outMask = %o\n",
                inpMask, outMask);

    } else if (inpMask == SECOND) {
        double sum;

        if (outMask & ~(THIRD | FOURTH | FIFTH | SIXTH))
            printf("Illegal call to conLiq_CO2_H2O with inpMask = %o and outMask = %o\n",
                inpMask, outMask);

        for (i=0, sum=0.0; i<NA; i++) sum += m[i];

        if (outMask & THIRD) {
            /* Converts a vector of moles of end-member components (m) into a vector
         of independent compositional variables (r) required as input for the
         remaining public functions.
         The dependent variable is taken to be SiO2 (1st component), as this
         component will never have a mole fraction of zero.                   */

            for (i=0; i<NR; i++) r[i] = (sum != 0.0) ? m[i+1]/sum : 0.0;
        }

        if (outMask & FOURTH) {
            /* Converts a vector of moles of end-member components (m) into a vector
         of mole fractions of endmember components                            */

            for (i=0; i<NA; i++) x[i] = (sum != 0.0) ? m[i]/sum : 0.0;
        }

        if (outMask & FIFTH) {
            /* Calculates the matrix dr[i]/dm[j] using m[] as input                 */

            if (sum == 0.0) {
                for (i=0; i<NR; i++) for (j=0; j<NA; j++) dm[i][j] = 0.0;
            } else {
                for (i=0; i<NR; i++) {
                    for (j=0; j<NA; j++) {
                        dm[i][j] = (i+1 == j) ? (1.0-m[i+1]/sum)/sum : - m[i+1]/SQUARE(sum);
                    }
                }
            }
        }

        if (outMask & SIXTH) {
            /* Calculates the matrix d2r[i]/dm[j]dm[k] using m[] as input           */

            if (sum == 0.0) {
                for (i=0; i<NR; i++) {
                    for (j=0; j<NA; j++)  {
                        for (k=0; k<NA; k++) d2m[i][j][k] = 0.0;
                    }
                }
            } else {
                for (i=0; i<NR; i++) {
                    for (j=0; j<NA; j++) {
                        for (k=0; k<NA; k++) {
                            d2m[i][j][k]  = 2.0*m[i+1]/CUBE(sum);
                            d2m[i][j][k] -= (i+1 == j) ? 1.0/SQUARE(sum) : 0.0;
                            d2m[i][j][k] -= (i+1 == k) ? 1.0/SQUARE(sum) : 0.0;
                        }
                    }
                }
            }

        }

    } else if (inpMask == THIRD && outMask == FOURTH) {
        /* Converts a vector of independent compositional variables (r)
            into a vector of mole fractions of end-member components (x)            */

        for (i=0, x[0] = 1.0; i<NR; i++) { x[0] -= r[i]; x[i+1] = r[i]; }

    } else if (inpMask == THIRD && outMask == (FOURTH | EIGHTH)) {
        /* Converts a vector of independent compositional variables (r)
            into a vector of mole fractions of endmember species (x)                */

        MTHREAD_MUTEX_LOCK(&global_data_mutex);
        order(0, t, p, r,
                NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
        for (i=0; i<NE; i++) x[i] = xSpecies[i];
        MTHREAD_MUTEX_UNLOCK(&global_data_mutex);

    } else if (inpMask == THIRD && outMask == (FOURTH | NINTH)) {
        double s[NT];
        /* Converts a vector of independent compositional variables (r)
            into a vector of fractions of system in CN state * (x)                */

        MTHREAD_MUTEX_LOCK(&global_data_mutex);
        order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
        for (i=0; i<NY; i++) x[i] = s[NS+i];
        MTHREAD_MUTEX_UNLOCK(&global_data_mutex);

    } else {
        printf("Illegal call to conLiq_CO2_H2O with inpMask = %o and outMask = %o\n",
            inpMask, outMask);
    }

}

int
testLiq_CO2_H2O(int mask, double t, double p,
    int na,          /* Expected number of endmember components                 */
    int nr,          /* Expected number of independent compositional variables  */
    char **names,    /* array of strings of names of endmember oxides           */
    char **formulas, /* array of strings of formulas of endmember components    */
    double *r,       /* array of indepependent compos variables, check bounds   */
    double *m)       /* array of moles of endmember components, check bounds    */
{
    const char *phase = "liquid.c";
    int result = TRUE, i;

    MTHREAD_ONCE(&initThreadBlock, threadInit);

    if (mask & FIRST) {
        result = result && (na == NA);
        if (!result) printf("<<%s>> Wrong number of components!\n", phase);
    }
    if (mask & SECOND) {
        result = result && (nr == NR);
        if (!result) printf("<<%s>> Wrong number of indep variables!\n", phase);
    }
    if (mask & THIRD) {
        for (i=0; i<NA; i++) {
            result = result && (strcmp(names[i],bulkSystem[i].label) == 0);
            if (!result)
                printf("<<%s>> Oxide[%d] should be %s not %s.\n",
                    phase, i, bulkSystem[i].label, names[i]);
        }
    }
    if (mask & FOURTH) {
        for (i=0; i<NA; i++) {
            result = result && (strcmp(formulas[i],liquid[i].label) == 0);
            if (!result)
                printf("<<%s>> Component[%d] should have formula %s not %s.\n",
                    phase, i, liquid[i].label, formulas[i]);
        }
    }
    /* Check bounds on the independent compositional variables */
    if (mask & FIFTH) {
        double s[NT];
        MTHREAD_MUTEX_LOCK(&global_data_mutex);
        initialGuessOrdering(r, s);
        result = rANDsTOx (r, s);
        MTHREAD_MUTEX_UNLOCK(&global_data_mutex);
    }
    /* Check bounds on moles of endmember components */
    if (mask & SIXTH) {
        double rTemp[NR], s[NT], sum;
        for (i=0, sum=0.0; i<NA; i++) sum += m[i];
        for (i=0; i<NR; i++) rTemp[i] = (sum != 0.0) ? m[i+1]/sum : 0.0;
        if (sum > 0.0) {
            MTHREAD_MUTEX_LOCK(&global_data_mutex);
            initialGuessOrdering(rTemp, s);
            result = rANDsTOx (rTemp, s);
            MTHREAD_MUTEX_UNLOCK(&global_data_mutex);
        } else {
            result = FALSE;
        }
    }
    /* Check if ordering state calculation has converged */
    if (mask & SEVENTH) { /* This call is NOT thread safe, but it should never be called from a threaded app */
        result = convergedInOrder;
    }
    /* Check if EOS calculation is valid */
    if (mask & EIGHTH) { /* This call is NOT thread safe, but it should never be called from a threaded app */
        result = TRUE;
    }
    /* Check if number of CN states is corrcet */
    if (mask & NINTH) {
        result = TRUE;
    }

    return result;
}

void
dispLiq_CO2_H2O(int mask, double t, double p, double *x,
    char **formula            /* Mineral formula for interface display MASK: 1 */
    )
{
    double *r = x;

    MTHREAD_ONCE(&initThreadBlock, threadInit);

    if (mask & FIRST) {    /* assume maximum string length is 5 */
        char *string = (char *) malloc((unsigned) (7+NA*12+1)*sizeof(char));;
        double m[NA], oxVal[NA], oxSum;
        int i, j, n;

        (void) snprintf(string, 8, "wt%% ox:");

        for (i=0, m[0] = 1.0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }
        for (i=0, oxSum=0.0; i<NA; i++) {
            for (j=0, oxVal[i]=0.0; j<NA; j++) oxVal[i] += m[j]*(liquid[j].liqToOx)[i];
            oxVal[i] *= bulkSystem[i].mw;
            oxSum    += oxVal[i];
        }

        if (oxSum != 0.0) for (i=0, n=7; i<NA; i++)
            if (oxVal[i] != 0.0) {
                double w = 100.0*oxVal[i]/oxSum;
                int nn = snprintf(&string[n], 13, " %s %.2f", bulkSystem[i].label, w);
    n += (nn < 13) ? nn : 12;
            }

        *formula = string;
    }
}

static int returnMixingProperties = TRUE;
void setModeToMixingLiq_CO2_H2O(int flag) { returnMixingProperties = flag; }

void
actLiq_CO2_H2O(int mask, double t, double p, double *x,
    double *a,   /* (pointer to a[]) activities              BINARY MASK: 000001 */
    double *mu,  /* (pointer to mu[]) chemical potentials    BINARY MASK: 000010 */
    double **dx, /* (pointer to dx[][]) d(a[])/d(x[])        BINARY MASK: 000100 */
    double **dw  /* (pointer to dw[][]) d(mu[])/d(w[])       BINARY MASK: 001000 */
    )
{
    double *r = x;
    double s[NT], g, dgdr[NR];
    double fr[NA][NR];
    int i, j;

    MTHREAD_ONCE(&initThreadBlock, threadInit);

    MTHREAD_MUTEX_LOCK(&global_data_mutex);
    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    for(i=0; i<NA; i++) for (j=0; j<NR; j++) fr[i][j] = rsEndmembers[i][j] - r[j];

    g = fillG(r, s, t, p);
    fillDGDR (r, s, t, p, dgdr);

    if (!mask && a != NULL) {
        for (i=0; i<NT; i++) a[i] = s[i];
    }

    if (mask & FIRST) {
        for(i=0; i<NA; i++) {
       for (a[i]=g, j=0; j<NR; j++) a[i] += fr[i][j]*dgdr[j];
       if (returnMixingProperties) a[i] -= G(i);
       a[i] = exp(a[i]/(R*t));
        }
    }

    if (mask & SECOND) {
        for(i=0; i<NA; i++) {
            for (mu[i]=g, j=0; j<NR; j++) mu[i] += fr[i][j]*dgdr[j];
            if (returnMixingProperties) mu[i] -= G(i);
        }

    }

    if (mask & THIRD) {
        double d2gdr2[NR][NR];
        double d2gdrds[NR][NT];
        double d2gds2[NT][NT];

        double dsdr[NT][NR], dfrdr[NA][NR], gs[NA][NT], aref[NA], sum;
        int k, l;

        fillD2GDR2  (r, s, t, p,  d2gdr2);
        fillD2GDRDS (r, s, t, p,  d2gdrds);
        fillD2GDS2  (r, s, t, p,  d2gds2);

        /* fill Darken structures */
        for(i=0; i<NA; i++) {
            for (j=0;  j<NR; j++) dfrdr[i][j] = -1.0;
            for (j=0;  j<NS; j++) gs[i][j] = rsEndmembers[i][NR+j] - s[j];
            for (j=NS; j<NT; j++) gs[i][j] = 0.0;
        }

        order(SECOND, t, p, r,
                    NULL, dsdr, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        /* get reference activities */
        if (!(mask & FIRST)) {
            for(i=0; i<NA; i++) {
                for (aref[i]=g, j=0; j<NR; j++) aref[i] += fr[i][j]*dgdr[j];
                if (returnMixingProperties) aref[i] -= G(i);
                aref[i] = exp(aref[i]/(R*t));
            }
        } else for (i=0; i<NA; i++) aref[i] = a[i];

        /* Compute derivatives of the chemical potentials */
        for (i=0; i<NA; i++) {
            for (k=0; k<NR; k++) {
                sum = (1.0+dfrdr[i][k])*dgdr[k];
                for (j=0; j<NR; j++) {
                    sum += fr[i][j]*d2gdr2[j][k];
                    for (l=0; l<NT; l++) sum += fr[i][j]*d2gdrds[j][l]*dsdr[l][k];
                }
                for (j=0; j<NT; j++) {
                    sum += gs[i][j]*d2gdrds[k][j];
                    for (l=0; l<NT; l++) sum += gs[i][j]*d2gds2[j][l]*dsdr[l][k];
                }
                dx[i][k] = sum; /* This is d mu/d r */
            }
        }

        /* convert result to d a/d r */
        for (i=0; i<NA; i++) for (j=0; j<NR; j++) dx[i][j] *= aref[i]/(R*t);

    }

    if (mask & FOURTH) {
        double dgdw[3*NP];
        double d2gdrds[NR][NT];
        double d2gdrdw[NR][3*NP];
        double dsdw[NT][3*NP], sum;
        int k, l;

        fillDGDW    (r, s, t, p,  dgdw);
        fillD2GDRDS (r, s, t, p,  d2gdrds);
        fillD2GDRDW (r, s, t, p,  d2gdrdw);

        order(ELEVENTH, t, p, r,
                    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, dsdw, NULL);

        for (i=0; i<NA; i++) {
            for (j=0; j<NP; j++) {
                if (meltsAndCO2_H2OModelParameters[j].activeH) {
                    for (k=0, dw[i][j]=dgdw[j]; k<NR; k++) {
                        for (l=0, sum=d2gdrdw[k][j]; l<NT; l++) sum += d2gdrds[k][l]*dsdw[l][j];
                        dw[i][j] += fr[i][k]*sum;
                    }
                } else dw[i][j] = 0.0;

                if (meltsAndCO2_H2OModelParameters[j].activeS) {
                    for (k=0, dw[i][NP+j]=dgdw[NP+j]; k<NR; k++) {
                        for (l=0, sum=d2gdrdw[k][NP+j]; l<NT; l++) sum += d2gdrds[k][l]*dsdw[l][NP+j];
                        dw[i][NP+j] += fr[i][k]*sum;
                    }
                } else dw[i][NP+j] = 0.0;

                if (meltsAndCO2_H2OModelParameters[j].activeV) {
                    for (k=0, dw[i][2*NP+j]=dgdw[2*NP+j]; k<NR; k++) {
                        for (l=0, sum=d2gdrdw[k][2*NP+j]; l<NT; l++) sum += d2gdrds[k][l]*dsdw[l][2*NP+j];
                        dw[i][2*NP+j] += fr[i][k]*sum;
                    }
                } else dw[i][2*NP+j] = 0.0;
            }

            if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
                if (meltsAndCO2_H2OModelParameters[NW+i].activeH) dw[i][     NW+i] -= 1.0;
                if (meltsAndCO2_H2OModelParameters[NW+i].activeS) dw[i][  NP+NW+i] -= -t;
                if (meltsAndCO2_H2OModelParameters[NW+i].activeV) dw[i][2*NP+NW+i] -= (p-1.0);
            }
        }
    }

    /************************************************************************
   * This return is not-public and is used for testing of the derivatives *
   * of G with respect to w[] internal to test_liquid.                    *
   ************************************************************************/
    if (mask & FIFTH) {
        fillDGDW (r, s, t, p, mu);

        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            mu[     NW+0] -= 1.0;
            mu[  NP+NW+0] -= -t;
            mu[2*NP+NW+0] -= (p-1.0);
            for (i=0; i<NR; i++) {
                mu[     NW  +0] +=         r[i];
                mu[  NP+NW  +0] +=      -t*r[i];
                mu[2*NP+NW  +0] += (p-1.0)*r[i];
                mu[     NW+i+1] -=         r[i];
                mu[  NP+NW+i+1] -=      -t*r[i];
                mu[2*NP+NW+i+1] -= (p-1.0)*r[i];
            }
        }

        for (i=0; i<(3*NP); i++) if (fabs(mu[i]) < 10.0*DBL_EPSILON) mu[i] = 0.0;
    }

    /*****************************************************************************
   * This return is used by the preclb_support.c functions to return the       *
   * configurational part of the activity.  It is now obsolete and is included *
   * for backwards compatibility.                                              *
   * ***************************************************************************/
    if (mask & SIXTH) {
        for (i=0; i<NA; i++) {
            a[i] = x[i];
            a[i] = (i != NA-1)   ? (1.0 - x[NA-1])*a[i] : x[NA-1]*a[NA-1];
        }
    }

    /*****************************************************************************
   * This return is not-public and is used within preclb_slave for             *
   * computation of the parameter derivatives of dgdr which are returned in dw.*
   *****************************************************************************/
    if (mask & SEVENTH) {
        double d2gdrds[NR][NT];
        double d2gdrdw[NR][3*NP];
        double dsdw[NT][3*NP];
        int k;

        fillD2GDRDS (r, s, t, p,  d2gdrds);
        fillD2GDRDW (r, s, t, p,  d2gdrdw);

        order(ELEVENTH, t, p, r,
                    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, dsdw, NULL);

        for (i=0; i<NR; i++) {
            for (j=0; j<NP; j++) {
                if (meltsAndCO2_H2OModelParameters[j].activeH) {
                    for (k=0, dw[i][j]=d2gdrdw[i][j]; k<NT; k++) dw[i][j] += d2gdrds[i][k]*dsdw[k][j];
                } else dw[i][j] = 0.0;

                if (meltsAndCO2_H2OModelParameters[j].activeS) {
                    for (k=0, dw[i][NP+j]=d2gdrdw[i][NP+j]; k<NT; k++) dw[i][NP+j] += d2gdrds[i][k]*dsdw[k][NP+j];
                } else dw[i][NP+j] = 0.0;

                if (meltsAndCO2_H2OModelParameters[j].activeV) {
                    for (k=0, dw[i][2*NP+j]=d2gdrdw[i][2*NP+j]; k<NT; k++) dw[i][2*NP+j] += d2gdrds[i][k]*dsdw[k][2*NP+j];
                } else dw[i][2*NP+j] = 0.0;
            }

            if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties dr[i] += (G(0)-G(i+1)); */
                if (meltsAndCO2_H2OModelParameters[NW  +0].activeH) dw[i][     NW  +0] += 1.0;
                if (meltsAndCO2_H2OModelParameters[NW  +0].activeS) dw[i][  NP+NW  +0] += -t;
                if (meltsAndCO2_H2OModelParameters[NW  +0].activeV) dw[i][2*NP+NW  +0] += (p-1.0);
                if (meltsAndCO2_H2OModelParameters[NW+i+1].activeH) dw[i][     NW+i+1] -= 1.0;
                if (meltsAndCO2_H2OModelParameters[NW+i+1].activeS) dw[i][  NP+NW+i+1] -= -t;
                if (meltsAndCO2_H2OModelParameters[NW+i+1].activeV) dw[i][2*NP+NW+i+1] -= (p-1.0);
            }
        }
    }

    /*******************************************************************************
    * This return is non-public and is designed for xMELTS calibration purposes.   *
    * The negative of the partial molar entropy (i.e. d mu dT) is returned in mu[] *
    ********************************************************************************/
    if (mask & EIGHTH) {
        double d2gdrds[NR][NT];
        double d2gdrdt[NR];
        double d2gds2[NT][NT];
        double d2gdsdt[NT];
        double dsdt[NT], gs[NA][NT], sum, dgdt;
        int k;

        dgdt = -fillS(r, s, t, p);
        fillD2GDRDS (r, s, t, p,  d2gdrds);
        fillD2GDRDT (r, s, t, p,  d2gdrdt);
        fillD2GDS2  (r, s, t, p,  d2gds2);
        fillD2GDSDT (r, s, t, p,  d2gdsdt);

        /* fill Darken structures */
        for(i=0; i<NA; i++) {
            for (j=0;  j<NS; j++) gs[i][j] = rsEndmembers[i][NR+j] - s[j];
            for (j=NS; j<NT; j++) gs[i][j] = 0.0;
        }

        order(THIRD, t, p, r,
                    NULL, NULL, dsdt, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        /* Compute derivatives of the chemical potentials */
        for (i=0; i<NA; i++) {
            sum = dgdt;
            for (j=0; j<NR; j++) {
                sum += fr[i][j]*d2gdrdt[j];
                for (k=0; k<NT; k++) sum += fr[i][j]*d2gdrds[j][k]*dsdt[k];
            }
            for (j=0; j<NT; j++) {
                sum += gs[i][j]*d2gdsdt[j];
                for (k=0; k<NT; k++) sum += gs[i][j]*d2gds2[j][k]*dsdt[k];
            }
            mu[i] = sum; /* This is d mu/d t */
            if (returnMixingProperties) mu[i] -= -S(i);
        }

    }

    /*******************************************************************************
    * This return is non-public and is designed for xMELTS calibration purposes.   *
    * The model parameter derivative of the negative of the partial molar entropy  *
    * (i.e. d mu / dT dW) is returned in dw[][]                                    *
    ********************************************************************************/
    if (mask & NINTH) {
        double d2gdsdw[NT][3*NP];
        double d2gdtdw[3*NP];
        double d2gdrds[NR][NT];
        double d2gds2[NT][NT];
        double d2gdsdt[NT];
        double d3gdrds2[NR][NT][NT];
        double d3gdrdsdt[NR][NT];
        double d3gdrdsdw[NR][NT][3*NP];
        double d3gdrdtdw[NR][3*NP];
        double d3gds3[NT][NT][NT];
        double d3gds2dt[NT][NT];
        double d3gds2dw[NT][NT][3*NP];
        double d3gdsdtdw[NT][3*NP];
        double dsdt[NT], dsdw[NT][3*NP], d2sdtdw[NT][3*NP], gs[NA][NT];
        int k, l, m;

        fillD2GDSDW   (r, s, t, p,  d2gdsdw);
        fillD2GDTDW   (r, s, t, p,  d2gdtdw);
        fillD2GDRDS   (r, s, t, p,  d2gdrds);
        fillD2GDS2    (r, s, t, p,  d2gds2);
        fillD2GDSDT   (r, s, t, p,  d2gdsdt);
        fillD3GDRDS2  (r, s, t, p,  d3gdrds2);
        fillD3GDRDSDT (r, s, t, p,  d3gdrdsdt);
        fillD3GDRDSDW (r, s, t, p,  d3gdrdsdw);
        fillD3GDRDTDW (r, s, t, p,  d3gdrdtdw);
        fillD3GDS3    (r, s, t, p,  d3gds3);
        fillD3GDS2DT  (r, s, t, p,  d3gds2dt);
        fillD3GDS2DW  (r, s, t, p,  d3gds2dw);
        fillD3GDSDTDW (r, s, t, p,  d3gdsdtdw);

        /* fill Darken structures */
        for(i=0; i<NA; i++) {
            for (j=0;  j<NS; j++) gs[i][j] = rsEndmembers[i][NR+j] - s[j];
            for (j=NS; j<NT; j++) gs[i][j] = 0.0;
        }

        order(THIRD | ELEVENTH | TWELFTH, t, p, r,
                    NULL, NULL, dsdt, NULL, NULL, NULL, NULL, NULL, NULL, NULL, dsdw, d2sdtdw);

        for (i=0; i<NA; i++) {
            for (j=0; j<NP; j++) {
                if (meltsAndCO2_H2OModelParameters[j].activeH) {
                    dw[i][j] = d2gdtdw[j];
                    for (l=0; l<NT; l++) dw[i][j] += d2gdsdt[l]*dsdw[l][j];
                    for (k=0; k<NR; k++) {
                        dw[i][j] += fr[i][k]*d3gdrdtdw[k][j];
                        for (l=0; l<NT; l++) dw[i][j] += fr[i][k]*d3gdrdsdt[k][l]*dsdw[l][j];
                        for (m=0; m<NT; m++) {
                            dw[i][j] += fr[i][k]*(d3gdrdsdw[k][m][j]*dsdt[m]+d2gdrds[k][m]*d2sdtdw[m][j]);
                            for (l=0; l<NT; l++) if (dsdw[l][j] != 0.0) dw[i][j] += fr[i][k]*d3gdrds2[k][m][l]*dsdw[l][j]*dsdt[m];
                        }
                    }
                    for (k=0; k<NT; k++) {
                        dw[i][j] += gs[i][k]*d3gdsdtdw[k][j];
                        for (l=0; l<NT; l++) dw[i][j] += gs[i][k]*d3gds2dt[k][l]*dsdw[l][j];
                        for (m=0; m<NT; m++) {
                            dw[i][j] += gs[i][k]*(d3gds2dw[k][m][j]*dsdt[m]+d2gds2[k][m]*d2sdtdw[m][j]);
                            for (l=0; l<NT; l++) if (dsdw[l][j] != 0.0) dw[i][j] += gs[i][k]*d3gds3[k][m][l]*dsdw[l][j]*dsdt[m];
                        }
                    }
                } else dw[i][j] = 0.0;

                if (meltsAndCO2_H2OModelParameters[j].activeS) {
                    dw[i][NP+j] = d2gdtdw[NP+j];
                    for (l=0; l<NT; l++) dw[i][NP+j] += d2gdsdt[l]*dsdw[l][NP+j];
                    for (k=0; k<NR; k++) {
                        dw[i][NP+j] += fr[i][k]*d3gdrdtdw[k][NP+j];
                        for (l=0; l<NT; l++) dw[i][NP+j] += fr[i][k]*d3gdrdsdt[k][l]*dsdw[l][NP+j];
                        for (m=0; m<NT; m++) {
                            dw[i][NP+j] += fr[i][k]*(d3gdrdsdw[k][m][NP+j]*dsdt[m]+d2gdrds[k][m]*d2sdtdw[m][NP+j]);
                            for (l=0; l<NT; l++) if (dsdw[l][NP+j] != 0.0) dw[i][NP+j] += fr[i][k]*d3gdrds2[k][m][l]*dsdw[l][NP+j]*dsdt[m];
                        }
                    }
                    for (k=0; k<NT; k++) {
                        dw[i][NP+j] += gs[i][k]*d3gdsdtdw[k][NP+j];
                        for (l=0; l<NT; l++) dw[i][NP+j] += gs[i][k]*d3gds2dt[k][l]*dsdw[l][NP+j];
                        for (m=0; m<NT; m++) {
                            dw[i][NP+j] += gs[i][k]*(d3gds2dw[k][m][NP+j]*dsdt[m]+d2gds2[k][m]*d2sdtdw[m][NP+j]);
                            for (l=0; l<NT; l++) if (dsdw[l][NP+j] != 0.0) dw[i][NP+j] += gs[i][k]*d3gds3[k][m][l]*dsdw[l][NP+j]*dsdt[m];
                        }
                    }
                } else dw[i][NP+j] = 0.0;

                if (meltsAndCO2_H2OModelParameters[j].activeV) {
                    dw[i][2*NP+j] = d2gdtdw[2*NP+j];
                    for (l=0; l<NT; l++) dw[i][2*NP+j] += d2gdsdt[l]*dsdw[l][2*NP+j];
                    for (k=0; k<NR; k++) {
                        dw[i][2*NP+j] += fr[i][k]*d3gdrdtdw[k][2*NP+j];
                        for (l=0; l<NT; l++) dw[i][2*NP+j] += fr[i][k]*d3gdrdsdt[k][l]*dsdw[l][2*NP+j];
                        for (m=0; m<NT; m++) {
                            dw[i][2*NP+j] += fr[i][k]*(d3gdrdsdw[k][m][2*NP+j]*dsdt[m]+d2gdrds[k][m]*d2sdtdw[m][2*NP+j]);
                            for (l=0; l<NT; l++) if (dsdw[l][2*NP+j] != 0.0) dw[i][2*NP+j] += fr[i][k]*d3gdrds2[k][m][l]*dsdw[l][2*NP+j]*dsdt[m];
                        }
                    }
                    for (k=0; k<NT; k++) {
                        dw[i][2*NP+j] += gs[i][k]*d3gdsdtdw[k][2*NP+j];
                        for (l=0; l<NT; l++) dw[i][2*NP+j] += gs[i][k]*d3gds2dt[k][l]*dsdw[l][2*NP+j];
                        for (m=0; m<NT; m++) {
                            dw[i][2*NP+j] += gs[i][k]*(d3gds2dw[k][m][2*NP+j]*dsdt[m]+d2gds2[k][m]*d2sdtdw[m][2*NP+j]);
                            for (l=0; l<NT; l++) if (dsdw[l][2*NP+j] != 0.0) dw[i][2*NP+j] += gs[i][k]*d3gds3[k][m][l]*dsdw[l][2*NP+j]*dsdt[m];
                        }
                    }
                } else dw[i][2*NP+j] = 0.0;
            }

            if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
                if (meltsAndCO2_H2OModelParameters[NW+i].activeS) dw[i][  NP+NW+i] -= -1.0;
            }

        }
    }

    MTHREAD_MUTEX_UNLOCK(&global_data_mutex);
}

void
gmixLiq_CO2_H2O(int mask, double t, double p, double *x,
    double *gmix, /* Gibbs energy of mixing             BINARY MASK: 0001 */
    double *dx,   /* (pointer to dx[]) d(g)/d(x[])      BINARY MASK: 0010 */
    double **dx2  /* (pointer to dx2[][]) d2(g)/d(x[])2 BINARY MASK: 0100 */
    )
{
    double *r = x;
    double s[NT];
    int i;

    MTHREAD_ONCE(&initThreadBlock, threadInit);

    MTHREAD_MUTEX_LOCK(&global_data_mutex);
    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    if (mask & FIRST) {
        *gmix  = fillG (r, s, t, p);
        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            *gmix -= G(0);
            for (i=0; i<NR; i++) *gmix += r[i]*(G(0)-G(i+1));
        }
    }

    if(mask & SECOND) {
        fillDGDR (r, s, t, p, dx);
        if (returnMixingProperties) for (i=0; i<NR; i++) dx[i] += (G(0)-G(i+1));
    }

    if(mask & THIRD) {
        double d2gdr2[NR][NR];
        double d2gdrds[NR][NT];
        double d2gds2[NT][NT];
        double dsdr[NT][NR];
        int j, k, l;

        fillD2GDR2  (r, s, t, p,  d2gdr2);
        fillD2GDRDS (r, s, t, p,  d2gdrds);
        fillD2GDS2  (r, s, t, p,  d2gds2);

        order(SECOND, t, p, r,
                    NULL, dsdr, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            for (j=0; j<NR; j++) {
                dx2[i][j] = d2gdr2[i][j];
                for (k=0; k<NT; k++) {
                    dx2[i][j] += d2gdrds[i][k]*dsdr[k][j] + d2gdrds[j][k]*dsdr[k][i];
                    for (l=0; l<NT; l++) dx2[i][j] += d2gds2[k][l]*dsdr[k][i]*dsdr[l][j];
                }
            }
        }
    }

    if(mask & FOURTH) {
        double d3gdr3[NR][NR][NR];
        double d3gdr2ds[NR][NR][NT];
        double d3gdrds2[NR][NT][NT];
        double d3gds3[NT][NT][NT];
        double dx3[NR][NR][NR]; /* This should be passed to the function */
        double dsdr[NT][NR], d2sdr2[NT][NR][NR];
        int i, j, k, l, m, n;

        fillD3GDR3   (r, s, t, p,  d3gdr3);
        fillD3GDR2DS (r, s, t, p,  d3gdr2ds);
        fillD3GDRDS2 (r, s, t, p,  d3gdrds2);
        fillD3GDS3   (r, s, t, p,  d3gds3);

        order(SECOND, t, p, r,
                    NULL, dsdr, NULL, NULL, d2sdr2, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            for (j=0; j<NR; j++) {
                for (k=0; k<NR; k++) {
                    dx3[i][j][k] = d3gdr3[i][j][k];
                    for (l=0; l<NT; l++) {
                        dx3[i][j][k] += d3gdr2ds[i][j][l]*dsdr[l][k] +
                            d3gdr2ds[j][k][l]*dsdr[l][i] + d3gdr2ds[k][i][l]*dsdr[l][j];
                        for (m=0; m<NT; m++) {
                            dx3[i][j][k] +=
                                d3gdrds2[i][l][m]*dsdr[l][j]*dsdr[m][k] +
                                d3gdrds2[j][l][m]*dsdr[l][k]*dsdr[m][i] +
                                d3gdrds2[k][l][m]*dsdr[l][i]*dsdr[m][j];
                            for (n=0; n<NT; n++)
                                dx3[i][j][k] +=
                                    d3gds3[l][m][n]*dsdr[l][i]*dsdr[m][j]*dsdr[n][k];
                        }
                    }
                }
            }
        }
    }

    MTHREAD_MUTEX_UNLOCK(&global_data_mutex);
}

void
hmixLiq_CO2_H2O(int mask, double t, double p, double *x,
    double *hmix, /* Enthalpy of mixing BINARY MASK:                 01 */
    double *dw    /* (pointer to dw[]) d(H)/d(x[])      BINARY MASK: 10 */
    )
{
    double *r = x;
    double s[NT], dsdt[NT];
    int i;

    MTHREAD_ONCE(&initThreadBlock, threadInit);

    MTHREAD_MUTEX_LOCK(&global_data_mutex);
    order(FIRST | THIRD, t, p, r,
                s, NULL, dsdt, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    if (mask & FIRST) {
        *hmix = fillG (r, s, t, p) + t*fillS (r, s, t, p);                            /* was: *hmix = fillH (r, s, t, p);  */

        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            *hmix -= (G(0) + t*S(0));                                                   /* was: *hmix -= H(0);               */
            for (i=0; i<NR; i++) *hmix += r[i]*((G(0) + t*S(0))-(G(i+1) + t*S(i+1)));   /* was: *hmix += r[i]*(H(0)-H(i+1)); */
        }
    }

    if (mask & SECOND) {
        double dgdw[3*NP];
        double d2gdsdw[NT][3*NP];
        double d2gdtdw[3*NP];
        double d2gds2[NT][NT];
        double d2gdsdt[NT];
        double dsdw[NT][3*NP], dsdt[NT];
        int i, k, l;

        fillDGDW    (r, s, t, p, dgdw);
        fillD2GDSDW (r, s, t, p, d2gdsdw);
        fillD2GDTDW (r, s, t, p, d2gdtdw);
        fillD2GDS2  (r, s, t, p, d2gds2);
        fillD2GDSDT (r, s, t, p, d2gdsdt);

        order(THIRD | ELEVENTH, t, p, r,
                    NULL, NULL, dsdt, NULL, NULL, NULL, NULL, NULL, NULL, NULL, dsdw, NULL);

        for (i=0; i<NP; i++) {
            if (meltsAndCO2_H2OModelParameters[i].activeH) {
                dw[i] = d2gdtdw[i];
                for (k=0; k<NT; k++) {
                    dw[i] += d2gdsdw[k][i]*dsdt[k] + d2gdsdt[k]*dsdw[k][i];
                    for (l=0; l<NT; l++) dw[i] += d2gds2[k][l]*dsdt[k]*dsdw[l][i] ;
                }
                dw[i] = dgdw[i] - t*dw[i];
            } else dw[i] = 0.0;

            if (meltsAndCO2_H2OModelParameters[i].activeS) {
                dw[NP+i] = d2gdtdw[NP+i];
                for (k=0; k<NT; k++) {
                    dw[NP+i] += d2gdsdw[k][NP+i]*dsdt[k] + d2gdsdt[k]*dsdw[k][NP+i];
                    for (l=0; l<NT; l++) dw[NP+i] += d2gds2[k][l]*dsdt[k]*dsdw[l][NP+i] ;
                }
                dw[NP+i] = dgdw[NP+i] - t*dw[NP+i];
            } else dw[NP+i] = 0.0;

            if (meltsAndCO2_H2OModelParameters[i].activeV) {
                dw[2*NP+i] = d2gdtdw[2*NP+i];
                for (k=0; k<NT; k++) {
                    dw[2*NP+i] += d2gdsdw[k][2*NP+i]*dsdt[k] + d2gdsdt[k]*dsdw[k][2*NP+i];
                    for (l=0; l<NT; l++) dw[2*NP+i] += d2gds2[k][l]*dsdt[k]*dsdw[l][2*NP+i] ;
                }
                dw[2*NP+i] = dgdw[2*NP+i] - t*dw[2*NP+i];
            } else dw[2*NP+i] = 0.0;
        }

        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            if (meltsAndCO2_H2OModelParameters[NW+0].activeH) dw[     NW+0] -= 1.0;
            if (meltsAndCO2_H2OModelParameters[NW+0].activeV) dw[2*NP+NW+0] -= (p-1.0);
            for (i=0; i<NR; i++) {
                if (meltsAndCO2_H2OModelParameters[NW  +0].activeH) dw[     NW  +0] +=         r[i];
                if (meltsAndCO2_H2OModelParameters[NW  +0].activeV) dw[2*NP+NW  +0] += (p-1.0)*r[i];
                if (meltsAndCO2_H2OModelParameters[NW+i+1].activeH) dw[     NW+i+1] -=         r[i];
                if (meltsAndCO2_H2OModelParameters[NW+i+1].activeV) dw[2*NP+NW+i+1] -= (p-1.0)*r[i];
            }
        }

        for (i=0; i<(3*NP); i++) if (fabs(dw[i]) < 10.0*DBL_EPSILON) dw[i] = 0.0;
    }

    MTHREAD_MUTEX_UNLOCK(&global_data_mutex);
}

void
smixLiq_CO2_H2O(int mask, double t, double p, double *x,
    double *smix, /* Entropy of mixing                  BINARY MASK: 0001 */
    double *dx,   /* (pointer to dx[]) d(s)/d(x[])      BINARY MASK: 0010 */
    double **dx2, /* (pointer to dx2[][]) d2(s)/d(x[])2 BINARY MASK: 0100 */
    double *dw    /* (pointer to dw[]) d(s)/d(x[])      BINARY MASK: 1000 */
    )
{
    double *r = x;
    double s[NT];

    MTHREAD_ONCE(&initThreadBlock, threadInit);

    MTHREAD_MUTEX_LOCK(&global_data_mutex);
    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    if (mask & FIRST) {
        int i;

        *smix = fillS (r, s, t, p);

        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            *smix -= S(0);
            for (i=0; i<NR; i++) *smix += r[i]*(S(0)-S(i+1));
        }
    }

    if(mask & SECOND) {
        double d2gdrds[NR][NT];
        double d2gdrdt[NR];
        double d2gds2[NT][NT];
        double d2gdsdt[NT];
        double dsdr[NT][NR], dsdt[NT];
        int i, k, l;

        fillD2GDRDS (r, s, t, p, d2gdrds);
        fillD2GDRDT (r, s, t, p, d2gdrdt);
        fillD2GDS2  (r, s, t, p, d2gds2);
        fillD2GDSDT (r, s, t, p, d2gdsdt);

        order(SECOND | THIRD, t, p, r,
                    NULL, dsdr, dsdt, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            dx[i] = d2gdrdt[i];
            for (k=0; k<NT; k++) {
                dx[i] += d2gdrds[i][k]*dsdt[k] + d2gdsdt[k]*dsdr[k][i];
                for (l=0; l<NT; l++) dx[i] += d2gds2[k][l]*dsdt[k]*dsdr[l][i] ;
            }
            dx[i] *= -1.0;
        }

        /* Convert Solution Properties -> Mixing Properties */
        if (returnMixingProperties) for (i=0; i<NR; i++) dx[i] += (S(0)-S(i+1));
    }

    if(mask & THIRD) {
        double d2gdrds[NR][NT];
        double d2gds2[NT][NT];
        double d2gdsdt[NT];
        double d3gdr2ds[NR][NR][NT];
        double d3gdr2dt[NR][NR];
        double d3gdrds2[NR][NT][NT];
        double d3gdrdsdt[NR][NT];
        double d3gds3[NT][NT][NT];
        double d3gds2dt[NT][NT];
        double dsdr[NT][NR], dsdt[NT], d2sdr2[NT][NR][NR], d2sdrdt[NT][NR];
        int i, j, k, l, m;

        fillD2GDRDS   (r, s, t, p, d2gdrds);
        fillD2GDS2    (r, s, t, p, d2gds2);
        fillD2GDSDT   (r, s, t, p, d2gdsdt);
        fillD3GDR2DS  (r, s, t, p, d3gdr2ds);
        fillD3GDR2DT  (r, s, t, p, d3gdr2dt);
        fillD3GDRDS2  (r, s, t, p, d3gdrds2);
        fillD3GDRDSDT (r, s, t, p, d3gdrdsdt);
        fillD3GDS3    (r, s, t, p, d3gds3);
        fillD3GDS2DT  (r, s, t, p, d3gds2dt);

        order(SECOND | THIRD | FIFTH | SIXTH, t, p, r,
                    NULL, dsdr, dsdt, NULL, d2sdr2, d2sdrdt, NULL, NULL, NULL, NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            for (j=0; j<NR; j++) {
                dx2[i][j] = d3gdr2dt[i][j];
                for (k=0; k<NT; k++) {
                    dx2[i][j] += d3gdr2ds[i][j][k]*dsdt[k]
                     + d3gdrdsdt[i][k]*dsdr[k][j]
                     + d3gdrdsdt[j][k]*dsdr[k][i]
                     + d2gdsdt[k]*d2sdr2[k][i][j]
                     + d2gdrds[i][k]*d2sdrdt[k][j]
                     + d2gdrds[j][k]*d2sdrdt[k][i];
                    for (l=0; l<NT; l++) {
                        dx2[i][j] += d3gdrds2[i][k][l]*dsdr[k][j]*dsdt[l]
                       + d3gdrds2[j][k][l]*dsdr[k][i]*dsdt[l]
                       + d2gds2[k][l]*d2sdr2[k][i][j]*dsdt[l]
                       + d3gds2dt[k][l]*dsdr[k][i]*dsdr[l][j]
                       + d2gds2[k][l]*dsdr[k][i]*d2sdrdt[l][j]
                       + d2gds2[k][l]*dsdr[k][j]*d2sdrdt[l][i];
                        for (m=0; m<NT; m++)
                            dx2[i][j] += d3gds3[k][l][m]*dsdr[k][i]*dsdr[l][j]*dsdt[m];
                    }
                }
                dx2[i][j] *= -1.0;
            }
        }
    }

    if(mask & FOURTH) {
        double d2gdsdw[NT][3*NP];
        double d2gdtdw[3*NP];
        double d2gds2[NT][NT];
        double d2gdsdt[NT];
        double dsdw[NT][3*NP], dsdt[NT];
        int i, k, l;

        fillD2GDSDW (r, s, t, p, d2gdsdw);
        fillD2GDTDW (r, s, t, p, d2gdtdw);
        fillD2GDS2  (r, s, t, p, d2gds2);
        fillD2GDSDT (r, s, t, p, d2gdsdt);

        order(THIRD | ELEVENTH, t, p, r,
                    NULL, NULL, dsdt, NULL, NULL, NULL, NULL, NULL, NULL, NULL, dsdw, NULL);

        for (i=0; i<NP; i++) {
            if (meltsAndCO2_H2OModelParameters[i].activeH) {
                dw[i] = d2gdtdw[i];
                for (k=0; k<NT; k++) {
                    dw[i] += d2gdsdw[k][i]*dsdt[k] + d2gdsdt[k]*dsdw[k][i];
                    for (l=0; l<NT; l++) dw[i] += d2gds2[k][l]*dsdt[k]*dsdw[l][i] ;
                }
                dw[i] *= -1.0;
            } else dw[i] = 0.0;

            if (meltsAndCO2_H2OModelParameters[i].activeS) {
                dw[NP+i] = d2gdtdw[NP+i];
                for (k=0; k<NT; k++) {
                    dw[NP+i] += d2gdsdw[k][NP+i]*dsdt[k] + d2gdsdt[k]*dsdw[k][NP+i];
                    for (l=0; l<NT; l++) dw[NP+i] += d2gds2[k][l]*dsdt[k]*dsdw[l][NP+i] ;
                }
                dw[NP+i] *= -1.0;
            } else dw[NP+i] = 0.0;

            if (meltsAndCO2_H2OModelParameters[i].activeV) {
                dw[2*NP+i] = d2gdtdw[2*NP+i];
                for (k=0; k<NT; k++) {
                    dw[2*NP+i] += d2gdsdw[k][2*NP+i]*dsdt[k] + d2gdsdt[k]*dsdw[k][2*NP+i];
                    for (l=0; l<NT; l++) dw[2*NP+i] += d2gds2[k][l]*dsdt[k]*dsdw[l][2*NP+i] ;
                }
                dw[2*NP+i] *= -1.0;
            } else dw[2*NP+i] = 0.0;
        }

        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            if (meltsAndCO2_H2OModelParameters[NW+0].activeS) dw[NP+NW+0] -= 1.0;
            for (i=0; i<NR; i++) {
                if (meltsAndCO2_H2OModelParameters[NW  +0].activeS) dw[NP+NW  +0] += r[i];
                if (meltsAndCO2_H2OModelParameters[NW+i+1].activeS) dw[NP+NW+i+1] -= r[i];
            }
        }

        for (i=0; i<(3*NP); i++) if (fabs(dw[i]) < 10.0*DBL_EPSILON) dw[i] = 0.0;
    }

    MTHREAD_MUTEX_UNLOCK(&global_data_mutex);
}

void
cpmixLiq_CO2_H2O(int mask, double t, double p, double *x,
    double *cpmix, /* Heat capacity of mixing               BINARY MASK: 001 */
    double *dt,    /* d(cp)/d(t)                            BINARY MASK: 010 */
    double *dx     /* d(cp)/d(x[])d(t)                      BINARY MASK: 100 */
    )
{
    double *r = x;
    double d2gdsdt[NT];
    double d2gds2[NT][NT];
    double s[NT], dsdt[NT];
    double d2gdt2;
    int i, j;

    MTHREAD_ONCE(&initThreadBlock, threadInit);

    MTHREAD_MUTEX_LOCK(&global_data_mutex);
    order(FIRST | THIRD, t, p, r,
                s,    NULL, dsdt, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    fillD2GDSDT (r, s, t, p, d2gdsdt);
    fillD2GDS2  (r, s, t, p, d2gds2);
    d2gdt2 = fillD2GDT2  (r, s, t, p);

    if (mask & FIRST) {
        *cpmix = d2gdt2;
        for (i=0; i<NT; i++) {
            *cpmix += 2.0*d2gdsdt[i]*dsdt[i];
            for (j=0; j<NT; j++) *cpmix += d2gds2[i][j]*dsdt[i]*dsdt[j];
        }
        *cpmix *= -t;

        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            *cpmix -= CP(0);
            for (i=0; i<NR; i++) *cpmix += r[i]*(CP(0)-CP(i+1));
        }
    }

    if(mask & SECOND) {
        double d3gds3[NT][NT][NT];
        double d3gds2dt[NT][NT];
        double d3gdsdt2[NT];
        double d3gdt3 = fillD3GDT3   (r, s, t, p);
        double d2sdt2[NT], temp;
        int k;

        fillD3GDS3   (r, s, t, p, d3gds3);
        fillD3GDS2DT (r, s, t, p, d3gds2dt);
        fillD3GDSDT2 (r, s, t, p, d3gdsdt2);

        order(EIGHTH, t, p, r,
                    NULL, NULL, NULL, NULL, NULL, NULL, NULL, d2sdt2, NULL, NULL, NULL, NULL);

        /* compute d2gdt2 */
        temp = d2gdt2;
        for (i=0; i<NT; i++) {
            temp += 2.0*d2gdsdt[i]*dsdt[i];
            for (j=0; j<NT; j++) temp += d2gds2[i][j]*dsdt[i]*dsdt[j];
        }

        *dt = d3gdt3;
        for (i=0; i<NT; i++) {
            *dt += 3.0*d3gdsdt2[i]*dsdt[i] + 3.0*d2gdsdt[i]*d2sdt2[i];
            for (j=0; j<NT; j++) {
                *dt += 3.0*d2gds2[i][j]*dsdt[i]*d2sdt2[j]
             + 3.0*d3gds2dt[i][j]*dsdt[i]*dsdt[j];
                for (k=0; k<NT; k++) *dt += d3gds3[i][j][k]*dsdt[i]*dsdt[j]*dsdt[k];
            }
        }
        *dt = -t*(*dt) - temp;

        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            *dt -= DCPDT(0);
            for (i=0; i<NR; i++) *dt += r[i]*(DCPDT(0)-DCPDT(i+1));
        }
    }

    if(mask & THIRD) {
        double d3gds3[NT][NT][NT];
        double d3gdrds2[NR][NT][NT];
        double d3gdrdsdt[NR][NT];
        double d3gds2dt[NT][NT];
        double d2gdrds[NR][NT];
        double d3gdrdt2[NR];
        double d3gdsdt2[NT];
        double dsdr[NT][NR], d2sdrdt[NT][NR], d2sdt2[NT];
        int k, l;

        fillD3GDS3    (r, s, t, p, d3gds3);
        fillD3GDRDS2  (r, s, t, p, d3gdrds2);
        fillD3GDRDSDT (r, s, t, p, d3gdrdsdt);
        fillD3GDS2DT  (r, s, t, p, d3gds2dt);
        fillD2GDRDS   (r, s, t, p, d2gdrds);
        fillD3GDRDT2  (r, s, t, p, d3gdrdt2);
        fillD3GDSDT2  (r, s, t, p, d3gdsdt2);

        order(SECOND | SIXTH | EIGHTH, t, p, r,
                    NULL, dsdr, NULL, NULL, NULL, d2sdrdt, NULL, d2sdt2, NULL, NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            for (j=0,dx[i]=d3gdrdt2[i]; j<NT; j++) {
                dx[i] += d3gdsdt2[j]*dsdr[j][i] + 2.0*d2gdsdt[j]*d2sdrdt[j][i] +
                 2.0*d3gdrdsdt[i][j]*dsdt[j] + d2gdrds[i][j]*d2sdt2[j];
                for (k=0; k<NT; k++) {
                    dx[i] += d3gdrds2[i][j][k]*dsdt[j]*dsdt[k] +
                   2.0*d2gds2[j][k]*dsdt[j]*d2sdrdt[k][i] +
                   2.0*d3gds2dt[j][k]*dsdr[j][i]*dsdt[k] +
                   d2gds2[j][k]*dsdr[j][i]*d2sdt2[k];
                    for (l=0; l<NT; l++)
                        dx[i] += d3gds3[j][k][l]*dsdr[j][i]*dsdt[k]*dsdt[l];
                }
            }
            dx[i] *= -t;
        }

        if (returnMixingProperties) for (i=0; i<NR; i++) dx[i] += (CP(0)-CP(i+1));
    }

    MTHREAD_MUTEX_UNLOCK(&global_data_mutex);
}

void
vmixLiq_CO2_H2O(int mask, double t, double p, double *x,
    double *vmix, /* Volume of mixing                BINARY MASK: 00000000001 */
    double *dx,   /* (pointer to dx[]) d(v)/d(x[])   BINARY MASK: 00000000010 */
    double **dx2, /* (point to dx2[][]) d(v)/d(x[])2 BINARY MASK: 00000000100 */
    double *dt,   /* d(v)/d(t)                       BINARY MASK: 00000001000 */
    double *dp,   /* d(v)/d(p)                       BINARY MASK: 00000010000 */
    double *dt2,  /* d2(v)/d(t)2                     BINARY MASK: 00000100000 */
    double *dtdp, /* d2(v)/d(t)d(p)                  BINARY MASK: 00001000000 */
    double *dp2,  /* d2(v)/d(p)2                     BINARY MASK: 00010000000 */
    double *dxdt, /* d2(v)/d(x[])d(t)                BINARY MASK: 00100000000 */
    double *dxdp, /* d2(v)/d(x[])d(p)                BINARY MASK: 01000000000 */
    double *dw    /* (pointer to dw[]) d(v)/d(x[])   BINARY MASK: 10000000000 */
    )
{
    double *r = x;
    double s[NT];

    MTHREAD_ONCE(&initThreadBlock, threadInit);

    MTHREAD_MUTEX_LOCK(&global_data_mutex);
    order(FIRST, t, p, r,
                s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    if (mask & FIRST) {
        int i;
        *vmix = fillV (r, s, t, p);
        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            *vmix -= V(0);
            for (i=0; i<NR; i++) *vmix += r[i]*(V(0)-V(i+1));
        }
    }

    if(mask & SECOND) {
        double d2gdrds[NR][NT];
        double d2gdrdp[NR];
        double d2gds2[NT][NT];
        double d2gdsdp[NT];
        double dsdr[NT][NR], dsdp[NT];
        int i, j, k;

        fillD2GDRDS (r, s, t, p, d2gdrds);
        fillD2GDRDP (r, s, t, p, d2gdrdp);
        fillD2GDS2  (r, s, t, p, d2gds2);
        fillD2GDSDP (r, s, t, p, d2gdsdp);

        order(SECOND | FOURTH, t, p, r,
                    NULL, dsdr, NULL, dsdp, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            dx[i] = d2gdrdp[i];
            for (j=0; j<NT; j++) {
                dx[i] += d2gdrds[i][j]*dsdp[j] + d2gdsdp[j]*dsdr[j][i];
                for (k=0; k<NT; k++) dx[i] += d2gds2[j][k]*dsdp[j]*dsdr[k][i];
            }
        }

        if (returnMixingProperties) for (i=0; i<NR; i++) dx[i] += (V(0)-V(i+1));
    }

    if(mask & THIRD) {
        double d2gdrds[NR][NT];
        double d2gds2[NT][NT];
        double d2gdsdp[NT];
        double d3gdr2ds[NR][NR][NT];
        double d3gdr2dp[NR][NR];
        double d3gdrds2[NR][NT][NT];
        double d3gdrdsdp[NR][NT];
        double d3gds3[NT][NT][NT];
        double d3gds2dp[NT][NT];
        double dsdr[NT][NR], dsdp[NT], d2sdr2[NT][NR][NR], d2sdrdp[NT][NR];
        int i, j, k, l, m;

        fillD2GDRDS   (r, s, t, p, d2gdrds);
        fillD2GDS2    (r, s, t, p, d2gds2);
        fillD2GDSDP   (r, s, t, p, d2gdsdp);
        fillD3GDR2DS  (r, s, t, p, d3gdr2ds);
        fillD3GDR2DP  (r, s, t, p, d3gdr2dp);
        fillD3GDRDS2  (r, s, t, p, d3gdrds2);
        fillD3GDRDSDP (r, s, t, p, d3gdrdsdp);
        fillD3GDS3    (r, s, t, p, d3gds3);
        fillD3GDS2DP  (r, s, t, p, d3gds2dp);

        order(SECOND | FOURTH | FIFTH | SEVENTH, t, p, r,
                    NULL, dsdr, NULL, dsdp, d2sdr2, NULL, d2sdrdp,  NULL, NULL, NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            for (j=0; j<NR; j++) {
                dx2[i][j] = d3gdr2dp[i][j];
                for (k=0; k<NT; k++) {
                    dx2[i][j] += d3gdr2ds[i][j][k]*dsdp[k]
                     + d3gdrdsdp[i][k]*dsdr[k][j]
                     + d3gdrdsdp[j][k]*dsdr[k][i]
                     + d2gdsdp[k]*d2sdr2[k][i][j]
                     + d2gdrds[i][k]*d2sdrdp[k][j]
                     + d2gdrds[j][k]*d2sdrdp[k][i];
                    for (l=0; l<NT; l++) {
                        dx2[i][j] += d3gdrds2[i][k][l]*dsdr[k][j]*dsdp[l]
                       + d3gdrds2[j][k][l]*dsdr[k][i]*dsdp[l]
                       + d2gds2[k][l]*d2sdr2[k][i][j]*dsdp[l]
                       + d3gds2dp[k][l]*dsdr[k][i]*dsdr[l][j]
                       + d2gds2[k][l]*dsdr[k][i]*d2sdrdp[l][j]
                       + d2gds2[k][l]*dsdr[k][j]*d2sdrdp[l][i];
                        for (m=0; m<NT; m++)
                            dx2[i][j] += d3gds3[k][l][m]*dsdr[k][i]*dsdr[l][j]*dsdp[m];
                    }
                }
            }
        }

    }

    if(mask & FOURTH) {
        double d2gds2[NT][NT];
        double d2gdsdt[NT];
        double d2gdsdp[NT];
        double d2gdtdp = fillD2GDTDP (r, s, t, p);
        double dsdt[NT], dsdp[NT];
        int i, j;

        fillD2GDS2  (r, s, t, p, d2gds2);
        fillD2GDSDT (r, s, t, p, d2gdsdt);
        fillD2GDSDP (r, s, t, p, d2gdsdp);

        order(THIRD | FOURTH, t, p, r,
                    NULL, NULL, dsdt, dsdp, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        *dt = d2gdtdp;
        for (i=0; i<NT; i++) {
            *dt += d2gdsdt[i]*dsdp[i] + d2gdsdp[i]*dsdt[i];
            for (j=0; j<NT; j++) *dt += d2gds2[i][j]*dsdt[i]*dsdp[j];
        }

        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            *dt -= DVDT(0);
            for (i=0; i<NR; i++) *dt += r[i]*(DVDT(0)-DVDT(i+1));
        }
    }

    if(mask & FIFTH) {
        double d2gds2[NT][NT];
        double d2gdsdp[NT];
        double d2gdp2  = fillD2GDP2  (r, s, t, p);
        double dsdp[NT];
        int i,j;

        fillD2GDS2  (r, s, t, p, d2gds2);
        fillD2GDSDP (r, s, t, p, d2gdsdp);

        order(FOURTH, t, p, r,
                    NULL, NULL, NULL, dsdp, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        *dp = d2gdp2;
        for (i=0; i<NT; i++) {
            *dp += 2.0*d2gdsdp[i]*dsdp[i];
            for (j=0; j<NT; j++) *dp += d2gds2[i][j]*dsdp[i]*dsdp[j];
        }

        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            *dp -= DVDP(0);
            for (i=0; i<NR; i++) *dp += r[i]*(DVDP(0)-DVDP(i+1));
        }
    }

    if(mask & SIXTH) {
        double d2gds2[NT][NT];
        double d2gdsdt[NT];
        double d2gdsdp[NT];
        double d3gds3[NT][NT][NT];
        double d3gds2dp[NT][NT];
        double d3gds2dt[NT][NT];
        double d3gdsdtdp[NT];
        double d3gdsdt2[NT];
        double d3gdt2dp  = fillD3GDT2DP  (r, s, t, p);
        double dsdt[NT], dsdp[NT], d2sdt2[NT], d2sdtdp[NT];
        int i, j, k;

        fillD2GDS2    (r, s, t, p, d2gds2);
        fillD2GDSDT   (r, s, t, p, d2gdsdt);
        fillD2GDSDP   (r, s, t, p, d2gdsdp);
        fillD3GDS3    (r, s, t, p, d3gds3);
        fillD3GDS2DP  (r, s, t, p, d3gds2dp);
        fillD3GDS2DT  (r, s, t, p, d3gds2dt);
        fillD3GDSDTDP (r, s, t, p, d3gdsdtdp);
        fillD3GDSDT2  (r, s, t, p, d3gdsdt2);

        order(THIRD | FOURTH | EIGHTH | NINTH, t, p, r,
                    NULL, NULL, dsdt, dsdp, NULL, NULL, NULL, d2sdt2, d2sdtdp, NULL, NULL, NULL);

        *dt2 = d3gdt2dp;
        for (i=0; i<NT; i++) {
            *dt2 += d3gdsdt2[i]*dsdp[i] + 2.0*d2gdsdt[i]*d2sdtdp[i]
                        + d2gdsdp[i]*d2sdt2[i] + 2.0*d3gdsdtdp[i]*dsdt[i];
            for (j=0; j<NT; j++) {
                *dt2 += 2.0*d3gds2dt[i][j]*dsdt[i]*dsdp[j]
                            + d2gds2[i][j]*d2sdt2[i]*dsdp[j]
                            + 2.0*d2gds2[i][j]*dsdt[i]*d2sdtdp[j]
                            + d3gds2dp[i][j]*dsdt[i]*dsdt[j];
                for (k=0; k<NT; k++) *dt2 += d3gds3[i][j][k]*dsdt[i]*dsdt[j]*dsdp[k];
            }
        }

        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            *dt2 -= D2VDT2(0);
            for (i=0; i<NR; i++) *dt2 += r[i]*(D2VDT2(0)-D2VDT2(i+1));
        }
    }

    if(mask & SEVENTH) {
        double d2gds2[NT][NT];
        double d2gdsdt[NT];
        double d2gdsdp[NT];
        double d3gds3[NT][NT][NT];
        double d3gds2dt[NT][NT];
        double d3gds2dp[NT][NT];
        double d3gdsdtdp[NT];
        double d3gdsdp2[NT];
        double    d3gdtdp2  = fillD3GDTDP2  (r, s, t, p);
        double dsdt[NT], dsdp[NT], d2sdtdp[NT], d2sdp2[NT];
        int i, j, k;

        fillD2GDS2    (r, s, t, p, d2gds2);
        fillD2GDSDT   (r, s, t, p, d2gdsdt);
        fillD2GDSDP   (r, s, t, p, d2gdsdp);
        fillD3GDS3    (r, s, t, p, d3gds3);
        fillD3GDS2DT  (r, s, t, p, d3gds2dt);
        fillD3GDS2DP  (r, s, t, p, d3gds2dp);
        fillD3GDSDTDP (r, s, t, p, d3gdsdtdp);
        fillD3GDSDP2  (r, s, t, p, d3gdsdp2);

        order(THIRD | FOURTH | NINTH | TENTH, t, p, r,
                    NULL, NULL, dsdt, dsdp, NULL, NULL, NULL, NULL, d2sdtdp, d2sdp2, NULL, NULL);

        *dtdp = d3gdtdp2;
        for (i=0; i<NT; i++) {
            *dtdp += 2.0*d3gdsdtdp[i]*dsdp[i] + d2gdsdt[i]*d2sdp2[i]
             + 2.0*d2gdsdp[i]*d2sdtdp[i] + d3gdsdp2[i]*dsdt[i];
            for (j=0; j<NT; j++) {
                *dtdp += 2.0*d3gds2dp[i][j]*dsdt[i]*dsdp[j]
               + d2gds2[i][j]*dsdt[i]*d2sdp2[j]
               + 2.0*d2gds2[i][j]*d2sdtdp[i]*dsdp[j]
               + d3gds2dt[i][j]*dsdp[i]*dsdp[j];
                for (k=0; k<NT; k++) *dtdp += d3gds3[i][j][k]*dsdt[i]*dsdp[j]*dsdp[k];
            }
        }

        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            *dtdp -= D2VDTDP(0);
            for (i=0; i<NR; i++) *dtdp += r[i]*(D2VDTDP(0)-D2VDTDP(i+1));
        }
    }

    if(mask & EIGHTH) {
        double d2gds2[NT][NT];
        double d2gdsdp[NT];
        double d3gds3[NT][NT][NT];
        double d3gds2dp[NT][NT];
        double d3gdsdp2[NT];
        double    d3gdp3   = fillD3GDP3   (r, s, t, p);
        double dsdp[NT], d2sdp2[NT];
        int i, j, k;

        fillD2GDS2   (r, s, t, p, d2gds2);
        fillD2GDSDP  (r, s, t, p, d2gdsdp);
        fillD3GDS3   (r, s, t, p, d3gds3);
        fillD3GDS2DP (r, s, t, p, d3gds2dp);
        fillD3GDSDP2 (r, s, t, p, d3gdsdp2);

        order(FOURTH | TENTH, t, p, r,
                    NULL, NULL, NULL, dsdp, NULL, NULL, NULL, NULL, NULL, d2sdp2, NULL, NULL);

        *dp2 = d3gdp3;
        for (i=0; i<NT; i++) {
            *dp2 += 3.0*d3gdsdp2[i]*dsdp[i] + 3.0*d2gdsdp[i]*d2sdp2[i];
            for (j=0; j<NT; j++) {
                *dp2 += 3.0*d2gds2[i][j]*dsdp[i]*d2sdp2[j]
                            + 3.0*d3gds2dp[i][j]*dsdp[i]*dsdp[j];
                for (k=0; k<NT; k++) *dp2 += d3gds3[i][j][k]*dsdp[i]*dsdp[j]*dsdp[k];
            }
        }

        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            *dp2 -= D2VDP2(0);
            for (i=0; i<NR; i++) *dp2 += r[i]*(D2VDP2(0)-D2VDP2(i+1));
        }
    }

    if(mask & NINTH) {
        double d3gds3[NT][NT][NT];
        double d3gdrds2[NR][NT][NT];
        double d3gdrdsdt[NR][NT];
        double d3gds2dp[NT][NT];
        double d2gdrds[NR][NT];
        double d3gdrdtdp[NR];
        double d3gdsdtdp[NT];
        double d2gds2[NT][NT];
        double d2gdsdt[NT];
        double d3gdrdsdp[NR][NT];
        double d2gdsdp[NT];
        double d3gds2dt[NT][NT];
        double dsdt[NT], dsdp[NT], dsdr[NT][NR], d2sdrdt[NT][NR], d2sdrdp[NT][NR],
            d2sdtdp[NT];
        int i, j, k, l;

        fillD3GDS3    (r, s, t, p, d3gds3);
        fillD3GDRDS2  (r, s, t, p, d3gdrds2);
        fillD3GDRDSDT (r, s, t, p, d3gdrdsdt);
        fillD3GDS2DP  (r, s, t, p, d3gds2dp);
        fillD2GDRDS   (r, s, t, p, d2gdrds);
        fillD3GDRDTDP (r, s, t, p, d3gdrdtdp);
        fillD3GDSDTDP (r, s, t, p, d3gdsdtdp);
        fillD2GDS2    (r, s, t, p, d2gds2);
        fillD2GDSDT   (r, s, t, p, d2gdsdt);
        fillD3GDRDSDP (r, s, t, p, d3gdrdsdp);
        fillD2GDSDP   (r, s, t, p, d2gdsdp);
        fillD3GDS2DT  (r, s, t, p, d3gds2dt);

        order(SECOND | THIRD | FOURTH | SIXTH | SEVENTH | NINTH, t, p, r,
                    NULL, dsdr, dsdt, dsdp, NULL, d2sdrdt, d2sdrdp, NULL, d2sdtdp, NULL, NULL, NULL);

        for (i=0; i<NR; i++) {
            for (j=0,dxdt[i]=d3gdrdtdp[i]; j<NT; j++) {
                dxdt[i] += d3gdsdtdp[j]*dsdr[j][i] + d2gdsdt[j]*d2sdrdp[j][i] +
                   d3gdrdsdt[i][j]*dsdp[j] + d2gdrds[i][j]*d2sdtdp[j] +
                   d3gdrdsdp[i][j]*dsdt[j] + d2gdsdp[j]*d2sdrdt[j][i];
                for (k=0; k<NT; k++) {
                    dxdt[i] += d3gdrds2[i][j][k]*dsdt[j]*dsdp[k] +
                     d2gds2[j][k]*dsdt[j]*d2sdrdp[k][i] +
                     d2gds2[j][k]*dsdp[j]*d2sdrdt[k][i] +
                     d3gds2dt[j][k]*dsdr[j][i]*dsdp[k] +
                     d3gds2dp[j][k]*dsdr[j][i]*dsdt[k] +
                     d2gds2[j][k]*dsdr[j][i]*d2sdtdp[k];
                    for (l=0; l<NT; l++)
                        dxdt[i] += d3gds3[j][k][l]*dsdr[j][i]*dsdt[k]*dsdp[l];
                }
            }
        }

        if (returnMixingProperties) for (i=0; i<NR; i++) dxdt[i] += (DVDT(0)-DVDT(i+1));
    }

    if(mask & TENTH) {
        double d3gds3[NT][NT][NT];
        double d3gdrds2[NR][NT][NT];
        double d3gds2dp[NT][NT];
        double d2gdrds[NR][NT];
        double d3gdrdsdp[NR][NT];
        double d3gdrdp2[NR];
        double d3gdsdp2[NT];
        double d2gdsdp[NT];
        double d2gds2[NT][NT];
        double dsdr[NT][NR], dsdp[NT], d2sdrdp[NT][NR], d2sdp2[NT];
        int i, j, k, l;

        fillD3GDS3    (r, s, t, p, d3gds3);
        fillD3GDRDS2  (r, s, t, p, d3gdrds2);
        fillD3GDS2DP  (r, s, t, p, d3gds2dp);
        fillD2GDRDS   (r, s, t, p, d2gdrds);
        fillD3GDRDSDP (r, s, t, p, d3gdrdsdp);
        fillD3GDRDP2  (r, s, t, p, d3gdrdp2);
        fillD3GDSDP2  (r, s, t, p, d3gdsdp2);
        fillD2GDSDP   (r, s, t, p, d2gdsdp);
        fillD2GDS2    (r, s, t, p, d2gds2);

        order(SECOND | FOURTH | SEVENTH | TENTH, t, p, r,
                    NULL, dsdr, NULL, dsdp, NULL, NULL, d2sdrdp, NULL, NULL, d2sdp2, NULL, NULL);

        for (i=0; i<NR; i++) {
            for (j=0,dxdp[i]=d3gdrdp2[i]; j<NT; j++) {
                dxdp[i] += d3gdsdp2[j]*dsdr[j][i] + d2gdsdp[j]*d2sdrdp[j][i] +
                   2.0*d3gdrdsdp[i][j]*dsdp[j] + d2gdrds[i][j]*d2sdp2[j] +
                   d2gdsdp[j]*d2sdrdp[j][i];
                for (k=0; k<NT; k++) {
                    dxdp[i] += d3gdrds2[i][j][k]*dsdp[j]*dsdp[k] +
                     2.0*d2gds2[j][k]*dsdp[j]*d2sdrdp[k][i] +
                     2.0*d3gds2dp[j][k]*dsdr[j][i]*dsdp[k] +
                     d2gds2[j][k]*dsdr[j][i]*d2sdp2[k];
                    for (l=0; l<NT; l++)
                        dxdp[i] += d3gds3[j][k][l]*dsdr[j][i]*dsdp[k]*dsdp[l];
                }
            }
        }

        if (returnMixingProperties) for (i=0; i<NR; i++) dxdp[i] += (DVDP(0)-DVDP(i+1));
    }

    if(mask & ELEVENTH) {
        double d2gdsdw[NT][3*NP];
        double d2gdpdw[3*NP];
        double d2gds2[NT][NT];
        double d2gdsdp[NT];
        double dsdw[NT][3*NP], dsdp[NT];
        int i, k, l;

        fillD2GDSDW (r, s, t, p, d2gdsdw);
        fillD2GDPDW (r, s, t, p, d2gdpdw);
        fillD2GDS2  (r, s, t, p, d2gds2);
        fillD2GDSDP (r, s, t, p, d2gdsdp);

        order(FOURTH | ELEVENTH, t, p, r,
                    NULL, NULL, NULL, dsdp, NULL, NULL, NULL, NULL, NULL, NULL, dsdw, NULL);

        for (i=0; i<NP; i++) {
            if (meltsAndCO2_H2OModelParameters[i].activeH) {
                dw[i] = d2gdpdw[i];
                for (k=0; k<NT; k++) {
                    dw[i] += d2gdsdw[k][i]*dsdp[k] + d2gdsdp[k]*dsdw[k][i];
                    for (l=0; l<NT; l++) dw[i] += d2gds2[k][l]*dsdp[k]*dsdw[l][i] ;
                }
                dw[i] *= -1.0;
            } else dw[i] = 0.0;

            if (meltsAndCO2_H2OModelParameters[i].activeS) {
                dw[NP+i] = d2gdpdw[NP+i];
                for (k=0; k<NT; k++) {
                    dw[NP+i] += d2gdsdw[k][NP+i]*dsdp[k] + d2gdsdp[k]*dsdw[k][NP+i];
                    for (l=0; l<NT; l++) dw[NP+i] += d2gds2[k][l]*dsdp[k]*dsdw[l][NP+i] ;
                }
                dw[NP+i] *= -1.0;
            } else dw[NP+i] = 0.0;

            if (meltsAndCO2_H2OModelParameters[i].activeV) {
                dw[2*NP+i] = d2gdpdw[2*NP+i];
                for (k=0; k<NT; k++) {
                    dw[2*NP+i] += d2gdsdw[k][2*NP+i]*dsdp[k] + d2gdsdp[k]*dsdw[k][2*NP+i];
                    for (l=0; l<NT; l++) dw[2*NP+i] += d2gds2[k][l]*dsdp[k]*dsdw[l][2*NP+i] ;
                }
                dw[2*NP+i] *= -1.0;
            } else dw[2*NP+i] = 0.0;
        }

        if (returnMixingProperties) { /* Convert Solution Properties -> Mixing Properties */
            if (meltsAndCO2_H2OModelParameters[NW+0].activeS) dw[NP+NW+0] -= 1.0;
            for (i=0; i<NR; i++) {
                if (meltsAndCO2_H2OModelParameters[NW  +0].activeS) dw[NP+NW  +0] += r[i];
                if (meltsAndCO2_H2OModelParameters[NW+i+1].activeS) dw[NP+NW+i+1] -= r[i];
            }
        }

        for (i=0; i<(3*NP); i++) if (fabs(dw[i]) < 10.0*DBL_EPSILON) dw[i] = 0.0;
    }

    MTHREAD_MUTEX_UNLOCK(&global_data_mutex);
}

/* ============================================================================
   In the following public routine:
   m = m[i] = moles of the ith component in the liquid, and
   mu*O2    = mu O2 - mu0 O2, defined from the vector o[]
   This routine should be deprecated as it is inconsistent with conLiq_CO2_H2O
   ==========================================================================*/

void
muO2Liq_CO2_H2O(int mask, double t, double p, double *m,
    double *muO2, /* muO2      = mu*O2                 BINARY MASK: 0000000001 */
    double *dm,   /* dm[i]     = d mu*O2/dm[i]         BINARY MASK: 0000000010 */
    double *dt,   /* dt        = d mu*O2/d T           BINARY MASK: 0000000100 */
    double *dp,   /* dp        = d mu*O2/d P           BINARY MASK: 0000001000 */
    double **d2m, /* d2m[i][j] = d mu*O2/dm[i][j]      BINARY MASK: 0000010000 */
    double *d2mt, /* d2mt[i]   = d mu*O2/dm[i]dt       BINARY MASK: 0000100000 */
    double *d2mp, /* d2mp[i]   = d mu*O2/dm[i]dp       BINARY MASK: 0001000000 */
    double *d2t2, /* d2t2      = d mu*O2/dt2           BINARY MASK: 0010000000 */
    double *d2tp, /* d2tp      = d mu*O2/dtdp          BINARY MASK: 0100000000 */
    double *d2p2) /* d2p2      = d mu*O2/dp2           BINARY MASK: 1000000000 */
{
    /*-------------------------------------------------------------------------
        The algorithm used is that given by:
        Kress, VC and Carmichael, ISE (1991) The compressibility of silicate
            liquids containing Fe2O3 and the effect of composition, temperature,
            oxygen fugacity and pressure on their redox states.
            Contributions to Mineralogy and Petrology (in press)
        Coefficients for the oxides are initialized in LIQ_STRUCT_DATA.H
    --------------------------------------------------------------------------*/
    static int indexFeO = -1, indexFe2O3 = -1;
    static const double t0 = 1673.15,                      /* K       */
                       a =    0.196,
                       b =    1.1492e4,                  /* K       */
                       c =   -6.675,
                       e =   -3.364,
                       f =   -7.01e-7  * 1.0e5,          /* K/bar   */
                       g =   -1.54e-10 * 1.0e5,          /* 1/bar   */
                       h =    3.85e-17 * 1.0e5 * 1.0e5;  /* K/bar^2 */
    double mOx[NA], vTemp[NA], mTemp[NA][NA];
    int i, j;
    double sum;

    MTHREAD_ONCE(&initThreadBlock, threadInit);

    MTHREAD_MUTEX_LOCK(&global_data_mutex);
    if ((indexFeO == -1) || (indexFe2O3 == -1)) {
        for (i=0; i<nc; i++) {
            if (bulkSystem[i].type == FEO)   indexFeO   = i;
            if (bulkSystem[i].type == FE2O3) indexFe2O3 = i;
        }
        if (indexFeO == 0 || indexFe2O3 == 0) {
       printf("Fatal error in muO2Liq_CO2_H2O (LIQUID.C)\n");
       printf("The oxides FeO and Fe2O3 cannot be identified.\n");
       return;
        }
    }
    MTHREAD_MUTEX_UNLOCK(&global_data_mutex);

    for (i=0; i<nc; i++)
        for (j=0, mOx[i]=0.0; j<nlc; j++) mOx[i] += (liquid[j].liqToOx)[i]*m[j];

    for (i=0, sum=0.0; i<nc; i++) { sum += mOx[i]; } sum += mOx[indexFe2O3];
    if (sum == 0.0 || mOx[indexFeO] == 0.0 || mOx[indexFe2O3] == 0.0) {
        if (mask & FIRST)   *muO2 = 0.0;
        if (mask & SECOND)  for (i=0; i<nlc; i++) dm[i] = 0.0;
        if (mask & THIRD)   *dt = 0.0;
        if (mask & FOURTH)  *dp = 0.0;
        if (mask & FIFTH)   {
            for (i=0; i<nlc; i++) for (j=0; j<nlc; j++) d2m[i][j] = 0.0;
        }
        if (mask & SIXTH)   for (i=0; i<nlc; i++) d2mt[i] = 0.0;
        if (mask & SEVENTH) for (i=0; i<nlc; i++) d2mp[i] = 0.0;
        if (mask & EIGHTH)  *d2t2 = 0.0;
        if (mask & NINTH)   *d2tp = 0.0;
        if (mask & TENTH)   *d2p2 = 0.0;
        return;
    }

    /*-------------------------------------------------------------------------*/

    if (mask & FIRST) {
        double temp;
        if (p < 50000.0) temp = b/t + c + e*(1.0 - t0/t - log(t/t0)) + f*p/t + g*(t-t0)*p/t + h*SQUARE(p)/t;
        else             temp = b/t + c + e*(1.0 - t0/t - log(t/t0)) + f*50000.0/t + g*(t-t0)*50000.0/t + h*SQUARE(50000.0)/t
                                                    - a*log(10.0)*(608.966*p/10000.0-608.966*5.0)/t;
        for (i=0; i<nc; i++) temp += bulkSystem[i].coeff*mOx[i]/sum;
        temp += 2.0*bulkSystem[indexFeO].coeff*mOx[indexFe2O3]/sum - bulkSystem[indexFe2O3].coeff*mOx[indexFe2O3]/sum;
        *muO2 = R*t*(log(mOx[indexFe2O3]/mOx[indexFeO]) - temp)/a;
    }

    if (mask & SECOND) {
        for (j=0; j<nc; j++) {
            double factor = (j == indexFe2O3) ? 2.0 : 1.0;
            for (i=0, vTemp[j]=0.0; i<nc; i++)
                vTemp[j] -= (i == j) ?
                    bulkSystem[i].coeff*(1.0-factor*mOx[i]/sum)/sum :
                    - bulkSystem[i].coeff*factor*mOx[i]/SQUARE(sum);
            vTemp[j] += - (factor*mOx[indexFe2O3]/SQUARE(sum))
                *(bulkSystem[indexFe2O3].coeff-2.0*bulkSystem[indexFeO].coeff);
            if      (j == indexFeO) vTemp[j] += -1.0/mOx[indexFeO];
            else if (j == indexFe2O3) {
                vTemp[j] += 1.0/mOx[indexFe2O3];
                vTemp[j] += (1.0/sum)
                    *(bulkSystem[indexFe2O3].coeff-2.0*bulkSystem[indexFeO].coeff);
            }
            vTemp[j] *= R*t/a;
        }
        for (i=0; i<nlc; i++)
            for (j=0, dm[i]=0.0; j<nc; j++) dm[i] += vTemp[j]*(liquid[i].liqToOx)[j];
    }

    if (mask & THIRD) {
        double temp;
        if (p < 50000.0) temp = b/t + c + e*(1.0 - t0/t - log(t/t0)) + f*p/t + g*(t-t0)*p/t + h*SQUARE(p)/t;
        else             temp = b/t + c + e*(1.0 - t0/t - log(t/t0)) + f*50000.0/t + g*(t-t0)*50000.0/t + h*SQUARE(50000.0)/t
                                                    - a*log(10.0)*(608.966*p/10000.0-608.966*5.0)/t;
        for (i=0; i<nc; i++) temp += bulkSystem[i].coeff*mOx[i]/sum;
        temp += 2.0*bulkSystem[indexFeO].coeff*mOx[indexFe2O3]/sum - bulkSystem[indexFe2O3].coeff*mOx[indexFe2O3]/sum;
        *dt = R*(log(mOx[indexFe2O3]/mOx[indexFeO]) - temp)/a + R*t*(b/SQUARE(t) - e*(t0/t-1.0)*(1.0/t))/a;
        if (p < 50000.0) *dt += R*t*(f*p/SQUARE(t) - g*(t0/t)*(p/t) + h*SQUARE(p/t))/a;
        else             *dt += R*t*(f*50000.0/SQUARE(t) - g*(t0/t)*(50000.0/t) + h*SQUARE(50000.0/t)
                                 - a*log(10.0)*(608.966*p/10000.0-608.966*5.0)/SQUARE(t))/a;
    }

    if (mask & FOURTH) {
        if (p < 50000.0) *dp = R*t*(-f/t - g*(t-t0)/t - 2.0*h*p/t)/a;
        else             *dp = R*t*(a*log(10.0)*608.966/10000.0/t)/a;
    }

    if (mask & FIFTH) {
        int k, l;
        for (k=0; k<nc; k++) {
            double factorK = (k == indexFe2O3) ? 2.0 : 1.0;
            for (j=0; j<nc; j++) {
                double factorJ = (j == indexFe2O3) ? 2.0 : 1.0;
                for (i=0, mTemp[k][j]=0.0; i<nc; i++) {
                    mTemp[k][j] -=
                        2.0*factorJ*factorK*bulkSystem[i].coeff*mOx[i]/CUBE(sum);
                    if (i == j) mTemp[k][j] -= - factorK*bulkSystem[i].coeff/SQUARE(sum);
                    if (i == k) mTemp[k][j] -= - factorJ*bulkSystem[i].coeff/SQUARE(sum);
                }
                mTemp[k][j] += 2.0*(factorJ*factorK*mOx[indexFe2O3]/CUBE(sum))
                    *(bulkSystem[indexFe2O3].coeff-2.0*bulkSystem[indexFeO].coeff);
                if      (j == indexFeO && k == indexFeO)
                    mTemp[k][j] += 1.0/SQUARE(mOx[indexFeO]);
                else if (j == indexFe2O3 && k == indexFe2O3) {
                    mTemp[k][j] += -1.0/SQUARE(mOx[indexFe2O3]);
                    mTemp[k][j] += -((factorJ+factorK)/SQUARE(sum))
                        *(bulkSystem[indexFe2O3].coeff-2.0*bulkSystem[indexFeO].coeff);
                }
                else if (j == indexFe2O3) mTemp[k][j] += -(factorK/SQUARE(sum))
                        *(bulkSystem[indexFe2O3].coeff-2.0*bulkSystem[indexFeO].coeff);
                else if (k == indexFe2O3) mTemp[k][j] += -(factorJ/SQUARE(sum))
                        *(bulkSystem[indexFe2O3].coeff-2.0*bulkSystem[indexFeO].coeff);
                mTemp[k][j] *= R*t/a;
            }
        }
        for (i=0; i<nlc; i++) for (j=0; j<nlc; j++)
            for (k=0, d2m[i][j]=0.0; k<nc; k++) for (l=0; l<nc; l++)
                d2m[i][j] += mTemp[k][l]*(liquid[i].liqToOx)[k]*(liquid[j].liqToOx)[l];
    }

    if (mask & SIXTH) {
        for (j=0; j<nc; j++) {
            double factor = (j == indexFe2O3) ? 2.0 : 1.0;
            for (i=0, vTemp[j]=0.0; i<nc; i++)
                vTemp[j] -= (i == j) ?
                    bulkSystem[i].coeff*(1.0-factor*mOx[i]/sum)/sum :
                    - bulkSystem[i].coeff*factor*mOx[i]/SQUARE(sum);
            vTemp[j] += - (factor*mOx[indexFe2O3]/SQUARE(sum))
                *(bulkSystem[indexFe2O3].coeff-2.0*bulkSystem[indexFeO].coeff);
            if      (j == indexFeO) vTemp[j] += -1.0/mOx[indexFeO];
            else if (j == indexFe2O3) {
                vTemp[j] += 1.0/mOx[indexFe2O3];
                vTemp[j] += (1.0/sum)
                    *(bulkSystem[indexFe2O3].coeff-2.0*bulkSystem[indexFeO].coeff);
            }
            vTemp[j] *= R/a;
        }
        for (i=0; i<nlc; i++) for (j=0, d2mt[i]=0.0; j<nc; j++)
            d2mt[i] += vTemp[j]*(liquid[i].liqToOx)[j];
    }

    if (mask & SEVENTH) {
        for (i=0; i<nlc; i++) d2mp[i] = 0.0;
    }

    if (mask & EIGHTH) {

        *d2t2 = 2.0*R*(b/SQUARE(t) - e*(t0/t-1.0)*(1.0/t))/a
                    + R*t*(-2.0*b/CUBE(t) - e*(1.0 - 2.0*t0/t)/SQUARE(t))/a;
    }

    if (mask & NINTH) {
        if (p < 50000.0) *d2tp = -R*(f/t + g*(t-t0)/t + 2.0*h*p/t)/a + R*t*(f/SQUARE(t) - g*(t0/t)/t + 2.0*h*p/SQUARE(t))/a;
        else             *d2tp = -R*(-a*log(10.0)*608.966/10000.0/t)/a + R*t*(-a*log(10.0)*608.966/10000.0/SQUARE(t))/a;
    }

    if (mask & TENTH) {
        if (p < 50000.0) *d2p2 = - (R*t*2.0*h/t)/a;
        else             *d2p2 = 0.0;
    }

}

void
visLiq_CO2_H2O(int mask, double t, double p, double *r,
    double *viscosity  /* log(10) viscosity            BINARY MASK: 00000001 */
    )
{
    double coeff[NA], factor[NA], m[NA], x[NA], sum;
    int nSiO2 = -1, i, j;

    MTHREAD_ONCE(&initThreadBlock, threadInit);

    struct _shawModel {
        char   *oxide;
        double coeff;
        double factor;
    } shawModel[] = {
        { "TiO2",    4.5, 1.0 }, { "Al2O3",  6.7, 2.0 },
        { "Fe2O3",  3.4, 2.0 }, { "FeO",    3.4, 1.0 },
        { "MgO",    3.4, 1.0 }, { "CaO",    4.5, 1.0 },
        { "Na2O",    2.8, 1.0 }, { "K2O",    2.8, 1.0 },
        { "H2O",    2.0, 1.0 }
    };
    const int nShaw = (sizeof shawModel / sizeof(struct _shawModel));

    for (j=0; j<NA; j++) { coeff[j] = 0.0; factor[j] = 0.0; }
    for (i=0; i<nShaw; i++) {
        for (j=0; j<NA; j++) if (strcmp(shawModel[i].oxide, bulkSystem[j].label) == 0) {
            coeff[j]  = shawModel[i].coeff;
            factor[j] = shawModel[i].factor;
            break;
        }
    }
    for (i=0; i<NA; i++) if (strcmp("SiO2", bulkSystem[i].label) == 0) { nSiO2 = i; break; }

    if (nSiO2 == -1) { *viscosity = 0.0; return; }

    /* m[0] --> m[NA-1] is an array of mole fractions of liquid components      */
    for (m[0]=1.0, i=0; i<NR; i++) { m[0] -= r[i]; m[i+1] = r[i]; }

    /* convert m[] -> x[] : mole fractions of liquid comp -> moles of oxides    */
    for (i=0; i<NA; i++)
        for (j=0, x[i]=0.0; j<NA; j++) x[i] += (liquid[j].liqToOx)[i]*m[j];

    /* Convert to the Shaw mole fractions                                       */
    for (i=0, sum=0.0; i<NA; i++) {
        if (factor[i] > 0.0) x[i] *= factor[i];
        sum += x[i];
    }
    for (i=0; i<NA; i++) x[i] /= (sum != 0.0) ? sum : 1.0;

    if (mask & FIRST) {
        for (i=0, *viscosity=0.0; i<NA; i++) *viscosity += coeff[i]*x[nSiO2]*x[i];
        *viscosity /= (x[nSiO2] < 1.0) ? 1.0 - x[nSiO2] : 1.0;
        *viscosity  = (*viscosity)*(10000.0/t - 1.50)  - 6.40;
        *viscosity /= log(10.0);
    }
}

/* end of file LIQUID.C */
