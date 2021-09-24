#ifndef AREPO_HELPER_LIBS_EOS_H
#define AREPO_HELPER_LIBS_EOS_H

#define EOS_MAXITER 40

/* basic constants in cgs units */
#define EOS_EPS 1.0e-13
#define PROTONMASS 1.6726e-24
#define BOLTZMANN 1.3806e-16
#define RADCONST 7.65677e-15 /* 4 * Stefan Boltzmann constant / speed of light */
#define ELECTRONCHARGE 4.80320427e-10
#define AVOGADRO 6.0222e23

#define EOS_A1 1.244e15
#define EOS_A2 7.56e-15
#define EOS_A3 3.024e-14

typedef struct {
    int nspecies;
    double *nuclearmasses;
    double *nuclearcharges;
    char *datafile;
    int ntemp, nrho, nye;
    double tempMin, tempMax, ltempMin, ltempMax, ltempDelta, ltempDeltaI;
    double rhoMin, rhoMax, lrhoMin, lrhoMax, lrhoDelta, lrhoDeltaI;
    double yeMin, yeMax, yeDelta, yeDeltaI;
    double *ltemp, *lrho, *ye;

    double *p;           /* pressure */
    double *e;           /* energy per mass */
    double *dedt;        /* derivative energy with temperature */
    double *dpdt;        /* derivative pressure with temperature */
    double *dpdr;        /* derivative pressure with density */
} t_eos_table;

int eos_init( t_eos_table* eos_table, char* datafile, char* speciesfile );
int eos_init_old( t_eos_table* eos_table, char* datafile, char* speciesfile );
void eos_deinit( t_eos_table* eos_table );
int eos_calc_egiven( t_eos_table* eos_table, double rho, const double *xnuc, double e, double *temp, double *p, double *dpdr );
int eos_calc_tgiven( t_eos_table* eos_table, double rho, const double *xnuc, double temp, double *e, double *dedt, double *p );
void eos_trilinear_e( t_eos_table* eos_table, double temp, double rho, double ye, double *e, double *dedt, double *p );
void eos_trilinear( t_eos_table* eos_table, double temp, double rho, double ye, double *e, double *dedt, double *p, double *dpdr, double *dpdt );
void eos_radiation( double temp, double rho, double *e, double *dedt, double *p, double *dpdt );
void eos_coulCorr( t_eos_table* eos_table, double temp, double rho, double *n, double ni, double ne, double zeff, double *e, double *dedt, double *p, double *dpdr, double *dpdt );

double eos_SwapDouble( double Val );
int eos_SwapInt( int Val );
void eos_checkswap( char* fname, int *swap );

#endif //AREPO_HELPER_LIBS_EOS_H
