#ifndef HELM_EOS_H
#define HELM_EOS_H

#define IMAX 271
#define JMAX 101

typedef double helm_eos_table_entry[IMAX][JMAX];

typedef struct {
    // number of entries
    int ntemp, nrho;

    // density and temperature ranges
    double ltempMin, ltempMax;
    double lrhoMin, lrhoMax;
    double ltempDelta;
    double lrhoDelta;
    double tempDelta;
    double rhoDelta;
    double tempMin, tempMax;
    double rhoMin, rhoMax;
    double temp[JMAX];
    double rho[IMAX];


    // d means derivative w. r. t. density; t means derivative w. r. t. temperature

    // Helmholtz free energy
    helm_eos_table_entry f;
    helm_eos_table_entry fd;
    helm_eos_table_entry ft;
    helm_eos_table_entry fdd;
    helm_eos_table_entry ftt;
    helm_eos_table_entry fdt;
    helm_eos_table_entry fddt;
    helm_eos_table_entry fdtt;
    helm_eos_table_entry fddtt;

    // pressure derivative w. r. t. density
    helm_eos_table_entry dpdf;
    helm_eos_table_entry dpdfd;
    helm_eos_table_entry dpdft;
    helm_eos_table_entry dpdfdd;
    helm_eos_table_entry dpdftt;
    helm_eos_table_entry dpdfdt;

    // chemical potential
    helm_eos_table_entry ef;
    helm_eos_table_entry efd;
    helm_eos_table_entry eft;
    helm_eos_table_entry efdd;
    helm_eos_table_entry eftt;
    helm_eos_table_entry efdt;

    // number density
    helm_eos_table_entry xf;
    helm_eos_table_entry xfd;
    helm_eos_table_entry xft;
    helm_eos_table_entry xfdd;
    helm_eos_table_entry xftt;
    helm_eos_table_entry xfdt;


    // species information
    int nspecies;
    double *na;
    double *nai;
    double *nz;

    int useCoulombCorrections;
    int verbose;
} t_helm_eos_table;

struct eos_value {
    double v; /* value */
    double drho; /* derivative with density */
    double dtemp; /* derivative with temperature */
    double dabar; /* derivative with abar */
    double dzbar; /* derivative with zbar */
};

struct eos_result {
    double temp;
    struct eos_value p; /*  pressure */
    struct eos_value e; /*  specific energy */
    struct eos_value s; /*  specific entropy */
    struct eos_value etaele; /*  degeneracy parameter: electron chemical potential / (k_B * T) */
    struct eos_value nep; /*  electron + positron number density */
    double cv, cp; /*  specific heat at constant volume, pressure */
    double chit, chid; /*  temperature and density exponents from Cox & Giuli */
    double gamma_1, gamma_2, gamma_3; /*  gammas from Cox & Giuli */
    double nabla_ad; /*  nabla adiabatic */
    double sound; /*  relativistic speed of sound */

    double abar, zbar; /*  mean mass number and charge number */
};

/* cached values that only depend on composition and density
 * huge time saver for the egiven iterations
 */
struct helm_eos_cache {
    double abar, zbar; /* mean nuclear mass and charge */
    double ye; /* zbar / abar */
    double din, ldin; /* rho * ye */
    double ytot; /* 1.0 / abar; */
    double xni; /* GSL_CONST_NUM_AVOGADRO * ytot * rho */
    double dxnidd; /* GSL_CONST_NUM_AVOGADRO * ytot */
    double dxnida; /* - xni * ytot */
    double lswot15; /* 1.5 * log((2.0 * M_PI * GSL_CONST_CGS_UNIFIED_ATOMIC_MASS * GSL_CONST_CGS_BOLTZMANN) / (GSL_CONST_CGS_PLANCKS_CONSTANT_H * GSL_CONST_CGS_PLANCKS_CONSTANT_H) * temp) without the temp factor */
    double ywot; /* log(abar * abar * sqrt(abar) / (rho * GSL_CONST_NUM_AVOGADRO)) */
};


int eos_init( t_helm_eos_table* helm_eos_table, const char* datafile, const char* speciesfile, int useCoulombCorrections, int verbose );
void eos_deinit( t_helm_eos_table* helm_eos_table );

int eos_calc_tgiven( t_helm_eos_table* helm_eos_table, double rho, const double xnuc[], double temp, struct eos_result *res );
int eos_calc_tgiven_onlye( t_helm_eos_table* helm_eos_table, double rho, const double xnuc[], double temp, struct eos_result *res );
int eos_calc_tgiven_azbar( t_helm_eos_table* helm_eos_table, double rho, const struct helm_eos_cache *cache, double temp, struct eos_result *res, int only_e );
int eos_calc_egiven( t_helm_eos_table* helm_eos_table, double rho, const double xnuc[], double e, double *tempguess, struct eos_result *res );
int eos_calc_pgiven(t_helm_eos_table* helm_eos_table, double rho, const double xnuc[], double p, double *tempguess, struct eos_result *res);
int eos_calc_egiven_y( t_helm_eos_table* helm_eos_table, double rho, const double y[], double e, double *tempguess, struct eos_result *res );
int eos_calc_ptgiven( t_helm_eos_table* helm_eos_table, double p, const double xnuc[], double temp, double *rho, struct eos_result *res );
int helm_eos_update_cache(double rho, double abar, double zbar, struct helm_eos_cache *cache );

#endif /* HELM_EOS_H */
