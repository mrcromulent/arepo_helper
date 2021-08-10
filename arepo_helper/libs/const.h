#ifndef AREPO_HELPER_LIBS_CONST_H
#define AREPO_HELPER_LIBS_CONST_H

#include <cstddef>

extern const double msol;
extern const double G;
extern const double electron_charge_esu;

#define NAMES \
X(ACCELERATION, "Acceleration") \
X(COORDINATES, "Coordinates") \
X(GRAVITYINTERACTIONS, "GravityInteractions") \
X(MASSES, "Masses") \
X(PARTICLEIDS, "ParticleIDs") \
X(POTENTIAL, "Potential") \
X(SOFTENINGS, "Softenings") \
X(TIMESTEP, "TimeStep") \
X(VELOCITIES, "Velocities") \
X(GFMMETALLICITY, "GFM_Metallicity") \
X(GFMMETALS, "GFM_Metals") \
X(METALLICITY, "Metallicity") \
X(ALLOWREFINEMENT, "AllowRefinement") \
X(BFIELDGRADIENT, "BfieldGradient") \
X(CELLSPIN, "CellSpin") \
X(CENTEROFMASS, "CenterOfMass") \
X(CHEMICALABUNDANCES, "ChemicalAbundances") \
X(COOLTIME, "CoolTime") \
X(COSMICRAYSPECIFICENERGY, "CosmicRaySpecificEnergy") \
X(COSMICRAYSTREAMINGREGULARIZATION, "CosmicRayStreamingRegularization") \
X(CRPRESSUREGRADIENT, "CRPressureGradient") \
X(CURLB, "CurlB") \
X(DENSITY, "Density") \
X(DENSITYGRADIENT, "DensityGradient") \
X(DIVBCLEENING, "DivBCleening") \
X(ELECTRONABUNDANCE, "ElectronAbundance") \
X(ENERGYDISSIPATION, "EnergyDissipation") \
X(GFMAGNRADIATION, "GFM_AGNRadiation") \
X(GFMCOOLINGRATE, "GFM_CoolingRate") \
X(HIGHRESGASMASS, "HighResGasMass") \
X(INTERNALENERGY, "InternalEnergy") \
X(MACHNUMBER, "Machnumber") \
X(MAGNETICFIELD, "MagneticField") \
X(MAGNETICFIELDDIVERGENCE, "MagneticFieldDivergence") \
X(MAGNETICVECTORPOTENTIAL, "MagneticVectorPotential") \
X(MOLECULARWEIGHT, "MolecularWeight") \
X(NEUTRALHYDROGENABUNDANCE, "NeutralHydrogenAbundance") \
X(NUCLEARCOMPOSITION, "NuclearComposition") \
X(NUCLEARENERGYGENERATIONRATE, "NuclearEnergyGenerationRate") \
X(PASSIVESCALARS, "PassiveScalars") \
X(PRESSUREGRADIENT, "PressureGradient") \
X(PRESSURE, "Pressure") \
X(RATEOFCHANGEOFMAGNETICFIELD, "RateOfChangeOfMagneticField") \
X(SFPROBABILITY, "SFProbability") \
X(SMOOTHEDMAGNETICFIELD, "SmoothedMagneticField") \
X(SOUNDSPEED, "SoundSpeed") \
X(STARFORMATIONRATE, "StarFormationRate") \
X(TEMPERATURE, "Temperature") \
X(TURBULENTENERGY, "TurbulentEnergy") \
X(VELOCITYDIVERGENCE, "VelocityDivergence") \
X(VELOCITYGRADIENT, "VelocityGradient") \
X(VERTEXVELOCITY, "VertexVelocity") \
X(VOLUME, "Volume") \
X(VORTICITY, "Vorticity") \
X(RADIUS, "Radius") \
X(MR, "Mr") \
X(DM, "Dm") \
X(DEDT, "dedT") \
X(DPDT, "dpdT") \
X(GAMMA1, "gamma1") \
X(GAMMA2, "gamma2") \
X(GAMMA3, "gamma3") \
X(CV, "cv") \
X(CP, "cp")

extern const char* const f[];

#define X(a, b) a,
enum N : size_t
{
    NAMES
};
#undef X

#endif //AREPO_HELPER_LIBS_CONST_H
