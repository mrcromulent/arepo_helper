from collections import OrderedDict


class BaseEnum:
    """Base class for the Enums used in the project."""

    @classmethod
    def names(cls):
        """Returns a list of user-specified attributes names in the class.

        :return: List of attributes
        :rtype: list
        """
        return [i for i in cls.__dict__.keys() if i[:1] != "_"]

    @classmethod
    def values(cls):
        """Returns a list of user-specified attribute values in the class.

        :return: List of attribute values
        :rtype: list
        """
        return [cls.__dict__[i] for i in cls.__dict__.keys() if i[:1] != "_"]


class ArepoHeader(BaseEnum):
    """Fields which appear in the Snapshot/IC header"""

    BOXSIZE = "BoxSize"
    COMPOSITIONVECTORLENGTH = "Composition_vector_length"
    FLAGCOOLING = "Flag_Cooling"
    FLAGDOUBLEPRECISION = "Flag_DoublePrecision"
    FLAGENTROPYICS = "Flag_Entropy_ICs"
    FLAGFEEDBACK = "Flag_Feedback"
    FLAGICINFO = "Flag_IC_info"
    FLAGMETALS = "Flag_Metals"
    FLAGSFR = "Flag_Sfr"
    FLAGSTELLARAGE = "Flag_StellarAge"
    GITCOMMIT = "Git_commit"
    GITDATE = "Git_date"
    HUBBLEPARAM = "HubbleParam"
    MASSTABLE = "MassTable"
    NGROUPSTHISFILE = "Ngroups_ThisFile"
    NGROUPSTOTAL = "Ngroups_Total"
    NIDSTHISFILE = "Nids_ThisFile"
    NIDSTOTAL = "Nids_Total"
    NSUBGROUPSTHISFILE = "Nsubgroups_ThisFile"
    NSUBGROUPSTOTAL = "Nsubgroups_Total"
    NUMFILES = "NumFiles"
    NUMFILESPERSNAPSHOT = "NumFilesPerSnapshot"
    NUMPARTTHISFILE = "NumPart_ThisFile"
    NUMPARTTOTAL = "NumPart_Total"
    NUMPARTTOTALHIGHWORD = "NumPart_Total_HighWord"
    OMEGA0 = "Omega0"
    OMEGABARYON = "OmegaBaryon"
    OMEGALAMBDA = "OmegaLambda"
    REDSHIFT = "Redshift"
    SNAPSHOTFILEBASE = "SnapshotFileBase"
    TIME = "Time"
    UNITLENGTHINCM = "UnitLength_in_cm"
    UNITMASSING = "UnitMass_in_g"
    UNITVELOCITYINCMPERS = "UnitVelocity_in_cm_per_s"


class ArepoRoot(BaseEnum):
    """Root groups in the snapshot/IC files"""
    HEADER = "Header"
    PARAMETERS = "Parameters"
    PARTTYPE = "PartType"
    TREE = "Tree"
    GROUP = "Group"
    SUBHALO = "Subhalo"
    IDS = "IDs"


class ArepoGasFields(BaseEnum):
    """Particle type fields which tend to appear for gas particles (Type 0)"""

    # Fields common to all particle types
    ACCELERATION = "Acceleration"
    COORDINATES = "Coordinates"
    GRAVITYINTERACTIONS = "GravityInteractions"
    MASSES = "Masses"
    PARTICLEIDS = "ParticleIDs"
    POTENTIAL = "Potential"
    SOFTENINGS = "Softenings"
    TIMESTEP = "TimeStep"
    VELOCITIES = "Velocities"

    # Fields common to baryonic particle types
    GFMMETALLICITY = "GFM_Metallicity"
    GFMMETALS = "GFM_Metals"
    METALLICITY = "Metallicity"

    # Fields common to gas only
    ALLOWREFINEMENT = "AllowRefinement"
    BFIELDGRADIENT = "BfieldGradient"
    CELLSPIN = "CellSpin"
    CENTEROFMASS = "CenterOfMass"
    CHEMICALABUNDANCES = "ChemicalAbundances"
    COOLTIME = "CoolTime"
    COSMICRAYSPECIFICENERGY = "CosmicRaySpecificEnergy"
    COSMICRAYSTREAMINGREGULARIZATION = "CosmicRayStreamingRegularization"
    CRPRESSUREGRADIENT = "CRPressureGradient"
    CURLB = "CurlB"
    DENSITY = "Density"
    DENSITYGRADIENT = "DensityGradient"
    DIVBCLEENING = "DivBCleening"
    ELECTRONABUNDANCE = "ElectronAbundance"
    ENERGYDISSIPATION = "EnergyDissipation"
    GFMAGNRADIATION = "GFM_AGNRadiation"
    GFMCOOLINGRATE = "GFM_CoolingRate"
    HIGHRESGASMASS = "HighResGasMass"
    INTERNALENERGY = "InternalEnergy"
    MACHNUMBER = "Machnumber"
    MAGNETICFIELD = "MagneticField"
    MAGNETICFIELDDIVERGENCE = "MagneticFieldDivergence"
    MAGNETICVECTORPOTENTIAL = "MagneticVectorPotential"
    MOLECULARWEIGHT = "MolecularWeight"
    NEUTRALHYDROGENABUNDANCE = "NeutralHydrogenAbundance"
    NUCLEARCOMPOSITION = "NuclearComposition"
    NUCLEARENERGYGENERATIONRATE = "NuclearEnergyGenerationRate"
    PASSIVESCALARS = "PassiveScalars"
    PRESSUREGRADIENT = "PressureGradient"
    PRESSURE = "Pressure"
    RATEOFCHANGEOFMAGNETICFIELD = "RateOfChangeOfMagneticField"
    SFPROBABILITY = "SFProbability"
    SMOOTHEDMAGNETICFIELD = "SmoothedMagneticField"
    SOUNDSPEED = "SoundSpeed"
    STARFORMATIONRATE = "StarFormationRate"
    TEMPERATURE = "Temperature"
    TURBULENTENERGY = "TurbulentEnergy"
    VELOCITYDIVERGENCE = "VelocityDivergence"
    VELOCITYGRADIENT = "VelocityGradient"
    VERTEXVELOCITY = "VertexVelocity"
    VOLUME = "Volume"
    VORTICITY = "Vorticity"


class ArepoStarFields(BaseEnum):
    """Particle type fields which tend to appear for star particles (Type 4)"""
    # Fields common to all particle types
    ACCELERATION = "Acceleration"
    COORDINATES = "Coordinates"
    GRAVITYINTERACTIONS = "GravityInteractions"
    MASSES = "Masses"
    PARTICLEIDS = "ParticleIDs"
    POTENTIAL = "Potential"
    SOFTENINGS = "Softenings"
    TIMESTEP = "TimeStep"
    VELOCITIES = "Velocities"

    # Fields common to baryonic particle types
    GFMMETALLICITY = "GFM_Metallicity"
    GFMMETALS = "GFM_Metals"
    METALLICITY = "Metallicity"

    # Fields common only to stars
    BIRTHPOS = "BirthPos"
    BIRTHVEL = "BirthVel"
    FEEDBACKDONE = "FeedbackDone"
    GFMINITIALMASS = "GFM_InitialMass"
    GFMMASSRELEASED = "GFM_MassReleased"
    GFMMETALRELEASED = "GFM_MetalReleased"
    GFMMETALSRELEASED = "GFM_MetalsReleased"
    GFMSTELLARFORMATIONTIME = "GFM_StellarFormationTime"


class ArepoBHFields(BaseEnum):
    """Particle type fields which tend to appear for Black hole particles (Type 5)"""
    # Fields common to all particle types
    ACCELERATION = "Acceleration"
    COORDINATES = "Coordinates"
    GRAVITYINTERACTIONS = "GravityInteractions"
    MASSES = "Masses"
    PARTICLEIDS = "ParticleIDs"
    POTENTIAL = "Potential"
    SOFTENINGS = "Softenings"
    TIMESTEP = "TimeStep"
    VELOCITIES = "Velocities"

    # Fields common only to black holes
    BHCUMMASSGROWTHQM = "BH_CumMassGrowth_QM"
    BHCUMMASSGROWTHRM = "BH_CumMassGrowth_RM"
    BHDENSITY = "BH_Density"
    BHHSML = "BH_Hsml"
    BHMASS = "BH_Mass"
    BHMASSMETALS = "BH_MassMetals"
    BHMDOT = "BH_Mdot"
    BHMDOTBONDI = "BH_MdotBondi"
    BHMDOTEDDINGTON = "BH_MdotEddington"
    BHMDOTQUASAR = "BH_Mdot_Quasar"
    BHMDOTRADIO = "BH_Mdot_Radio"
    BHPRESSURE = "BH_Pressure"
    BHTIMESTEP = "BH_TimeStep"
    BHU = "BH_U"


class ArepoTreeFields(BaseEnum):
    """Probably unused Tree fields"""
    DESCENDANT = "Descendant"
    FILENR = "FileNr"
    FIRSTHALOINFOFGROUP = "FirstHaloInFOFGroup"
    FIRSTPROGENITOR = "FirstProgenitor"
    GROUPMCRIT200 = "Group_M_Crit200"
    GROUPMCRIT500 = "Group_M_Crit500"
    GROUPMMEAN200 = "Group_M_Mean200"
    NEXTHALOINFOFGROUP = "NextHaloInFOFGroup"
    NEXTPROGENITOR = "NextProgenitor"
    SUBHALOBFLDDISK = "SubhaloBfldDisk"
    SUBHALOBFLDHALO = "SubhaloBfldHalo"
    SUBHALOBHMASS = "SubhaloBHMass"
    SUBHALOBHMDOT = "SubhaloBHMdot"
    SUBHALOCM = "SubhaloCM"
    SUBHALOGASMETALFRACTIONSHALFRAD = "SubhaloGasMetalFractionsHalfRad"
    SUBHALOGASMETALFRACTIONSMAXRAD = "SubhaloGasMetalFractionsMaxRad"
    SUBHALOGASMETALFRACTIONSSFR = "SubhaloGasMetalFractionsSfr"
    SUBHALOGASMETALFRACTIONSSFRWEIGHTED = "SubhaloGasMetalFractionsSfrWeighted"
    SUBHALOGASMETALFRACTIONS = "SubhaloGasMetalFractions"
    SUBHALOGASMETALLICITYHALFRAD = "SubhaloGasMetallicityHalfRad"
    SUBHALOGASMETALLICITYMAXRAD = "SubhaloGasMetallicityMaxRad"
    SUBHALOGASMETALLICITYSFR = "SubhaloGasMetallicitySfr"
    SUBHALOGASMETALLICITYSFRWEIGHTED = "SubhaloGasMetallicitySfrWeighted"
    SUBHALOGASMETALLICITY = "SubhaloGasMetallicity"
    SUBHALOGRNR = "SubhaloGrNr"
    SUBHALOHALFMASSRAD = "SubhaloHalfmassRad"
    SUBHALOHALFMASSRADTYPE = "SubhaloHalfmassRadType"
    SUBHALOIDMOSTBOUND = "SubhaloIDMostbound"
    SUBHALOLEN = "SubhaloLen"
    SUBHALOLENTYPE = "SubhaloLenType"
    SUBHALOMASSINHALFRAD = "SubhaloMassInHalfRad"
    SUBHALOMASSINHALFRADTYPE = "SubhaloMassInHalfRadType"
    SUBHALOMASSINMAXRAD = "SubhaloMassInMaxRad"
    SUBHALOMASSINMAXRADTYPE = "SubhaloMassInMaxRadType"
    SUBHALOMASSINRAD = "SubhaloMassInRad"
    SUBHALOMASSINRADTYPE = "SubhaloMassInRadType"
    SUBHALOMASS = "SubhaloMass"
    SUBHALOMASSTYPE = "SubhaloMassType"
    SUBHALONUMBER = "SubhaloNumber"
    SUBHALOOFFSETTYPE = "SubhaloOffsetType"
    SUBHALOPARENT = "SubhaloParent"
    SUBHALOPOS = "SubhaloPos"
    SUBHALOSFRINHALFRAD = "SubhaloSFRinHalfRad"
    SUBHALOSFRINMAXRAD = "SubhaloSFRinMaxRad"
    SUBHALOSFRINRAD = "SubhaloSFRinRad"
    SUBHALOSFR = "SubhaloSFR"
    SUBHALOSPIN = "SubhaloSpin"
    SUBHALOSTARMETALFRACTIONSHALFRAD = "SubhaloStarMetalFractionsHalfRad"
    SUBHALOSTARMETALFRACTIONSMAXRAD = "SubhaloStarMetalFractionsMaxRad"
    SUBHALOSTARMETALFRACTIONS = "SubhaloStarMetalFractions"
    SUBHALOSTARMETALLICITYHALFRAD = "SubhaloStarMetallicityHalfRad"
    SUBHALOSTARMETALLICITYMAXRAD = "SubhaloStarMetallicityMaxRad"
    SUBHALOSTARMETALLICITY = "SubhaloStarMetallicity"
    SUBHALOSTELLARPHOTOMETRICSMASSINRAD = "SubhaloStellarPhotometricsMassInRad"
    SUBHALOSTELLARPHOTOMETRICSRAD = "SubhaloStellarPhotometricsRad"
    SUBHALOSTELLARPHOTOMETRICS = "SubhaloStellarPhotometrics"
    SUBHALOVELDISP = "SubhaloVelDisp"
    SUBHALOVEL = "SubhaloVel"
    SUBHALOVMAXRAD = "SubhaloVmaxRad"
    SUBHALOVMAX = "SubhaloVmax"
    SUBHALOWINDMASS = "SubhaloWindMass"


class ArepoGroupFields(BaseEnum):
    """Probably unused group fields"""
    GROUPBHMASS = "GroupBHMass"
    GROUPBHMDOT = "GroupBHMdot"
    GROUPCM = "GroupCM"
    GROUPDENSGASANGMOMENTUM = "GroupDensGasAngMomentum"
    GROUPFIRSTSUB = "GroupFirstSub"
    GROUPFUZZOFFSETTYPE = "GroupFuzzOffsetType"
    GROUPGASMETALFRACTIONS = "GroupGasMetalFractions"
    GROUPGASMETALLICITY = "GroupGasMetallicity"
    GROUPLEN = "GroupLen"
    GROUPLENTYPE = "GroupLenType"
    GROUPMASS = "GroupMass"
    GROUPMASSTYPE = "GroupMassType"

    GROUPMTOPHAT200 = "Group_M_TopHat200"
    GROUPNSUBS = "GroupNsubs"
    GROUPPOS = "GroupPos"
    GROUPRCRIT200 = "Group_R_Crit200"
    GROUPRCRIT500 = "Group_R_Crit500"
    GROUPRMEAN200 = "Group_R_Mean200"
    GROUPRTOPHAT200 = "Group_R_TopHat200"
    GROUPSFR = "GroupSFR"
    GROUPSTARMETALFRACTIONS = "GroupStarMetalFractions"
    GROUPSTARMETALLICITY = "GroupStarMetallicity"
    GROUPVEL = "GroupVel"
    GROUPWINDMASS = "GroupWindMass"


class ArepoNames(BaseEnum):

    def __init__(self, list_of_enums):

        for en in list_of_enums:
            for name in en.names():
                setattr(self, name, getattr(en, name))

    def names(self):
        return [i for i in self.__dict__.keys() if i[:1] != "_"]

    def values(self):
        return [self.__dict__[i] for i in self.__dict__.keys() if i[:1] != "_"]


n = ArepoNames([ArepoRoot,
                ArepoHeader, ArepoGasFields, ArepoStarFields, ArepoBHFields,
                ArepoTreeFields, ArepoGroupFields])

d = OrderedDict({ArepoGasFields: {"Root": ArepoRoot.PARTTYPE, "Attr": False},
                 ArepoStarFields: {"Root": ArepoRoot.PARTTYPE, "Attr": False},
                 ArepoHeader: {"Root": ArepoRoot.HEADER, "Attr": True},
                 ArepoGroupFields: {"Root": ArepoRoot.GROUP, "Attr": False},
                 ArepoTreeFields: {"Root": ArepoRoot.TREE, "Attr": False}})


def path(field, ptype=0):
    """Returns the formatted path name of the field given

    :param field: String name of field which is desired
    :type field: str
    :param ptype: Particle type (almost always zero for gas particles)
    :type ptype: int

    :raises ValueError: If field is not found

    :return: Field path and a bool specifying if it's an attribute
    :rtype: (str, bool)
    """

    if field not in n.values():
        raise ValueError(f"Attempting to access non-existant field {field}")
    else:

        for key in d:
            if field in key.values():

                root    = d[key]["Root"]
                isattr  = d[key]["Attr"]

                if isattr:
                    return f"/{root}/", True
                else:
                    return f"/{root}{ptype}/", False

        raise ValueError(f"{field} not found in any known root groups.")
