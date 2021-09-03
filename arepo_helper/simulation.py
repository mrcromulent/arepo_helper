from utilities import arepo_git_versions
from sim_config import Param, Config, J, Paths
from names import n, ArepoHeader as h
from pyhelm_eos import loadhelm_eos
from h5_file import ArepoICs
import numpy as np
import os


class ArepoSimulation(object):
    """Top level class for creating and managing AREPO simulations. """

    def __init__(self, proj_name, proj_dir, js, config_eo, param_eo, version="dissipative"):
        """Creates an Arepo Simulation object using Config and Param objects and additional version information.

        :param proj_name: Project name
        :type proj_name: str
        :param proj_dir: Directory to which simulation should be written
        :type proj_dir: str
        :param js: Jobscript dict for NCI
        :type js: dict
        :param config_eo: Explicit options which will override defaults in the Config object
        :type config_eo: dict
        :param param_eo: Explicit options which will override defaults in the Param object
        :type param_eo: dict
        :param version: Arepo version string
        :type version: str

        :raises ValueError: If unknown version is requested

        :return: ArepoSimulation object
        :rtype: ArepoSimulation
        """

        self.project_name = proj_name
        self.js = js
        self.ignored = []

        if version not in arepo_git_versions:
            raise ValueError(f"Unknown version requested: {version}")
        else:
            self.version = version

        #
        self.param = Param(explicit_options=param_eo, version=self.version)
        self.config = Config(explicit_options=config_eo, version=self.version)
        self.sense_checks()

        # Create a project and set paths
        self.paths = Paths(proj_dir, proj_name, version=version)

        # Create the Makefiles
        self.make_systype_file()
        self.write_arepo_compiler()
        self.modify_makefile()

        #
        self.write_jobscript()
        self.param.write_file(self.paths.param, self.ignored)
        self.config.write_file(self.paths.config, self.ignored)

    def check_for_incompatibilities(self, arepo_input):
        """Checks for incompatible options

        :param arepo_input: Config or Param object
        :type arepo_input: Config|Param
        :rtype: None
        """

        d = arepo_input.data["data"]
        empty_value = arepo_input.data["default_value"]

        for group in d:
            for s in group["items"]:

                name = s["name"]
                required = s["required"]
                value = s["value"]
                ignored_if = s["ignored_if"]
                requirements = s["requirements"]

                # Required settings must have the non-empty value
                if required and value == empty_value:
                    raise ValueError(f"{name} is required and was not set.")

                if value is not empty_value:
                    for requirement in requirements:
                        met = self.check_requirement_met(requirement)
                        assert met, f"{name}: Requirement not met: {requirement}"

                    # Check whether any ignored conditions are true
                    ignored = False
                    for ic in ignored_if:
                        if self.check_requirement_met(ic):
                            print(f"{name}: Ignored due to condition: {ic}")
                            ignored = True

                    if ignored:
                        self.ignored.append(name)

    def check_requirement_met(self, requirement):
        """Checks if a particular requirement is met as specified in the Config/Param JSON file.

        :param requirement: Dict specifying requirement
        :type requirement: dict

        :raises ValueError: If requirement operator is known

        :return: Boolean specifying if requirement is met
        :rtype: bool
        """

        req_name = requirement["name"]
        op = requirement["operator"]
        req_value = requirement["value"]

        if req_name.isupper():
            req_actual_value = self.config.get(req_name)
        else:
            req_actual_value = self.param.get(req_name)

        if op == "EQ":
            return req_actual_value == req_value
        elif op == "NEQ":
            return req_actual_value != req_value
        elif op == "MAX":
            return req_actual_value <= req_value
        elif op == "MIN":
            return req_actual_value >= req_value
        else:
            raise ValueError(f"Unknown operator: {op}")

    def sense_checks(self):
        """Performs sense checks to see if you have anything incompatible in your simulation specifications."""

        self.check_for_incompatibilities(self.param)
        self.check_for_incompatibilities(self.config)

        refinement = self.config.get("REFINEMENT_SPLIT_CELLS") or self.config.get("REFINEMENT_MERGE_CELLS")
        ev = self.param.data["default_value"]

        if refinement:
            assert self.param.get("ReferenceGasPartMass") is not ev
            assert self.param.get("TargetGasMassFactor") is not ev
            assert self.param.get("RefinementCriterion") is not ev
            assert self.param.get("DerefinementCriterion") is not ev
        else:
            assert self.param.get("ReferenceGasPartMass") is ev
            assert self.param.get("TargetGasMassFactor") is ev
            assert self.param.get("RefinementCriterion") is ev
            assert self.param.get("DerefinementCriterion") is ev

        if not self.config.get("REGULARIZE_MESH_FACE_ANGLE"):
            assert self.param.get("CellShapingFactor") is not ev

        # Relax Runtime requires a long time scale
        if self.config.get("RELAX_RUNTIME"):
            assert self.param.get("TimeMax") >= 40.0

    def copy_files_to_input(self, file_list):
        """Copies files to the input folder for your simulation (e.g. species.txt files)

        :param file_list: List of files to copy
        :type file_list: list
        """
        if file_list is not None:
            for file in file_list:
                base = os.path.basename(file)
                target = os.path.join(self.paths.input_dir, base)
                os.system(f"cp {file} {target}")

    def write_arepo_compiler(self):
        """Writes the AREPO compiler file."""

        lines = [
            "# !/bin/bash",
            " ",
            "module load python3",
            "module load python3-as-python",
            "module load gsl",
            "module load gcc",
            "module load fftw3",
            "module load hdf5/1.10.5",
            "module load openmpi",
            " ",
            "make"
        ]

        with open(f"{self.paths.arepo_compiler}", "w") as f:
            f.write("\n".join(lines))

        # Make the compiler executable
        os.chmod(self.paths.arepo_compiler, 0o755)

    def write_jobscript(self):
        """Writes the jobscript file. """

        lines = [
            "#!/bin/bash",
            " ",
            f"#PBS -N {self.js[J.NAME]}",
            f"#PBS -P {self.js[J.PROJECT_CODE]}",
            f"#PBS -q {self.js[J.QUEUE_TYPE]}",
            f"#PBS -l {J.WALLTIME}={self.js[J.WALLTIME]}",
            f"#PBS -l {J.MEMORY}={self.js[J.MEMORY]}",
            f"#PBS -l {J.N_CPUS}={self.js[J.N_CPUS]}",
            f"#PBS -l {self.js[J.DIRECTORY]}",
            f"#PBS -lstorage=scratch/{self.js[J.PROJECT_CODE]}+gdata/{self.js[J.PROJECT_CODE]}",
            f"#PBS -o ./{self.paths.PBS}",
            f"#PBS -e ./{self.paths.PBS}",
            f"#PBS -m abe -M {self.js[J.EMAIL]}",
            "\n",
            "module load openmpi",
            "module load hdf5/1.10.5",
            "module load gsl",
            "module load python3-as-python",
            "module load fftw3",
            "\n",
            f"mpirun -np $PBS_NCPUS ./{self.paths.CODE}/Arepo {self.paths.PARAM}"
        ]

        with open(f"{self.paths.jobscript}", "w") as f:
            f.write("\n".join(lines))

    def make_systype_file(self):
        """Writes the Makefile.systype file."""

        systype = ""
        if self.version in ["ivo_2016", "dissipative"]:
            systype = "Raijin"

        elif self.version in ["public_2021"]:
            systype = "Ubuntu"

        with open(self.paths.makefile_systype, "w") as f:
            f.write(f"SYSTYPE=\"{systype}\"")

    def modify_makefile(self):
        """Removes incorrect specifications (ipo, -xCORE-AVX2) from the Raijin specifications in the Makefile."""

        if self.version in ["ivo_2016", "dissipative"]:
            with open(self.paths.makefile, "r") as f:
                data = f.read()
                incorrect_string = "OPTIMIZE =  -std=c11 -g -O3 -ipo -m64 -Wall -xCORE-AVX2"
                correct_string = "OPTIMIZE =  -std=c11 -g -O3 -m64 -Wall"
                data = data.replace(incorrect_string, correct_string)

            with open(self.paths.makefile, "w") as f:
                f.write(data)

    def validate_input_file(self, arepo_ic_file, helm_file, species_file):
        """Performs a series of checks to determine whether your simulation will be immediately rejected.

        .. notes::
            - Checks whether all particles are within the boxsize,
            - Checks the calculated temperature using EOS_DEGENERATE

        :param arepo_ic_file: Location of Arepo Initial Conditions file
        :type arepo_ic_file: str
        :param helm_file: Helmholtz EOS file name
        :type helm_file: str
        :param species_file: Species.txt file
        :type species_file: str
        """

        ic_h5 = ArepoICs(arepo_ic_file)

        boxsize     = self.param.get(h.BOXSIZE)
        n_particles = ic_h5.num_particles()

        #
        coords      = ic_h5.get_from_h5(n.COORDINATES)
        density     = ic_h5.get_from_h5(n.DENSITY)
        masses      = ic_h5.get_from_h5(n.MASSES)
        xnuc        = ic_h5.get_from_h5(n.NUCLEARCOMPOSITION)
        intern_ener = ic_h5.get_from_h5(n.INTERNALENERGY)

        #
        eos = loadhelm_eos(helm_file, species_file, True)

        # This option implies that density figures are stored in the masses field
        if self.config.get("MESHRELAX_DENSITY_IN_INPUT"):
            assert np.all(masses) < 1e8

        #
        for i in range(n_particles):

            # Assert that all positions are within [0, boxsize]
            for j in range(0, 3):
                assert coords[i, j] >= 0.0
                assert coords[i, j] <= boxsize

            #
            if self.config.get("EOS_DEGENERATE"):

                min_temp = 1e03
                max_temp = 1e10

                temp_calculated, p_calculated = eos.egiven(density[i], xnuc[i], intern_ener[i])

                if temp_calculated <= 0.0:
                    raise ValueError(f"{temp_calculated=}")
                else:
                    if temp_calculated > max_temp:
                        temp_calculated = max_temp
                    elif temp_calculated < min_temp:
                        temp_calculated = min_temp

                    e_calculated, dedt, p_calculated, csnd = eos.tgiven(density[i], xnuc[i], temp_calculated)

                if temp_calculated < min_temp or temp_calculated > max_temp:
                    egy_injection = e_calculated - intern_ener[i]
                    print(f"{e_calculated=}, {egy_injection=}")

        # Volume
        volume = np.divide(masses, density)
        print(f"{max(volume)=}. {min(volume)=}")
