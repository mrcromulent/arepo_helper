from sim_config import Param, Config, J, Paths
from utilities import arepo_git_versions
import os


class ArepoSimulation(object):

    def __init__(self, proj_name, proj_dir, js, config_eo, param_eo, version="dissipative"):

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

    def check_for_incompatibilities(self, container):

        d = container.data["data"]
        empty_value = container.data["default_value"]

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
        req_name = requirement["name"]
        op = requirement["operator"]
        req_value = requirement["value"]

        if req_name.isupper():
            req_actual_value = self.config.get(req_name)["value"]
        else:
            req_actual_value = self.param.get(req_name)["value"]

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

        self.check_for_incompatibilities(self.param)
        self.check_for_incompatibilities(self.config)

        # This option causes numerical errors when running with AREPO
        if self.version in ["ivo_2016", "dissipative"]:
            assert not self.config.get("MESHRELAX_DENSITY_IN_INPUT")["value"]

        refinement = self.config.get("REFINEMENT_SPLIT_CELLS")["value"] or \
                     self.config.get("REFINEMENT_MERGE_CELLS")["value"]
        ev = self.param.data["default_value"]

        if refinement:
            assert self.param.get("ReferenceGasPartMass")["value"] is not ev
            assert self.param.get("TargetGasMassFactor")["value"] is not ev
            assert self.param.get("RefinementCriterion")["value"] is not ev
            assert self.param.get("DerefinementCriterion")["value"] is not ev
        else:
            assert self.param.get("ReferenceGasPartMass")["value"] is ev
            assert self.param.get("TargetGasMassFactor")["value"] is ev
            assert self.param.get("RefinementCriterion")["value"] is ev
            assert self.param.get("DerefinementCriterion")["value"] is ev

        if not self.config.get("REGULARIZE_MESH_FACE_ANGLE")["value"]:
            assert self.param.get("CellShapingFactor")["value"] is not ev

    def copy_files_to_input(self, file_list):
        if file_list is not None:
            for file in file_list:
                base = os.path.basename(file)
                target = os.path.join(self.paths.input_dir, base)
                os.system(f"cp {file} {target}")

    def write_arepo_compiler(self):

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

        with open(f"{self.paths.arepo_compiler}", 'w') as f:
            f.write('\n'.join(lines))

        os.chmod(self.paths.arepo_compiler, 0o755)

    def write_jobscript(self):

        js = self.js

        lines = [
            "#!/bin/bash",
            " ",
            f"#PBS -N {js[J.NAME]}",
            f"#PBS -P {js[J.PROJECT_CODE]}",
            f"#PBS -q {js[J.QUEUE_TYPE]}",
            f"#PBS -l {J.WALLTIME}={js[J.WALLTIME]}",
            f"#PBS -l {J.MEMORY}={js[J.MEMORY]}",
            f"#PBS -l {J.N_CPUS}={js[J.N_CPUS]}",
            f"#PBS -l {js[J.DIRECTORY]}",
            f"#PBS -lstorage=scratch/{js[J.PROJECT_CODE]}+gdata/{js[J.PROJECT_CODE]}",
            f"#PBS -o ./{self.paths.PBS}",
            f"#PBS -e ./{self.paths.PBS}",
            f"#PBS -m abe -M {js[J.EMAIL]}",
            "\n",
            "module load openmpi",
            "module load hdf5/1.10.5",
            "module load gsl",
            "module load python3-as-python",
            "module load fftw3",
            "\n",
            f"mpirun -np $PBS_NCPUS ./{self.paths.CODE}/Arepo {self.paths.PARAM}"
        ]

        with open(f"{self.paths.jobscript}", 'w') as f:
            f.write('\n'.join(lines))

    def make_systype_file(self):

        systype = ""
        if self.version in ["ivo_2016", "dissipative"]:
            systype = "Raijin"

        elif self.version in ["public_2021"]:
            systype = "Ubuntu"

        with open(self.paths.makefile_systype, "w") as f:
            f.write(f"SYSTYPE=\"{systype}\"")

    def modify_makefile(self):

        if self.version in ["ivo_2016", "dissipative"]:
            with open(self.paths.makefile, "r") as f:
                data = f.read()
                incorrect_string = "OPTIMIZE =  -std=c11 -g -O3 -ipo -m64 -Wall -xCORE-AVX2"
                correct_string = "OPTIMIZE =  -std=c11 -g -O3 -m64 -Wall"
                data = data.replace(incorrect_string, correct_string)

            with open(self.paths.makefile, "w") as f:
                f.write(data)
