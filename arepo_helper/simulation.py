from sim_config import Param, Config, J, Paths
import os


class ArepoSimulation(object):
    SYSTYPE = "Raijin"

    def __init__(self, proj_name, proj_dir, js, config_eo, param_eo, version=2019):

        self.project_name   = proj_name
        self.js             = js
        self.version        = version

        # Create a project and set paths
        self.paths = Paths(proj_dir, project_name)

        # Create the Makefiles
        self.make_systype_file(self.SYSTYPE)
        self.write_arepo_compiler()
        self.modify_makefile()

        # Write jobscript, param and config file
        self.write_jobscript()
        self.param = Param(explicit_options=param_eo, version=self.version)
        self.param.write_file(self.paths.param)
        self.config = Config(explicit_options=config_eo, version=self.version)
        self.config.write_file(self.paths.config)

        self.sense_checks()

    def sense_checks(self):

        assert self.param.get("BoxSize")["value"] is not None
        assert self.param.get("InitCondFile")["value"] is not None
        assert self.param.get("OutputDir")["value"] is not None
        assert self.param.get("TimeMax")["value"] is not None
        assert self.param.get("ReferenceGasPartMass")["value"] is not None
        assert self.param.get("MaxVolume")["value"] is not None

        for group in self.param.data["data"]:
            for setting in group["items"]:
                if setting["value"] is not None:
                    for requirement in setting["requires"]:
                        print(f"Checking for existence of requirement: {requirement}")
                        assert self.config.get(requirement)["value"]

                    for incomp in setting["incompatibilities"]:
                        print(f"Checking for incompatible option: {incomp}")
                        assert not (self.config.get(incomp)["value"])

        for group in self.config.data["data"]:
            for setting in group["items"]:
                if setting["value"] is not False:
                    for requirement in setting["requires"]:
                        print(f"Checking for existence of requirement: {requirement}")
                        assert self.config.get(requirement)["value"]

                    for incomp in setting["incompatibilities"]:
                        print(f"Checking for incompatible option: {incomp}")
                        assert not (self.config.get(incomp)["value"])

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
            f"# PBS -N {js[J.NAME]}",
            f"# PBS -P {js[J.PROJECT_CODE]}",
            f"# PBS -q {js[J.QUEUE_TYPE]}",
            f"# PBS -l {J.WALLTIME}={js[J.WALLTIME]}",
            f"# PBS -l {J.MEMORY}={js[J.MEMORY]}",
            f"# PBS -l {J.N_CPUS}={js[J.N_CPUS]}",
            f"# PBS -l {js[J.DIRECTORY]}",
            f"# PBS -lstorage=scratch/{js[J.PROJECT_CODE]}+gdata/{js[J.PROJECT_CODE]}",
            f"# PBS -o {self.paths.pbs_dir}/{js[J.NAME]}.o$PBS_JOBID",
            f"# PBS -e {self.paths.pbs_dir}/{js[J.NAME]}.e$PBS_JOBID",
            f"# PBS -m abe -M {js[J.EMAIL]}",
            " ",
            "module load openmpi",
            "module load hdf5/1.10.5",
            "module load gsl",
            "module load python3-as-python",
            "module load fftw3",
            " ",
            f"mpirun -np $PBS_NCPUS {self.paths.code_dir}/Arepo {self.paths.param}"
        ]

        with open(f"{self.paths.jobscript}", 'w') as f:
            f.write('\n'.join(lines))

    def make_systype_file(self, systype):
        with open(self.paths.makefile_systype, "w") as f:
            f.write(f"SYSTYPE=\"{systype}\"")

    def modify_makefile(self):
        with open(self.paths.makefile, "r") as f:
            data = f.read()
            incorrect_string = "OPTIMIZE =  -std=c11 -g -O3 -ipo -m64 -Wall -xCORE-AVX2"
            correct_string = "OPTIMIZE =  -std=c11 -g -O3 -m64 -Wall"
            data = data.replace(incorrect_string, correct_string)

        with open(self.paths.makefile, "w") as f:
            f.write(data)


if __name__ == '__main__':
    project_name = "dissipative"
    boxsize = 1e10

    js_explicit_options = {
        J.NAME: project_name,
        J.PROJECT_CODE: "y08",
        J.QUEUE_TYPE: "normal",
        J.WALLTIME: "23:59:00",
        J.EMAIL: "uri.burmester@anu.edu.au",
        J.MEMORY: "512gb",
        J.N_CPUS: "240",
        J.DIRECTORY: "wd",
    }

    config_explicit_options = {
        "OUTPUT_PRESSURE": True,
    }
    param_explicit_options = {
        "BoxSize": boxsize,
        "InitCondFile": f"{Paths.INPUT}/bin.dat.ic",
        "OutputDir": f"{Paths.OUTPUT}",
        "TimeMax": 5.0,
        "ReferenceGasPartMass": 2e27,
        "MaxVolume": boxsize ** 3
    }

    simulation = ArepoSimulation(
        project_name,
        "/home/pierre/",
        js_explicit_options,
        config_explicit_options,
        param_explicit_options,
        version=2019)

    simulation.copy_files_to_input([
        "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/species05.txt",
        "/home/pierre/PycharmProjects/arepo_helper/arepo_helper/data/eostable/helm_table.dat"])
