from definitions import DATA_DIR, AREPO_SRC_DIR
from utilities import arepo_git_versions
from typing import Any
import textwrap
import shutil
import pprint
import json
import git
import os


class J:
    """Jobscript fields"""
    NAME = "name"
    PROJECT_CODE = "project_code"
    QUEUE_TYPE = "queue_type"
    WALLTIME = "walltime"
    EMAIL = "email"
    MEMORY = "mem"
    N_CPUS = "ncpus"
    DIRECTORY = "directory"


class Paths:
    """Directories and files which are created an modified in an AREPO simulation."""

    CODE = "code"
    CONFIG = "Config.sh"
    OUTPUT = "output"
    INPUT = "input"
    PBS = "pbs"
    PARAM = "param.txt"
    JOBSCRIPT = "jobscript"
    MAKEFILE = "Makefile"
    MAKEFILE_SYSTYPE = "Makefile.systype"

    def __init__(self,
                 parent_dir: str,
                 proj_name: str,
                 version: str = "dissipative"):
        """Constructor which creates the folders and files required for an AREPO simulation.

        :param parent_dir: Directory where the simulation directory should be created.
        :type parent_dir: str
        :param proj_name: Name of the simulation directory (also used in some files)
        :type proj_name: str
        :param version: Version string specifying the version of the code to be used
        :type version: str

        :raises OSError: Raised if the directory already exists

        :return: Paths object
        :rtype: Paths
        """

        self.version = version
        self.project_dir = os.path.join(parent_dir, proj_name)
        self.input_dir = os.path.join(self.project_dir, Paths.INPUT)
        self.output_dir = os.path.join(self.project_dir, Paths.OUTPUT)
        self.pbs_dir = os.path.join(self.project_dir, Paths.PBS)
        self.code_dir = os.path.join(self.project_dir, Paths.CODE)
        self.jobscript = os.path.join(self.project_dir, self.JOBSCRIPT)
        self.param = os.path.join(self.project_dir, self.PARAM)
        self.makefile = os.path.join(self.code_dir, self.MAKEFILE)
        self.config = os.path.join(self.code_dir, self.CONFIG)
        self.makefile_systype = os.path.join(self.code_dir, self.MAKEFILE_SYSTYPE)
        self.arepo_compiler = os.path.join(self.code_dir, "make_arepo.sh")

        if os.path.exists(self.project_dir):
            raise OSError
        else:
            self.make_dirs()

    def make_dirs(self) -> None:
        """Creates the relevant directories."""

        os.mkdir(self.project_dir)
        os.mkdir(self.input_dir)
        os.mkdir(self.output_dir)
        os.mkdir(self.pbs_dir)

        self.get_code()

    def get_code(self) -> None:
        """Clones code from Github or copies it from the local folder.

        :raises ValueError: If unknown version requested

        :rtype: None
        """

        if self.version not in arepo_git_versions:
            raise ValueError(f"Unknown version requested: {self.version}")

        else:
            wa  = arepo_git_versions[self.version]["url"]
            c_id = arepo_git_versions[self.version]["commit_id"]

            if wa is not None:
                # Get source code from web
                repo = git.Repo.clone_from(wa, self.code_dir, no_checkout=True)
                repo.git.checkout(c_id)
            else:
                # Get source code from local source
                local_source = os.path.join(AREPO_SRC_DIR, c_id)
                shutil.copytree(local_source, self.code_dir)


class ArepoInput:
    """Base class for Param and Config classes"""
    comment_char = None
    length_limit = 200
    input_name   = None

    def __init__(self,
                 explicit_options: dict,
                 version: str) -> None:
        """Constructor for base class

        :param explicit_options: Dict of options to override defaults
        :type explicit_options: dict
        :param version: String version name to specify AREPO version
        :type version: str

        :return: ArepoInput
        :rtype: ArepoInput
        """
        self.explicit_options = explicit_options
        self.version = version
        self.data = None

    def comment(self,
                text_string: str) -> str:
        """Returns the text in text_string formatted as a comment.

        :param text_string: Text to be commented
        :type text_string: str

        :return: Commented text
        :rtype: str
        """
        return self.comment_char + text_string + "\n"

    def indent_and_wrap(self,
                        text: str) -> str:
        """Indents and wraps long line of text.

        :param text: Text to be wrapped
        :type text: str

        :return: Wrapped text
        :rtype: str
        """
        return textwrap.indent(textwrap.fill(text, width=self.length_limit), self.comment_char) + "\n"

    def get(self,
            key: str) -> Any:
        """Gets data from a particular element of the Config/Param

        :param key: Name of element
        :type key: str

        :raises ValueError: If element is not found

        :return: Value associated with element
        :rtype: Any
        """

        for group in self.data["data"]:
            for setting in group["items"]:
                curr_key = setting["name"]
                if key == curr_key:
                    return setting["value"]

        raise ValueError(f"{key} not found.")

    def load(self) -> None:
        """Loads default values from the JSON file associated with self.version."""

        # Load the default Config from file
        c_id = arepo_git_versions[self.version]["commit_id"]
        json_file = os.path.join(DATA_DIR, f"{self.input_name}/{c_id}.json")
        with open(json_file, "r") as f:
            self.data = json.load(f)

        # Overwrite default values with those in explicit_options
        for group in self.data["data"]:
            for setting in group["items"]:
                curr_option = setting["name"]
                if curr_option in self.explicit_options:
                    setting["value"] = self.explicit_options[curr_option]

    def __str__(self) -> str:

        d = dict()
        for group in self.data["data"]:
            for setting in group["items"]:
                curr_option = setting["name"]
                if setting["value"]:
                    d[curr_option] = setting["value"]

        s = f"{self.input_name}: \n " \
            f"{self.version=}. \n" \
            f"Non-trivial settings: \n " \
            f"{pprint.pformat(d)}"

        return s


class Param(ArepoInput):
    """Class to hold information associated with param.txt files."""
    comment_char = "% "
    input_name = "param"

    def __init__(self,
                 explicit_options: dict = None,
                 version: str = "dissipative") -> None:
        """Constructor for param.txt files

        :param explicit_options: Explicit param.txt options
        :type explicit_options: dict
        :param version: Version of AREPO to be used
        :type version: str

        :return: Param object
        :rtype: Param
        """

        super(Param, self).__init__(explicit_options, version)
        self.load()

    def write_file(self,
                   filename: str,
                   ignored: list[str]) -> None:
        """Writes the param.txt file to filename, specifying which values have been ignored.

        :param filename: Location to which file is written
        :type filename: str
        :param ignored: Ignored options
        :type ignored: list[str]
        """
        delimiter       = self.comment_char + 50 * "-" + " " + "\n"
        ev              = self.data["default_value"]
        double_newline  = "\n\n"

        with open(filename, "w") as f:

            heading = self.indent_and_wrap(self.data["heading"])

            f.write(delimiter)
            f.write(heading)
            f.write(delimiter)
            f.write(double_newline)

            # Write info from each group
            for group in self.data["data"]:

                f.write(delimiter)
                f.write(self.comment(group["heading"]))
                f.write(delimiter)
                f.write(double_newline)

                for s in group["items"]:
                    name = s["name"]
                    docs = s["docs"]
                    value = s["value"]

                    # Write
                    if value is not ev:
                        f.write(self.comment(name))
                        if docs is not None:
                            docs = self.indent_and_wrap(docs)
                            f.write(docs)

                        if name in ignored:
                            text = name + " " + str(value) \
                                   + "\t" + self.comment_char + "IGNORED" + double_newline
                        else:
                            text = name + " " + str(value) + double_newline
                        f.write(text)


class Config(ArepoInput):
    """Class to hold information associated with Config.sh files."""
    comment_char = "# "
    input_name = "config"

    def __init__(self,
                 explicit_options: dict = None,
                 version: str = "dissipative") -> None:
        """Constructor for Config.sh files

        :param explicit_options: Explicit Config.sh options
        :type explicit_options: dict
        :param version: Version of AREPO to be used
        :type version: str

        :return: Config object
        :rtype: Config
        """

        super(Config, self).__init__(explicit_options, version)
        self.load()

    def write_file(self,
                   filename: str,
                   ignored: list[str]) -> None:
        """Writes the Config.sh file to filename, specifying which values have been ignored.

        :param filename: Location to which file is written
        :type filename: str
        :param ignored: Ignored options
        :type ignored: list[str]
        """
        delimiter = self.comment_char + 50 * "-" + " " + "\n"
        ev = self.data["default_value"]
        double_newline = "\n\n"

        with open(filename, "w") as f:

            f.write("#!/bin/bash")
            f.write(double_newline)

            # Write info from each group
            for group in self.data["data"]:

                f.write(delimiter)
                f.write(self.comment(group["heading"]))
                f.write(delimiter)
                f.write(double_newline)

                for s in group["items"]:
                    name = s["name"]
                    docs = s["docs"]
                    value = s["value"]

                    # Write
                    if value is not ev:
                        if docs is None:
                            docs = ""

                        text = name + ": " + docs
                        text = self.indent_and_wrap(text)
                        f.write(text)

                        value_text = ""

                        if value is not True:
                            value_text += name + "=" + str(value)
                        else:
                            value_text += name

                        if name in ignored:
                            value_text += "\t" + self.comment_char + "IGNORED"

                        value_text += double_newline
                        f.write(value_text)
