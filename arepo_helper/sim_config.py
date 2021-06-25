from definitions import ROOT_DIR
import textwrap
import pprint
import json
import git
import os


class J:
    NAME = "name"
    PROJECT_CODE = "project_code"
    QUEUE_TYPE = "queue_type"
    WALLTIME = "walltime"
    EMAIL = "email"
    MEMORY = "mem"
    N_CPUS = "ncpus"
    DIRECTORY = "directory"


class Paths:
    arepo_web_address = "https://github.com/boywert/dissipative/"
    default_folder_name = "dissipative"

    CODE = "code"
    CONFIG = "Config.sh"
    OUTPUT = "output"
    INPUT = "input"
    PBS = "pbs"
    PARAM = "param.txt"
    JOBSCRIPT = "jobscript"
    MAKEFILE = "Makefile"
    MAKEFILE_SYSTYPE = "Makefile.systype"

    def __init__(self, parent_dir, proj_name):

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

    def make_dirs(self):

        os.mkdir(self.project_dir)
        os.mkdir(self.input_dir)
        os.mkdir(self.output_dir)
        os.mkdir(self.pbs_dir)

        git.Git(self.project_dir).clone(self.arepo_web_address)
        os.rename(os.path.join(self.project_dir, self.default_folder_name), self.code_dir)


class Param:

    comment_char    = "% "

    def __init__(self,
                 explicit_options=None,
                 version=2019):

        self.explicit_options   = explicit_options
        self.version            = version
        self.data               = None

        self.load()

    def load(self):
        json_file = os.path.join(ROOT_DIR, f"arepo_helper/data/param/{self.version}.json")

        with open(json_file, 'r') as f:
            self.data = json.load(f)

        for group in self.data["data"]:
            for setting in group["items"]:
                curr_option = setting["name"]
                if curr_option in self.explicit_options:
                    setting["value"] = self.explicit_options[curr_option]

    def comment(self, text_string):
        return self.comment_char + text_string + '\n'
    
    def indent_and_wrap(self, text):
        length_limit = 200

        return textwrap.indent(textwrap.fill(text, width=length_limit), self.comment_char) + '\n'

    def write_file(self, filename, ignored):
        delimiter       = self.comment_char + 50 * '-' + ' ' + '\n'
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
                f.write(self.comment(group['heading']))
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
                            text = name + ' ' + str(value) \
                                   + '\t' + self.comment_char + "IGNORED" + double_newline
                        else:
                            text = name + ' ' + str(value) + double_newline
                        f.write(text)

    def get(self, key):
        for group in self.data["data"]:
            for setting in group["items"]:
                curr_key = setting["name"]
                if key == curr_key:
                    return setting

        raise ValueError(f"{key} not found.")

    def __str__(self):

        d = dict()
        for group in self.data["data"]:
            for setting in group["items"]:
                curr_option = setting["name"]
                if setting["value"]:
                    d[curr_option] = setting["value"]

        s = f"CONFIG: \n " \
            f"{self.version=}. \n" \
            f"Non-trivial settings: \n " \
            f"{pprint.pformat(d)}"

        return s


class Config:
    comment_char    = "# "

    def __init__(self,
                 explicit_options=None,
                 version=2019):

        self.version            = version
        self.explicit_options   = explicit_options
        self.data               = None
        self.load()

    def load(self):

        # Load the default Config from file
        json_file = os.path.join(ROOT_DIR, f"arepo_helper/data/config/{self.version}.json")
        with open(json_file, 'r') as f:
            self.data = json.load(f)

        # Overwrite default values with those in explicit_options
        for group in self.data["data"]:
            for setting in group["items"]:
                curr_option = setting["name"]

                if curr_option in self.explicit_options:
                    setting["value"] = self.explicit_options[curr_option]

    def comment(self, text_string):
        return self.comment_char + text_string + '\n'

    def indent_and_wrap(self, text):
        length_limit = 200

        return textwrap.indent(textwrap.fill(text, width=length_limit), self.comment_char) + '\n'

    def write_file(self, filename, ignored):
        delimiter = self.comment_char + 50 * '-' + ' ' + '\n'
        ev = self.data["default_value"]
        double_newline = '\n\n'

        with open(filename, "w") as f:

            f.write('#!/bin/bash')
            f.write(double_newline)

            # Write info from each group
            for group in self.data["data"]:

                f.write(delimiter)
                f.write(self.comment(group['heading']))
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

                        text = name + ': ' + docs
                        text = self.indent_and_wrap(text)
                        f.write(text)

                        value_text = ""

                        if value is not True:
                            value_text += name + '=' + str(value)
                        else:
                            value_text += name

                        if name in ignored:
                            value_text += '\t' + self.comment_char + "IGNORED"

                        value_text += double_newline
                        f.write(value_text)

    def get(self, key):
        for group in self.data["data"]:
            for setting in group["items"]:
                curr_key = setting["name"]
                if key == curr_key:
                    return setting

        raise ValueError(f"{key} not found.")

    def __str__(self):

        d = dict()
        for group in self.data["data"]:
            for setting in group["items"]:
                curr_option = setting["name"]
                if setting["value"]:
                    d[curr_option] = setting["value"]

        s = f"CONFIG: \n " \
            f"{self.version=}. \n" \
            f"Non-trivial settings: \n " \
            f"{pprint.pformat(d)}"

        return s
