"""
test
"""


class ArepoAnalyser(object):

    inner_boxsize       = None
    cutoff_table        = dict()
    highlights          = dict()
    slice_width         = 0.005
    line_radius         = 0.005

    def __init__(self, analysis_options=None):

        if analysis_options is not None:
            for key in analysis_options.keys():
                setattr(self, key, analysis_options[key])


if __name__ == "__main__":
    pass
