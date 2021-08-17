"""
test
"""


class ArepoAnalyser(object):

    inner_boxsize       = None
    cutoff_table        = dict()
    analysis_options    = dict()

    def __init__(self, analysis_options=None):

        if analysis_options is not None:
            self.analysis_options = analysis_options
            for key in analysis_options.keys():
                setattr(self, key, analysis_options[key])
        else:
            self.analysis_options = dict()
