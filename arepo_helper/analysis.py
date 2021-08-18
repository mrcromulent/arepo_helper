class ArepoAnalyser(object):

    inner_boxsize       = None
    cutoff_table        = dict()
    analysis_options    = dict()
    select_column       = None

    def __init__(self, analysis_options=None):
        """Arepo Analyser object

        :param analysis_options: Explicit options
        :type analysis_options: dict
        :return: Arepo Analyser object
        :rtype: ArepoAnalyser
        """

        if analysis_options is not None:
            self.analysis_options = analysis_options
            for key in analysis_options.keys():
                setattr(self, key, analysis_options[key])
        else:
            self.analysis_options = dict()
