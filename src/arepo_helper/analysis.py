from typing import Any, Optional, Dict


class ArepoAnalyser:

    inner_boxsize: float = 0.0
    cutoff_table: dict[str, list[float]] = dict()
    analysis_options: dict[str, Any] = dict()
    select_column: Optional[int] = None
    weight_by_mass: bool = False

    def __init__(self,
                 analysis_options: Optional[Dict[str, Any]] = None):
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
