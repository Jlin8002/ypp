import json

from steps.astrometry import AstrometryStep

# from steps.calibrate import CalibrateStep
from steps.crop import CropStep
from steps.match import MatchStep
from steps.sextract import SextractStep

STEP_MAP = {
    "astrometry": AstrometryStep,
    "crop": CropStep,
    "match": MatchStep,
    "sextract": SextractStep,
}


class Pipeline(object):
    """A pipeline object used to send data through multiple steps one after another.

    Parameters
    ----------
    steps : list
        A list of step configurations to be run in the pipeline.
    """

    def __init__(self, directory=None, config=None, files=None):
        self.steps = []
        self.files = []

        if directory is not None:
            self.directory = directory
        else:
            self.directory = "."

        if config is not None:
            self.load_config(config)

        if files is not None:
            self.load_files(files)

    def load_config(self, config):
        """Read in a pipeline configuration to the object.

        Parameters
        ----------
        config : Union[str, dict]
            The configuration to load. If config is a string, it must be a valid path.

        Raises
        ------
        TypeError
            If config is not a string or dict.
        NameError
            If a step name does not exist in STEP_MAP.
        """
        if type(config) == str:
            with open(config, "r") as f:
                self.config = json.load(f)
        elif type(config) == dict:
            self.config = config
        else:
            raise TypeError("config must be a string or dict.")

        for step in self.config:
            try:
                self.steps.append(STEP_MAP[step](self.directory, self.config[step]))
            except:
                raise NameError("Step {} not found".format(step))

    def load_files(self, files):
        """Load in a specific set of files to the pipeline.

        Parameters
        ----------
        files : Union[str, list]
            The files to load. If files is a string, it must be a valid path.
        """
        if type(files) == str:
            self.files = []
            with open(files, "r") as f:
                self.files.append(f.readline())
        elif type(files) == list:
            self.files = files
        else:
            raise TypeError("files must be a string or list.")

    def run(self):
        """Run the pipeline.
        """
        for file in self.files:
            result = file
            for step in self.steps:
                result = step.run(result)
