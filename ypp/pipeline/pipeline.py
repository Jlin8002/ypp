class Pipeline(object):
    def __init__(self, directory=None):
        self.directory = directory
        self.steps = []
