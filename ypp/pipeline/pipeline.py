class Pipeline(object):
    def __init__(self, directory=None):
        self.directory = directory
        self.steps = []
        self.files = []

    def select_files(self):
        pass

    def run(self):
        for files in self.files:
            pass
