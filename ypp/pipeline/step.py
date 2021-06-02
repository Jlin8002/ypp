class PipeStep(object):
    def __init__(self, directory=None):
        self.directory = directory
        self.input_file_ext = None
        self.input_step_ext = None
        self.output_file_ext = None
        self.output_step_ext = None
        self.input_data = None
        self.output_data = None

    def set_input_exts(self, file_ext, step_ext):
        self.input_file_ext = file_ext
        self.input_step_ext = step_ext

    def set_output_exts(self, file_ext, step_ext):
        self.output_file_ext = file_ext
        self.output_step_ext = step_ext

    def set_directory(self, directory):
        self.directory = directory

    def read_data(self):
        if self.input_file_ext == "fits" or self.input_file_ext == "fit":
            pass
        elif self.input_file_ext == "tif" or self.input_file_ext == "tiff":
            pass
        elif self.input_file_ext == "csv":
            pass
        else:
            raise "Unrecognized file extension."

    def write_data(self):
        if self.output_file_ext == "fits" or self.output_file_ext == "fit":
            pass
        elif self.output_file_ext == "tif" or self.output_file_ext == "tiff":
            pass
        elif self.output_file_ext == "csv":
            pass
        else:
            raise "Unrecognized file extension."

    def process_data(self):
        pass

    def run(self):
        self.read_data()
        self.process_data()
        self.write_data()

