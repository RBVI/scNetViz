import falcon
import zipfile
from urllib import parse
import pathlib

class SaveFile:
    def __init__(self, manager):
      self.manager = manager

    def on_post(self, req, resp, source, accession):
        print("scNetViz save file")

        if (source is None or accession is None):
            resp.code = falcon.HTTP_400_BAD_REQUEST
            resp.text = '{"error":"both source and accession must be specified"}'
            return

        print ("Unzipping %s"%matrixFile.file)
        matrix = zipfile.ZipFile(matrixFile.file, 'r')
        matrix.extractall(self.manager.get_fetch_path())

        # This will create the hdf5 cache for us
        self.manager.get_file(source, accession, resp)

        return

