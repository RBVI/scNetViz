import falcon
import zipfile
from urllib import parse
import pathlib
import mimetypes
import os

class LoadFile:
    def __init__(self, manager):
      self.manager = manager

    def on_get(self, req, resp, source, accession):
      return self.process_file(resp, source, accession)

    def on_post(self, req, resp, source, accession):
      return self.process_file(resp, source, accession)

    def process_file(self, resp, source, accession):
        print("scNetViz load file")

        if (source is None or accession is None):
            resp.code = falcon.HTTP_400_BAD_REQUEST
            resp.text = '{"error":"both source and accession must be specified"}'
            return

        print ("Zipping %s/%s"%(source,accession))
        zip_file = "/tmp/%s_%s.zip"%(source,accession)
        with zipfile.ZipFile(zip_file, 'w') as myzip:
          cache_path = self.manager.get_file_path(source,accession)
          myzip.write(cache_path+"/matrix.mtx")
          myzip.write(cache_path+"/barcodes.tsv")
          myzip.write(cache_path+"/genes.tsv")

        resp.code = falcon.HTTP_200
        resp.content_type = mimetypes.guess_type(zip_file)[0]
        resp.stream = open(zip_file, 'rb')
        resp.content_length = os.path.getsize(zip_file)
        return

