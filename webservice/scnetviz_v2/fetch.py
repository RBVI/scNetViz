import falcon
import json
import os
import pathlib
import scanpy as sc

import scnetviz_v2.gxa

class FetchFile:
  """ Fetch a file from the designated source and cache it """
  manager = None

  def __init__(self, mgr):
    self.manager = mgr

  def on_get(self, req, resp, source, accession):
    filePath = self.manager.get_file_path(source, accession)
    # First, see if we already have it cached
    if (os.path.exists(filePath+"/matrix.mtx")):
      self.manager.get_file(source, accession, resp)
      return;

    # Get the MTX file
    if source == "GXA":
      print("Fetching "+source+"/"+accession)
      scnetviz_v2.gxa.fetch_experiment(accession, filePath)
    else:
      resp.code = falcon.HTTP_400
      resp.text = '{"error":"Only supported for GXA"}'
      return

    self.manager.get_file(source, accession, resp)
    return

