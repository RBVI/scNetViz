import falcon
import json
import os
import pathlib
import scanpy as sc

class GetSources:
  """ Return a JSON list of sources we currently have information about """
  manager = None

  def __init__(self, mgr):
    self.manager = mgr

  def on_get(self, req, resp):
    dirs = os.listdir(self.manager.get_fetch_path())
    resp.code = falcon.HTTP_200
    resp.text = json.dumps(dirs)
    return

class GetAccessions:
  """ Return a JSON list of accessions for a given source we currently have information about """
  manager = None

  def __init__(self, mgr):
    self.manager = mgr

  def on_get(self, req, resp, source):
    dirs = os.listdir(self.manager.get_fetch_path())
    if source in dirs:
      accs = os.listdir(self.manager.get_fetch_path()+'/'+source)
      resp.code = falcon.HTTP_200
      resp.text = json.dumps(accs)
      return

    resp.code = falcon.HTTP_400
    resp.test = '{"error":"No files for source '+source+'"}'
    return

class GetInfo:
  """ Return some information about a current cached file """
  manager = None

  def __init__(self, mgr):
    self.manager = mgr

  def on_get(self, req, resp, source, accession):
    # Get the MTX file
    filePath = self.manager.get_fetch_path()+'/'+source+'/'+accession
    if not os.path.isdir(filePath):
      resp.code = falcon.HTTP_400
      resp.test = '{"error":"No data for '+source+'/'+accession+'"}'
      return

    adata = sc.read_10x_mtx(filePath, var_names="gene_symbols", cache=True)
    # Get the number of rows and columns
    genes = adata.n_vars
    cells = adata.n_obs
    # Return it
    resp.code = falcon.HTTP_200
    resp.text = '{"genes":'+str(genes)+',"cells":'+str(cells)+'}'
    return
