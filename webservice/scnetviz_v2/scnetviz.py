import falcon
import logging
import zipfile
import pathlib
import os

import scanpy as sc

import scnetviz_v2.gxa

FETCHPATH = "/tmp/scnetviz"
CACHEPATH = "/tmp/scanpy_cache"

class ScNetVizManager():
  def __init__(self, manager):
    self.manager = manager
    sc._settings.ScanpyConfig.cachedir = pathlib.Path(CACHEPATH)
    sc._settings.ScanpyConfig.logfile = "/var/tmp/scanpy.log"

  def get_fetch_path(self):
    return FETCHPATH

  def get_file_path(self, source, accession):
    return self.get_fetch_path()+'/'+source+'/'+accession

  def get_file(self, source, accession, resp):
    path = self.get_file_path(source, accession)

    # Make sure we have the cached file
    if not os.path.exists(path+"/matrix.mtx"):
      if source == "GXA":
        scnetviz_v2.gxa.fetch_experiment(accession, path)
      else:
        return None

    adata = sc.read_10x_mtx(path, var_names="gene_symbols", cache=True)

    if resp != None:
      genes = adata.n_vars
      cells = adata.n_obs

      resp.code = falcon.HTTP_200
      resp.text = '{"genes":'+str(genes)+',"cells":'+str(cells)+'}'

    return adata

  def get_param_as_string(self, req, param, default):
      if req.has_param(param):
          return str(req.get_param(param))
      else:
          return default

  def get_param_as_float(self, req, param, default):
      if req.has_param(param):
          return req.get_param_as_float(param)
      else:
          return default

  def get_param_as_int(self, req, param, default):
      if req.has_param(param):
          return req.get_param_as_int(param)
      else:
          return default

  def get_param_as_bool(self, req, param, default):
      if req.has_param(param):
          return req.get_param_as_bool(param)
      else:
          return default

  def preprocess(self, adata, n_neighbors=None, min_genes=100, min_cells=1, 
                 normalize=True, log1p=True, 
                 hvg=True, scale=True):

      # Filter options:
      # scanpy.pp.filter_cells(data, min_counts=None, min_genes=None, max_counts=None, 
      #                        max_genes=None, inplace=True, copy=False)
      # scanpy.pp.filter_genes(data, min_counts=None, min_cells=None, max_counts=None, 
      #                        max_cells=None, inplace=True, copy=False)
      #
      sc.pp.filter_cells(adata, min_genes=min_genes)
      sc.pp.filter_genes(adata, min_cells=min_cells)

      # normalization options
      #  scanpy.pp.normalize_total(adata, target_sum=None, exclude_highly_expressed=False, 
      #                            max_fraction=0.05, key_added=None, layers=None, layer_norm=None, 
      #                            inplace=True)
      if normalize is True:
        sc.pp.normalize_total(adata)

      if log1p is True:
          sc.pp.log1p(adata)

      #
      # scanpy.pp.highly_variable_genes(adata, min_disp=None, max_disp=None, min_mean=None, 
      #                                 max_mean=None, n_top_genes=None, n_bins=20, flavor='seurat', 
      #                                 subset=False, inplace=True, batch_key=None)
      if hvg is True:
          sc.pp.highly_variable_genes(adata)
          adata = adata[:, adata.var['highly_variable']]

      # Scal options
      # scanpy.pp.scale(data, zero_center=True, max_value=None, copy=False)
      if scale is True:
          sc.pp.scale(adata, max_value=10)

      if (n_neighbors != None):
          # neighbors options:
          #  scanpy.pp.neighbors(adata, n_neighbors=15, n_pcs=None, use_rep=None, knn=True, 
          #                      random_state=0, method='umap', metric='euclidean', metric_kwds={}, copy=False)
          try:
              sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=40)
          except Exception as e:
              logging.debug('preprocess: got exception calculating neighbors: '+str(e))
      else:
          # PCA Options
          #  scanpy.tl.pca(data, n_comps=50, zero_center=True, svd_solver='auto', random_state=0, 
          #                return_info=False, use_highly_variable=None, dtype='float32', copy=False, 
          #                chunked=False, chunk_size=None)
          sc.tl.pca(adata, svd_solver='arpack')

      return adata
