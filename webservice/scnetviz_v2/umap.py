import falcon
import scanpy as sc
import pandas as pd
import os
import sys
import shutil
import logging
import zipfile
from urllib import parse
import pathlib
#
import logging

class Umap():
    min_genes=100
    min_cells=1
    normalize=True
    log1p=True
    hvg=True

    def __init__(self, manager):
      self.manager = manager

    def on_get(self, req, resp, source, accession):
        logging.info("scNetViz post")
        self.process(req, resp, source, accession)

    def on_post(self, req, resp, source, accession):
        logging.info("scNetViz post")
        self.process(req, resp, source, accession)

    def process(self, req, resp, source, accession):
        if (source is None or accession is None):
            resp.code = falcon.HTTP_400_BAD_REQUEST
            resp.text = '{"error":"both source and accession must be specified"}'
            return

        print("umap: processing")

        self.min_genes = self.manager.get_param_as_int(req, 'min_genes', self.min_genes)
        self.min_cells = self.manager.get_param_as_int(req, 'min_cells', self.min_cells)
        self.normalize = self.manager.get_param_as_bool(req, 'normalize', True)
        self.log1p = self.manager.get_param_as_bool(req, 'log1p', True)
        self.hvg = self.manager.get_param_as_bool(req, 'hvg', True)
        self.scale = self.manager.get_param_as_bool(req, 'scale', True)

        adata = None
        try:
            adata = self.handle_umap(req, resp, source, accession)
        except Exception as e:
            logging.error("scNetViz error: "+repr(e))
            resp.status = falcon.HTTP_500
            resp.text = '{"error": "'+str(e)+'"}'
            return

        resp.status = falcon.HTTP_200

        # Note that we add the extra newline for compatability
        resp.text = adata.to_csv(header=None)+'\n'
        return

    def handle_umap(self, req, resp, source, accession):
        # umap arguments
        # scanpy.tl.umap(adata, min_dist=0.5, spread=1.0, n_components=2, maxiter=None, 
        #                alpha=1.0, gamma=1.0, negative_sample_rate=5, init_pos='spectral', 
        #                random_state=0, a=None, b=None, copy=False)
        #
        print('umap')
        logging.info("umap")
        n_neighbors = self.manager.get_param_as_int(req, 'n_neighbors', 10)
        min_dist = self.manager.get_param_as_float(req, 'min_dist', 0.5)

        adata = self.manager.get_file(source, accession, None)

        logging.info("preprocessing")
        adata = self.manager.preprocess(adata, n_neighbors=n_neighbors,
                                        min_genes=self.min_genes,
                                        min_cells=self.min_cells,
                                        normalize=self.normalize,
                                        log1p=self.log1p,
                                        hvg=self.hvg,
                                        scale=self.scale)
        print('calculating')
        logging.info("calculating")
        sc.tl.umap(adata, min_dist=min_dist)
        print('returning')
        logging.info("returning")
        return pd.DataFrame(adata.obsm['X_umap'], index=adata.obs_names)
