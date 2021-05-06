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

class Tsne():
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
            adata = self.handle_tsne(req, resp, source, accession)
        except Exception as e:
            logging.error("scNetViz error: "+repr(e))
            resp.status = falcon.HTTP_500
            resp.text = '{"error": "'+str(e)+'"}'
            return

        resp.status = falcon.HTTP_200

        # Note that we add the extra newline for compatability
        resp.text = adata.to_csv(header=None)+'\n'
        return

    def handle_tsne(self, req, resp, source, accession):
        # tSNE arguments
        #  scanpy.tl.tsne(adata, n_pcs=None, use_rep=None, perplexity=30, early_exaggeration=12, 
        #                 learning_rate=1000, random_state=0, use_fast_tsne=True, n_jobs=None, copy=False)
        #
        print('tsne')
        adata = self.manager.get_file(source, accession, None)
        adata = self.manager.preprocess(adata,
                                        min_genes=self.min_genes,
                                        min_cells=self.min_cells,
                                        normalize=self.normalize,
                                        log1p=self.log1p,
                                        hvg=self.hvg,
                                        scale=self.scale)
        n_pcs = self.manager.get_param_as_int(req, 'n_pcs', None)
        perplexity = self.manager.get_param_as_float(req, 'perplexity', 30.0)
        learning_rate = self.manager.get_param_as_float(req, 'learning_rate', 1000.0)
        early_exaggeration = self.manager.get_param_as_float(req, 'early_exaggeration', 12.0)

        sc.tl.tsne(adata, n_pcs=n_pcs, perplexity=perplexity, learning_rate=learning_rate, 
                   early_exaggeration=early_exaggeration)
        return pd.DataFrame(adata.obsm['X_tsne'], index=adata.obs_names)
