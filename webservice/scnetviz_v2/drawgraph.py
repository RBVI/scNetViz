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

class DrawGraph():
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
            adata = self.handle_drawgraph(req, resp, source, accession)
        except Exception as e:
            logging.error("scNetViz error: "+repr(e))
            resp.status = falcon.HTTP_500
            resp.text = '{"error": "'+str(e)+'"}'
            return

        resp.status = falcon.HTTP_200

        # Note that we add the extra newline for compatability
        resp.text = adata.to_csv(header=None)+'\n'
        return

    def handle_drawgraph(self, req, resp, source, accession):
        # draw_graph arguments
        #  scanpy.tl.draw_graph(adata, layout='fa', init_pos=None, root=None, random_state=0, 
        #                       n_jobs=None, adjacency=None, key_added_ext=None, copy=False, **kwds)
        print('draw_graph')
        n_neighbors = self.manager.get_param_as_int(req, 'n_neighbors', 10)
        layout = self.manager.get_param_as_string(req, 'layout', 'fr')

        adata = self.manager.get_file(source, accession, None)

        adata = self.manager.preprocess(adata, n_neighbors=n_neighbors,
                                        min_genes=self.min_genes,
                                        min_cells=self.min_cells,
                                        normalize=self.normalize,
                                        log1p=self.log1p,
                                        hvg=self.hvg,
                                        scale=self.scale)
        sc.tl.draw_graph(adata, layout=layout)
        return pd.DataFrame(adata.obsm['X_draw_graph_'+layout], index=adata.obs_names)
