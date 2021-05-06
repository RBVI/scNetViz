import falcon
from falcon_multipart.middleware import MultipartMiddleware
from multiprocessing import Manager
from scnetviz_v2.scnetviz import ScNetVizManager
from scnetviz_v2.fetch import FetchFile
from scnetviz_v2.load_file import LoadFile
from scnetviz_v2.save_file import SaveFile
from scnetviz_v2.get_info import GetInfo,GetSources,GetAccessions
from scnetviz_v2.umap import Umap
from scnetviz_v2.tsne import Tsne
from scnetviz_v2.drawgraph import DrawGraph
from scnetviz_v2.louvain import Louvain
from scnetviz_v2.leiden import Leiden

import logging
import sys

manager = None

class Test:
    def on_get(self, req, resp):
        resp.status = falcon.HTTP_200

        # Note that we add the extra newline for compatability
        resp.body = '{"result": "Well, foo to you, too"}\n'
        return

def create_app(mgr: Manager):
    manager = mgr
    app = falcon.App(middleware=[MultipartMiddleware()])
    logging.basicConfig(filename="/var/tmp/scnetviz_api.log",level=logging.INFO)
    print("scnetviz service initializing", file=sys.stderr)

    scNetViz = ScNetVizManager(mgr)

    # Fetch a file and cache it
    app.add_route('/fetch/{source}/{accession}', FetchFile(scNetViz))

    # Return a cached file
    app.add_route('/load/{source}/{accession}', LoadFile(scNetViz))

    # Upload a file as a compressed mtx and cache it
    app.add_route('/save/{source}/{accession}', SaveFile(scNetViz))

    # Return a list of cached sources
    app.add_route('/info', GetSources(scNetViz))

    # Return a list of cached files for a given source
    app.add_route('/info/{source}', GetAccessions(scNetViz))

    # Return the number of rows and columns in a cached file
    app.add_route('/info/{source}/{accession}', GetInfo(scNetViz))

    app.add_route('/umap/{source}/{accession}', Umap(scNetViz))
    app.add_route('/tsne/{source}/{accession}', Tsne(scNetViz))
    app.add_route('/drawgraph/{source}/{accession}', DrawGraph(scNetViz))
    app.add_route('/louvain/{source}/{accession}', Louvain(scNetViz))
    app.add_route('/leiden/{source}/{accession}', Leiden(scNetViz))
    #app.add_route('/foo', Test())

    print("scnetviz_v2 initializion done", file=sys.stderr)
    return app

def get_app():
    return create_app(Manager())
