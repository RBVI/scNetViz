import site
import sys
from importlib import import_module
from multiprocessing import Manager

if __name__ == '__main__' or __name__.startswith('_mod_wsgi'):
    site.addsitedir('/usr/local/www/webservices/wsgi-scripts')
    print("scnetviz_v2.wsgi initializing", file=sys.stderr)
    manager = Manager()
    print("scnetviz_v2.wsgi importing app", file=sys.stderr)
    app = import_module("scnetviz_v2.app")
    print("scnetviz_v2.wsgi creating app", file=sys.stderr)
    application = app.create_app(manager)
