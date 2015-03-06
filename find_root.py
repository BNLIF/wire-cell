import os
import os.path as osp
import waflib
import waflib.Utils
import waflib.Logs as msg
from waflib.Configure import conf
from waflib.TaskGen import feature, before_method, after_method, extension, after

_tooldir = osp.dirname(osp.abspath(__file__))


def options(opt):
    opt.load('find_system', tooldir=_tooldir)
    opt.add_option('--with-root', default=None,
                   help="Look for CERN ROOT System at the given path")
    return

def configure(cfg):
    path_list = list()
    for topdir in [getattr(cfg.options, 'with_root', None), os.getenv('ROOTSYS', None)]:
        if topdir:
            path_list.append(osp.join(topdir, "bin"))
        
    kwargs = dict(path_list=path_list)

    cfg.find_program('root-config', var='ROOT-CONFIG', **kwargs)
    cfg.check_cfg(path=cfg.env['ROOT-CONFIG'], uselib_store='ROOTSYS',
                  args = '--cflags --libs --ldflags', package='')
    cfg.find_program('rootcint', var='ROOTCINT', path_list=path_list)
    return

