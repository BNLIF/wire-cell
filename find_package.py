# -*- python -*-
'''This tool implements a source package following a few contentions.

Your source package may build any combination of the following:

 - shared libraries 
 - headers exposing an API to libraries
 - a ROOT dictionary for this API
 - main programs
 - test programs

This tool will produce various methods on the build context.  You can
avoid passing <name> to them if you set APPNAME in your wscript file.

'''

import os.path as osp
from waflib.Utils import to_list
from waflib.Configure import conf
import waflib.Context

_tooldir = osp.dirname(osp.abspath(__file__))

def options(opt):
    opt.load('compiler_cxx')
    opt.load('waf_unit_test')
    opt.load('find_root', tooldir=_tooldir)
def configure(cfg):
    cfg.load('compiler_cxx')
    cfg.load('waf_unit_test')
    cfg.load('find_root', tooldir=_tooldir)

    pass




def _includes(bld, includes):
    if not includes:
        includes = ['inc']
    return to_list(includes)

def _name(name=None):
    return name or getattr(waflib.Context.g_module, 'APPNAME')

def _headers(bld, headers, name=None):
    name = _name(name)
    print '_headers NAME',name
    if not headers:
        incdir = bld.path.find_dir('inc/%s' % name)
        return incdir.ant_glob('*.h')
    if type(headers) == type(""):
        return bld.path.ant_glob(headers)
    return headers


    


@conf
def shared_library(bld, name=None, srcdir='src', source = '', includes ='', use=''):
    '''
    Make a shared library named <name>

    Place source in src/ or set <srcdir> or explicitly give <source>.

    Place public API headers in inc/ or explicitly give <includes>
    '''
    if not source:
        srcdir = bld.path.find_dir(srcdir)
        source = srcdir.ant_glob('*.cxx')
    source = to_list(source)
    use = to_list(use) + ['ROOTSYS']

    bld.shlib(source = source, 
              target = _name(name), 
              includes = _includes(bld, includes), 
              use = use)


@conf
def api_headers(bld, name=None, headers = ''):
    '''
    Install headers for API <name>.

    Place headers in inc/<name>/*.h or explicitly set headers.
    '''
    name = _name(name)
    bld.install_files('${PREFIX}/include/%s' % name, 
                      _headers(bld, headers, name))


@conf
def root_dictionary(bld, name=None, linkdef=None,
                    headers = None, includes=None):
    '''
    Build a root dictionary library for package <name>.
    
    Make a dict/LinkDef.h or explicitly set <linkdef>.

    Use all API headers in inc/<name>/*.h or explicitly set headers
    '''
    name = _name(name)
    linkdef = linkdef or "dict/LinkDef.h"

    bld.gen_rootcling_dict(name, linkdef,
                           headers = _headers(bld, headers, name),
                           includes = _includes(bld, includes))

### fixme: make this configurable via options
    # generate rootcint dictionary and build a shared library
#    bld.gen_rootcint_dict(name, linkdef,
#                          headers = _headers(bld, headers, name),
#                          includes = _includes(bld, includes))



@conf
def test_program(bld, source = '', includes=''):
    '''
    Build and run a test program named <name>.

    Using source tests/test_<name>.cxx or explicitly as given by <source>

    '''

    from waflib.Tools import waf_unit_test
    bld.add_post_fun(waf_unit_test.summary)

    source = to_list(source)
    print source
    target = osp.splitext(str(source[0]))[0]

    bld.program(features='test', source=source, target=target,
                includes = _includes(bld, includes))
