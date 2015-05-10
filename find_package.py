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
from waflib.Logs import debug, info, error

_tooldir = osp.dirname(osp.abspath(__file__))

def options(opt):
    opt.load('compiler_cxx')
    opt.load('waf_unit_test')
    opt.load('find_root', tooldir=_tooldir)
    opt.load('find_eigen3', tooldir=_tooldir)
    opt.load('boost', tooldir=_tooldir)
    

def configure(cfg):
    cfg.load('compiler_cxx')
    cfg.load('waf_unit_test')
    cfg.load('find_root', tooldir=_tooldir)
    cfg.load('find_eigen3', tooldir=_tooldir)
    cfg.load('boost', tooldir=_tooldir)
    pass

def build(bld):
    from waflib.Tools import waf_unit_test
    bld.add_post_fun(waf_unit_test.summary)


@conf
def make_package(bld, name, use=''):
    use = to_list(use) + ['ROOTSYS']

    includes = []
    headers = []
    source = []

    incdir = bld.path.find_dir('inc')
    srcdir = bld.path.find_dir('src')
    dictdir = bld.path.find_dir('dict')
    testsrc = bld.path.ant_glob('test/test_*.cxx') + bld.path.ant_glob('tests/test_*.cxx')
    appsdir = bld.path.find_dir('apps')

    if incdir:
        headers += incdir.ant_glob(name + '/*.h')
        includes += ['inc']
        bld.env['INCLUDES_'+name] = [incdir.abspath()]

    if headers:
        bld.install_files('${PREFIX}/include/%s' % name, headers)

    if srcdir:
        source += srcdir.ant_glob('*.cxx')

    if incdir and srcdir:
        bld(features = 'cxx cxxshlib',
            name = name,
            source = source,
            target = name,
            includes = 'inc',
            export_includes = 'inc',
            use=use)

    if dictdir:
        if not headers:
            error('No header files for ROOT dictionary "%s"' % name)
        linkdef = dictdir.find_resource('LinkDef.h')
        bld.gen_rootcling_dict(name, linkdef,
                               headers = headers,
                               includes = includes, 
                               use = use)
    if testsrc:

        for test_main in testsrc:
            bld.program(features = 'test', 
                        source = [test_main], 
                        target = test_main.name.replace('.cxx',''),
                        install_path = None,
                        includes = 'inc',
                        use = use + [name])
    if appsdir:
        for app in appsdir.ant_glob('*.cxx'):
            bld.program(source = [app], 
                        target = app.name.replace('.cxx',''),
                        includes = 'inc',
                        use = use + [name])

    
