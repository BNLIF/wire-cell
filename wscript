#!/usr/bin/env python

TOP = '.'
APPNAME = 'WireCell'

subdirs = ['data',
           'nav',
           'sst',
           'tiling',
           'examples',
           'matrix', 
           '2dtoy',
]

def options(opt):
    opt.load('doxygen', tooldir='waf-tools')
    opt.load('find_package', tooldir='waf-tools')
    opt.add_option('--build-debug', default='-O2',
                   help="Build with debug symbols")
    opt.add_option('--doxygen-tarball', default=None,
                   help="Build Doxygen documentation to a tarball")

def configure(cfg):
    cfg.load('doxygen', tooldir='waf-tools')
    cfg.load('find_package', tooldir='waf-tools')
    
    cfg.env.CXXFLAGS += [cfg.options.build_debug]

def build(bld):
    bld.load('find_package', tooldir='waf-tools')
    bld.recurse(subdirs)
    if bld.env.DOXYGEN and bld.options.doxygen_tarball:
        bld(features="doxygen", doxyfile=bld.path.find_resource('Doxyfile'),
            doxy_tar = bld.options.doxygen_tarball)

