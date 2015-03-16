import os
import os.path as osp
import waflib
import waflib.Utils
import waflib.Logs as msg
from waflib.Configure import conf
from waflib.TaskGen import feature, before_method, after_method, extension, after

_tooldir = osp.dirname(osp.abspath(__file__))


def options(opt):
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
    cfg.find_program('rootcling', var='ROOTCLING', path_list=path_list)
    cfg.find_program('rootcint', var='ROOTCINT', path_list=path_list)
    cfg.find_program('rlibmap', var='RLIBMAP', path_list=path_list, mandatory=False)
    return

@conf
def gen_rootcling_dict(bld, name, linkdef, headers = '', includes = ''):
    '''
rootcling -f dictToto.cxx -rml libtoto.so -rmf libtoto.rootmap myHeader1.h myHeader2.h ... LinkDef.h
    '''
    headers = waflib.Utils.to_list(headers)
    incs = list()
    for maybe in waflib.Utils.to_list(includes):
        if maybe.startswith('/'):
            incs.append('-I%s' % maybe)
        else:
            incs.append('-I%s' % bld.path.find_dir(maybe).abspath())
    incs = ' '.join(incs)
    print 'INCS:',incs
    
    dict_src = name + 'Dict.cxx'
    dict_lib = 'lib' + name + 'Dict.so' # what for Mac OS X?
    dict_map = 'lib' + name + 'Dict.rootmap'

    rule = '${ROOTCLING} -f ${TGT[0]} -rml %s -rmf ${TGT[1]} %s ${SRC}' % (dict_lib, incs)
    print 'RULE:',rule
    bld(source = headers + [linkdef],
        target = [dict_src, dict_map],
        rule=rule)

    bld.shlib(source = dict_src,
              target = name+'Dict',
              includes = includes,
              use = 'ROOTSYS')

    bld.install_files('${PREFIX}/lib/', dict_map)


@conf
def gen_rootcint_dict(bld, name, linkdef, headers = '', includes=''):
    '''Generate a rootcint dictionary, compile it to a shared lib,
    produce its rootmap file and install it all.
    '''
    headers = waflib.Utils.to_list(headers)
    incs = ['-I%s' % bld.path.find_dir(x).abspath() for x in waflib.Utils.to_list(includes)]
    incs = ' '.join(incs)
    
    dict_src = name + 'Dict.cxx'

    bld(source = headers + [linkdef],
        target = dict_src,
        rule='${ROOTCINT} -f ${TGT} -c %s ${SRC}' % incs)

    bld.shlib(source = dict_src,
              target = name+'Dict',
              includes = includes,
              use = 'ROOTSYS')

    rootmap = 'lib%sDict.rootmap'%name
    bld(source = [linkdef], target=rootmap,
        rule='${RLIBMAP} -o ${TGT} -l lib%sDict.so -d lib%s.so -c ${SRC}' % (name, name))
    
    bld.install_files('${PREFIX}/lib/', rootmap)
                      
