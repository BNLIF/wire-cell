def options(opt):
    pass
    
def configure(cfg):
    # likely need to set environment variable: PKG_CONFIG_PATH
    cfg.check_cfg(package='eigen3',  uselib_store='EIGEN', args='--cflags --libs', mandatory=False)
    cfg.check(header_name="Eigen/Dense", use='EIGEN', mandatory=False)
