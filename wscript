#!/usr/bin/env python

TOP = '.'
APPNAME = 'WireCell'

subdirs = ['data','nav','sst']

def options(ctx):
    ctx.load('find_package')

def configure(ctx):
    ctx.load('find_package')
    # ctx.env.LIBPATH_WireCell = ctx.options.prefix + '/lib'
    # ctx.env.INCLUDES_WireCell = ctx.options.prefix + '/include'
 
def build(bld):
    bld.recurse(subdirs)


