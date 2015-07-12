#!/bin/bash

# A little helper to run an individual test without requiring a full
# install.  It runs the test from the build directory.
#
# You problably want to do:
#   cd wire-cell && alias run_test=`pwd`/run_test.sh

usage () {
    echo "run_test.sh testname [package]"
    echo "if run from a package directory last arg isn't needed"
    exit 1;
}


tst="$1" ; shift
pkg="$1" ; shift
if [ -z "$pkg" ] ; then
    pkg=$(basename $(pwd))
fi

top=$(dirname $(readlink -f $BASH_SOURCE))

dir="$top/build/$pkg"
if [ ! -d $dir ] ; then
    echo "No build directory: $dir"
    usage
fi

exe="$dir/test_$tst"
if [ ! -x $exe ] ; then
    echo "Test $tst not built"
    usage
fi

LD_LIBRARY_PATH=$dir:$LD_LIBRARY_PATH $exe


