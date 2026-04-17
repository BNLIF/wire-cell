#!/usr/bin/perl

# process the result.root files ... 
for (my $i = 0; $i!=35; $i ++){
    if ($i % 8 == 7){
        system("wire-cell-imaging-lmem-celltree ./input_data_files/ChannelWireGeometry_v2.txt celltree_intrinsic_nue/celltreeOVERLAY.1649126904.root $i -d1 -s2");
    }else{
        system("wire-cell-imaging-lmem-celltree ./input_data_files/ChannelWireGeometry_v2.txt celltree_intrinsic_nue/celltreeOVERLAY.1649126904.root $i -d1 -s2 &");
    }
}