#!/usr/bin/perl

# Get event number from command line if provided
my $targetEvent = $ARGV[0] if @ARGV;

# process the result.root files ... 
# for (my $i = 0; $i!=35; $i ++){
#     if ($i % 8 == 7){
#         system("wire-cell-imaging-lmem-celltree ../ChannelWireGeometry_v2.txt celltreeOVERLAY.root $i -d1 -s2");
#     }else{
#         system("wire-cell-imaging-lmem-celltree ../ChannelWireGeometry_v2.txt celltreeOVERLAY.root $i -d1 -s2 &");
#     }
# }


# process the clustering files ...
open(infile,"filelist_5384");
my $i = 0;
while(<infile>){
    my $filename = $_;
    chomp($filename);
    $filename =~ /nuselEval_(\d+)_\d+_(\d+)\.root/;
    my $runNo = $1;
    my $eventNo = $2;

    # Skip if target event is specified and doesn't match
    next if (defined $targetEvent && $eventNo != $targetEvent);

    if ($i%36 == 35){
        #system("dev-wire-cell-clustering-port ./input_data_files/ChannelWireGeometry_v2.txt $filename -b1 >$runNo\_$eventNo\.log");
        #system("prod-wire-cell-matching-nusel-port ./input_data_files/ChannelWireGeometry_v2.txt $filename -d1 >$runNo\_$eventNo\.log");
        #system("prod-wire-cell-matching-nusel ./input_data_files/ChannelWireGeometry_v2.txt $filename -d1 >$runNo\_$eventNo\.log");
        #system("wire-cell-prod-stm-port ./input_data_files/ChannelWireGeometry_v2.txt $filename 0 -d0 -o1 -g2 >$runNo\_$eventNo\.log");
        system("wire-cell-prod-nue-port ./input_data_files/ChannelWireGeometry_v2.txt $filename 0 -d0 -o1 -ginit_first_segment >$runNo\_$eventNo\.log");
    }else{
        #system("dev-wire-cell-clustering-port ./input_data_files/ChannelWireGeometry_v2.txt $filename -b1 >$runNo\_$eventNo\.log&");
        #system("prod-wire-cell-matching-nusel-port ./input_data_files/ChannelWireGeometry_v2.txt $filename -d1 >$runNo\_$eventNo\.log &");
        #system("prod-wire-cell-matching-nusel ./input_data_files/ChannelWireGeometry_v2.txt $filename -d1 >$runNo\_$eventNo\.log &");
        #system("wire-cell-prod-stm-port ./input_data_files/ChannelWireGeometry_v2.txt $filename 0 -d0 -o1 -g2 >$runNo\_$eventNo\.log &");
        #system("wire-cell-prod-nue ./input_data_files/ChannelWireGeometry_v2.txt $filename 0 -d0 -o1 -g2 >$runNo\_$eventNo\.log &");
        system("wire-cell-prod-nue-port ./input_data_files/ChannelWireGeometry_v2.txt $filename 0 -d0 -o1 -ginit_first_segment >$runNo\_$eventNo\.log &");
    }
    $i++;
}

# Exit if no matching event was found when target event was specified
if (defined $targetEvent && $i == 0) {
    print "No matching event found for event number: $targetEvent\n";
    exit 1;
}