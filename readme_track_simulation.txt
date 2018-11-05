README file for a MCS simulation ... 

WCP
1. wire-cell-mcs-sim
2. rm -f mcssim.json.bz2
3. bzip2 mcssim.json


Now go to WCT
4. wire-cell -c pgrapher/experiment/uboone/wct-jsondepo-sim-nf-sp.jsonnet


Now go back to WCP
5. wire-cell-imaging-lmem-celltree ../ChannelWireGeometry_v2.txt celltree.root -s2
6. dev-wire-cell-tracking-fit ../ChannelWireGeometry_v2.txt result_0_0_0.root
7. root -l plot_track_proj.C
8. cd bee; python dump_json.py ../tracking_0_0_0.root cluster simple charge deadarea


