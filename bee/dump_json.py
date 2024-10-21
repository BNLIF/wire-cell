#!/usr/bin/env python
import os, sys, glob, shutil

ALIAS = {
    'true' : 'truth',
    'simple'  : 'rec_simple',
    'charge' : 'rec_charge_blob',
    'deblob' : 'rec_charge_cell',
    'mc' : 'mc',
    'deadarea' : 'channel-deadarea',
    'flash' : 'op',
    'cluster' : 'cluster'
}

def main(filename, options):
    if (options == []):
        options =["rec_charge_blob", "rec_simple", "truth", "mc"]
    for i in range(len(options)):
        options[i] = ALIAS.get(options[i], options[i])
    # print options
    to_glob = filename[:filename.rfind('_')] + '_*.root'
    # print to_glob
    list_of_files = glob.glob(to_glob)
    list_of_files.sort(key=lambda x: int(x[x.rfind('_')+1:x.rfind('.')]))
    # print list_of_files
    # root -b -q 'WireCell2JSON.C("shower3D_data_0.root","rec_simple", "1.json")'
    try:
        os.mkdir('data')
    except OSError:
        shutil.rmtree('data')
        print('removing old data directory ....')
        os.mkdir('data')
    for i in range(len(list_of_files)):
        str_i = str(i)
        os.mkdir('data/'+str_i)
        filename = list_of_files[i]
        for rec in options:
            cmd = "root -b -q -l loadClasses.C 'WireCell2JSON.C(\"%s\", \"%s\", \"%s\")'" % (
                filename,
                rec,
                'data/'+str_i+'/'+str_i+'-'+rec+'.json')
            # print cmd
            os.system(cmd)
        # trackfilename = list_of_files[i].replace(
        #     'shower3D',
        #     'cluster'
        #     )
        # if (os.path.exists(trackfilename)):
        #     cmd = "root -b -q -l 'Track2JSON.C(\"%s\", \"%s\")'" % (
        #         trackfilename,
        #         'data/'+str_i+'/'+str_i+'-WireCell-charge-track.json')
        #     #print cmd
        #     #os.system(cmd)

    if (os.path.exists('to_upload.zip')):
        print('removing old to_upload.zip ...')
        os.system('rm to_upload.zip')

    os.system('zip -r to_upload data')

def usage():
    print("""
    python dump_json.py [filename] [alg1 alg2 ...]

        available algorithms: simple, charge, true, deblob, mc, deadarea, cluster
    """)

if __name__ == "__main__":
    if (len(sys.argv)<=1):
        usage()
    else:
        main(sys.argv[1], sys.argv[2:])
