import sys

inputfiles=''
args =''
for arg in sys.argv[1:]:
    if arg.startswith('--'):
        args+= arg+" "
    else :
        inputfiles+=arg+" "

#remove the base directory name
inputfiles = inputfiles.replace("/lustre/expphy/cache/clas12/rg-a/production/recon/fall2018/","")
with open('tmp.xml', 'w') as ofile:
    with open('make_tuples_template.xml', 'r') as ifile:
        s = ifile.read()
        s = s.replace("INPUT_FILES", inputfiles).replace("ARGS", args)
        s = s.replace("QADB", "/home/spaul/clasqaDB/qadb/qa.rga_outbending/qaTree.json" if 'torus+1' in inputfiles else "/home/spaul/clasqaDB/qadb/qa.rga_inbending/qaTree.json" )
        ofile.write(s)
