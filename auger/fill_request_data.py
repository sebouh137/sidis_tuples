import sys

inputfiles=''
args =''
for arg in sys.argv:
    if arg.startswith('--'):
        args+= arg+" "
    else :
        inputfiles+=arg+" "
        
with open('tmp.xml', 'w') as ofile:
    with open('make_tuples_template.xml', 'r') as ifile:
        s = ifile.read()
        s = s.replace("INPUT_FILES", inputfiles).replace("ARGS", args)
        
        ofile.write(s)
