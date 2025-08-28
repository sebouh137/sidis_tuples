import uproot3,pandas as pd

#create single-tree root file
def to_root(df,filename, treename):
    with uproot3.recreate(filename) as f:
        f[treename] = uproot3.newtree({col:df[col].dtype for col in df.columns})
        f[treename].extend(dict(df))

#create root file with multiple trees:
# d is a map of tree name to data frames {"tree1":df1, "tree2":df2 ... etc}
def to_root_multi(filename, d):
    with uproot3.recreate(filename) as f:
        for treename in d.keys():
            df = d[treename]
            f[treename] = uproot3.newtree({col:df[col].dtype for col in df.columns})
            f[treename].extend(dict(df))
pd.DataFrame.to_root = to_root

def read_root(filename,treename=None,N=None):
    if type(filename) != str:
        return pd.concat([read_root(f,treename,N) for f in filename])
        
    with uproot3.open(filename) as f:
        if treename is None:
            if len(f.keys()) == 1:
                treename = f.keys()[0]
            elif len(f.keys()) == 0:
                raise Exception("error: no tree names in file " + filename)
            else:
                raise Exception("error: treename must be specified; more than one tree in " + filename)
            
        df = f[treename].pandas.df(entrystop=N)
    return df
if __name__=='__main__':
    df = pd.DataFrame({'a':[1.0,2.2,3.3],'b':[0,1,2]})
    #to_root(df,'out.root', 'tree')
    df.to_root('out.root', 'tree')
    df = read_root('out.root')
    print("reading a one tree file without specifying tree name\n",df)

    df1 = pd.DataFrame({'a':[1.0,2.2,3.3],'b':[3,2,1]})
    df2 = pd.DataFrame({'c':[1.3,1.2,3.3],'d':[0,1,2]})
    to_root_multi('out.root', {'tree1':df1, 'tree2':df2})
    df1 = read_root('out.root',"tree1")
    print("reading tree1\n",df1)
    df2 = read_root('out.root',"tree2")
    print("reading tree2\n",df2)
