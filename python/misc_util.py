import pandas as pd


#mode options are
# df:  return the slice dataframes (default)
# len: return the numbers of entries in each slice
class BinIterator:
    def __init__(self,df,xvar,min,max,bins,mode='df'):
        self._df = df
        self._xvar = xvar
        self._min = min
        self._max = max
        self._bins = bins
        self._i = 0
        self._mode = mode
    def __iter__(self):
        self._i = 0
        return self
    def __next__(self):
        if(self._i >= self._bins):
            raise StopIteration
        dx = (self._max-self._min)/self._bins
        #case:  single dataframe
        if type(self._df)==pd.DataFrame:
            ret = self._df.query("%s >= %s and %s < %s"\
                       %(self._xvar,
                         self._min + self._i*dx,
                         self._xvar,
                         self._min +(self._i+1)*dx))
            if self._mode == 'len':
                ret = len(ret)
        #case:  list of dataframes
        else: 
            ret = [df.query("%s >= %s and %s < %s"\
                       %(self._xvar,
                         self._min + self._i*dx,
                         self._xvar,
                         self._min +(self._i+1)*dx)) for df in self._df]
            if self._mode == 'len':
                ret = [len(r) for r in ret]
        self._i+=1
        return self._min+(self._i-1/2)*dx, ret


def query_or_all(df,q):
    if q != "":
        return df.query(q)
    else:
        return df
