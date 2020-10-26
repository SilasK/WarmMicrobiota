import pandas as pd
import numpy as np

working_dir='../WD'

# get LAabels
DT= pd.read_table(f'{working_dir}/genomes/taxonomy/taxonomy_names.tsv')
DT= DT.drop_duplicates(['# bin']).sort_values('# bin')

DT.index= DT['# bin']

Tax= DT.loc[:,'superkingdom':].fillna('not classified')
Tax= Tax.applymap(lambda s: s.split(':')[0])
Tax= Tax.replace('not classified',np.nan)
Tax=Tax.fillna(method='ffill',axis=1)

Labels=Tax.species+' '+Tax.index

# Labels= checkmTax.fillna(method='ffill',axis=1).Genus+' '+checkmTax.index

#from adjustText import adjust_text

def MA_plot(abundance,change,sig=None,ax=None):
    
    if ax is None:
        ax = plt.subplot(111)
        
    abundance= pd.Series(abundance)
    change= pd.Series(change)

    sns.scatterplot(abundance,change,
                    color='g',ax=ax)
    if sig is not None:
        sns.scatterplot(abundance.loc[sig],change.loc[sig],color='r',ax=ax)
    
    ax.set_xlabel('Abundnace')
    ax.set_ylabel('FC [log2]')
    ax.hlines(0,*ax.get_xlim())
    
    #ax.fill_between(ax.get_xlim(),-1,1,alpha=0.5,color='grey')