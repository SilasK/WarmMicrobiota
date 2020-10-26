
import matplotlib.pylab as plt
import numpy as np

import pandas as pd
import seaborn as sns
import os,sys
import helper_scripts as hs
from helper_scripts import effect_plot as EP

from helper_scripts import DimensionalReduction

import yaml



def load_metadata():

    metadata= pd.read_csv('metadata.tsv',index_col=0,sep='\t')

    return metadata

name_mapping= {'SH-H':'34°C','SH-R':'RT','OVA-R':'ova RT','OVA-H':'ova 34°C',
                   'OVA_t-RT':'ova\ntranspl.\nRT','OVA_t-H':'ova\ntranspl.\n34°C','trPF':'trPF','PF':'PF',
                   'trH':'34°C\ntranspl.', 'trRT':'RT\ntranspl.', 'H':'34°C','RT':'RT'
                  }

def rename_metadata(data,metadata):




    for i,row in metadata.iterrows():
        metadata.loc[i,'SampleID']= row.SampleID.replace(row.Group,name_mapping[row.Group]+' ')

    metadata['Group'] = metadata.Group.map(name_mapping)


    data.index= data.index.map(metadata.SampleID)

    metadata.index = metadata.index.map(metadata.SampleID)


def take_subset(data,metadata,subset_groups=None):

    if not subset_groups is None:

        assert all([g in metadata.Group.unique() for g in subset_groups] )

        metadata= metadata.loc[metadata.Group.isin(subset_groups)].copy()


    intersection= metadata.index.intersection(data.index)
    assert len(intersection)>2

    metadata=metadata.loc[intersection].copy()
    data= data.loc[intersection].copy()

    return data,metadata


def load_aldex(aldex_file,get_sample_data=True):

    Aldex= pd.read_csv(aldex_file,index_col=0,sep='\t')

    sample_columns= Aldex.columns.str.startswith('rab.sample')

    Stats=  Aldex.loc[:,~sample_columns].copy()

    if get_sample_data:
        if not sample_columns.any():
            raise Exception('No sample data in Aldex stats table')

        data= Aldex.loc[:,sample_columns].copy()
        data=data.T
        data.index= data.index.str.replace('rab.sample.','')

        return Stats,data
    else:
        return Stats




#Tax= pd.read_table('../taxonomy/Silva.tsv',index_col=0)
Tax= pd.read_csv('../data/Taxonomy_silva_13_2.tsv',index_col=0,sep='\t')
Tax.columns= Tax.columns.str.lower()

OTUnames= pd.Series(data=hs.microbiota.gen_names_for_range(Tax.shape[0],'otu'),
               index=Tax.index)

Labels=  Tax.ffill(axis=1).genus +' '+OTUnames


#Tax= pd.read_table('../taxonomy/mapseq.tsv',index_col=0,skiprows=1)
#Tax.columns= Tax.columns.str.lower()
#Labels= Tax[['family','species']].fillna('').apply(lambda s: ' '.join(s),axis=1)


# Richness

import numpy as np
from numpy.random import RandomState
import skbio

def rarefaction(M, seed=None):
    prng = RandomState(seed) # reproducible results
    noccur = np.sum(M, axis=1) # number of occurrences for each sample
    nvar = M.shape[1] # number of variables
    depth = np.min(noccur) # sampling depth

    Mrarefied = np.empty_like(M)
    for i in range(M.shape[0]): # for each sample
        p = M[i] / float(noccur[i]) # relative frequency / probability
        choice = prng.choice(nvar, depth, p=p)
        Mrarefied[i] = np.bincount(choice, minlength=nvar)


    return Mrarefied



def abundance_FC_plot(Stats,sig=None,ax=None):
    if ax is None:
        ax= plt.subplot(111)

    if sig is not None:
        D= Stats.copy()
        D['label']='Not significant'
        D.loc[sig,'label']='Significant'
    else:
        D= Stats

    sns.set_palette(['orange','darkgrey'])

    sns.scatterplot(data=D,x='rab.all',y='effect',hue='label',ax=ax,hue_order=['Significant','Not significant'])

    x= np.array(ax.get_xlim())
    ax.plot(x,[0,0],'k-')
    #ax.plot(x ,[-1,-1],'k--')
    ax.set_xlabel("Median relative abundance [$\log_2$]")
    ax.set_ylabel("Effect size [$\log_2$]")
    ax.get_legend().remove()


def multiplot(V,n_rows,n_cols,subset=None,sharey=True,sharex=False,figsize=None,sig_labels_params=None,pannel_letters=True,**kws):
    """
    Plots multiple variables from the Viewpont as a multi pannel figure, define subset, n_col and n_rows
    """

    N=  n_cols*n_rows

    if subset is None:
        subset= V.Data.columns[:N]

    f,axe = plt.subplots(n_rows,n_cols,sharey=sharey,sharex=sharex,figsize=figsize,constrained_layout=True)
    if sig_labels_params is None:
        sig_labels_params={}
        if sharey:
            max_value= V.Data[subset].max().max()
            sig_labels_params=dict(y0=max_value,deltay=max_value/10.)

    axe= np.ravel(axe)

    for i in range(N):
        V.Boxplot(subset[i],ax= axe[i],sig_labels_params=sig_labels_params,**kws)


    if pannel_letters:
        import string
        for i,ax in enumerate(axe):
            annotate_subplot(ax,string.ascii_uppercase[i])
    return axe




def Effect_size_plot(Stats,ab_treshold=0, p_treshold = 0.05 ,effect='effect' ,p_value_method='we.eBH'):


    sig= (Stats[p_value_method]<p_treshold)#&(Stats['rab.all']>1)&(Stats[effect].abs()>1)


    assert effect=='effect', 'plot is made for effect size'
    abundance_FC_plot(Stats,sig.index[sig])

    ax= plt.gca()
    plt.vlines(ab_treshold,*ax.get_ylim(),linestyles='dashed')

    if 'BH' in p_value_method:
        Pvalue_method_name="P_{BH}"
    else:
        Pvalue_method_name='P'

    ax.text(0.7,0.1,f"${Pvalue_method_name} < {p_treshold}$",transform=ax.transAxes)



    sig = sig&(Stats['rab.all']>ab_treshold)

    sig = Stats.loc[sig].sort_values(effect,ascending=False).index

    return sig
