import numpy as np
import matplotlib
import matplotlib.pylab as plt
import pandas as pd
import seaborn as sns

from .plotting import format_p_value


from helper_scripts.plotting import format_p_value

def effect_plot(Values,P_values,colors=None,width=0.7,Labels=None, Title= None,ax=None):

    sns.set_palette(colors)

    from matplotlib.font_manager import FontProperties
    font=FontProperties(size='x-small')


    #assert Values.isnull().any().any()==False

    if colors is None:
        colors= sns.color_palette('deep',n_colors= Values.shape[1])




    ##
    d1,d2= Values.shape
    if ax is None:
        f,ax= plt.subplots(1,1,figsize=(10,d1/width/10*(d2+1)+1))

    assert P_values.shape[1] == d2
    Values= Values.loc[reversed(Values.index)]
    P_values= P_values.loc[Values.index]


    Values.plot.barh(width=width,ax=ax)

#

    Horizontal_alignement={1:'left',-1:'right'}

    for i in range(d1):
        for j in range(d2):

            y= Values.iloc[i,j]
            text= format_p_value(P_values.iloc[i,j],Stars=True)
            if not text =='n.s.':
                ax.annotate(text, xy=(y,i+width/d2*(0.5+j)- width/2),textcoords='offset points',xytext=(np.sign(y),-1),
                    horizontalalignment= Horizontal_alignement[np.sign(y)], verticalalignment='center',font_properties=font)

    if Labels is not None:
        if type(Labels) is pd.Series:
            ax.set_yticklabels(Labels.loc[Values.index])
        else:
            ax.set_yticklabels(Labels)
    if Title is not None:
        ax.set_title(Title)

    handles,labels = ax.get_legend_handles_labels()

    ax.legend(reversed(handles),reversed(labels))

    return ax





def aldex_plot(D):

    f,axe = plt.subplots(1,2,figsize=(15,5))

    sns.set_palette(['orange','darkgrey'])

    sns.scatterplot(data=D,x='rab.all',y='effect',hue='label',ax=axe[0],hue_order=['Significant','Not significant'])

    x= np.array(axe[0].get_xlim())
    axe[0].plot(x,[1,1],'k--')
    axe[0].plot(x ,[-1,-1],'k--')
    axe[0].set_xlabel("Median $\log_2$ relative abundance")
    axe[0].set_ylabel("Effect size [$\log_2$]")
    axe[0].get_legend().remove()
    sns.scatterplot(data=D,x='diff.win',y='diff.btw',hue='label',ax=axe[1],hue_order=['Significant','Not significant'])

    x= np.array(axe[1].get_xlim())
    axe[1].plot(x,x,'k--')
    axe[1].plot(x ,-x,'k--')
    axe[1].set_xlabel("Median $\log_2$ win condition diff")
    axe[1].set_ylabel("Median $\log_2$ btw condition diff")
    axe[1].legend(*axe[1].get_legend_handles_labels(),bbox_to_anchor=(1, 1))

#plt.savefig(name)
