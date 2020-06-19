
import matplotlib.pylab as plt
import matplotlib
import os
import numpy as np




def plotting_Setup(sns_contex='paper',font_scale=1.5, save_dpi=150):
    import matplotlib
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_context(sns_contex,font_scale=font_scale)
    matplotlib.rcParams['savefig.dpi']=save_dpi
    matplotlib.rcParams['pdf.fonttype']=42


def saveplot(name,figurefolder='../Figures/',formats=['.pdf','.png'],SAVEPLOT=True):
    matplotlib.rcParams['pdf.fonttype']=42 # to save as vector format
    if SAVEPLOT:
        if not os.path.exists(figurefolder):
            os.makedirs(figurefolder)
        if len(formats)==0:
            raise Exception('no format specified')

        for format in formats:
            plt.savefig(os.path.join(figurefolder,name+format),bbox_inches='tight')




    

