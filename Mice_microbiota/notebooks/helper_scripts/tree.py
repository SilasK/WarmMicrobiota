import ete3

import numpy as np


def load_tree(netwik_file):
    return ete3.Tree(netwik_file,quoted_node_names=True,format=1)

def root_tree_by_phyla(T,phyla):
    """ Root the tree next to the least frequent phyla if possible

    """


    Freq_pyla= phyla.value_counts()

    for p in reversed(Freq_pyla.index):
        LCA = T.get_common_ancestor(*tuple(phyla.index[phyla==p].values))

        if not T== LCA:
            T.set_outgroup(LCA)
            print(f"set {p} as outgroup for Tree rooting")
            break


    T.unroot()

def layout_black_circles(node):
    # If node is a leaf
    if node.is_leaf():
        node.img_style["fgcolor"]='k'
    else:
        node.img_style["size"]=0

def render_tree(T,out):
    from ete3 import TreeStyle

    ts = TreeStyle()
    ts.show_leaf_name= False
    ts.mode = "c"
    ts.scale=200
    ts.show_scale=False

    T.render(out,tree_style=ts,layout=layout_black_circles)

def tree2linkage(T,species=None,linkage_method='complete'):

    #import scipy.spatial as sp,
    import scipy.cluster.hierarchy as hc
    from itertools import combinations

    if species is None:
        species = T.get_leaf_names()


    dist=np.array([T.get_distance(*pair) for pair in combinations(species,2)])
    linkage = hc.linkage(dist, method=linkage_method)

    return linkage


def taxonomy2linkage(T,species=None,linkage_method='average'):

    import scipy.cluster.hierarchy as hc
    from itertools import combinations

    if species is None:
        species = T.index

    def dist_taxonomy(T,i,j):
        return T.shape[1]-sum(T.loc[i]==T.loc[j])



    dist=np.array([dist_taxonomy(T,*pair) for pair in combinations(species,2)])
    linkage = hc.linkage(dist, method=linkage_method)

    return linkage
