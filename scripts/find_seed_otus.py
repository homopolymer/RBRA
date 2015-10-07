#!/usr/bin/env python
"""
Find seed OTUs based on the phylogeny tree and abundance profile

Change Log
==========
Sep 25, 2015    Feng Zeng    Add coverage criterion

"""

import os
import sys
import argparse
import logging
import textwrap
import traceback
import time
import csv
import itertools
from sets import Set
from collections import defaultdict

from ete2 import Tree

import numpy as np
from scipy.cluster import hierarchy



gene_set_dt = np.dtype([('gene','S1024'),('abundance','f8')])



def load_gene_abundance(filename):
    '''load gene abundances from the given file

    Parameters
    ----------
    filename : gene abundance profile file

    Returns
    -------
    dict
        a dictionary object storing gene abundance

    '''
    gene_abun = defaultdict(float)
    with open(filename) as f:
        r = csv.reader(f,delimiter='\t')
        for row in r:
            gene_abun[row[0]] = float(row[3])
    return gene_abun

def load_gene_coverage(filename):
    '''load gene coverage from the given file

    Parameters
    ----------
    filename : gene abundance profile file

    Returns
    -------
    dict
        a dictionary object storing gene coverage

    '''
    gene_cover = defaultdict(float)
    with open(filename) as f:
        r = csv.reader(f,delimiter='\t')
        for row in r:
            gene_cover[row[0]] = float(row[4])
    return gene_cover


def load_phylogeny_tree(filename):
    '''load phylogeny tree from the given newick file

    Paramters
    ---------
    filename : phylogene tree file

    Returns
    -------
    tree
        a tree object

    '''
    tree = Tree(filename)
    
    # add a numeric id to tree node
    count = 0
    for node in tree.traverse('postorder'):
        node.add_feature('id',count)
        count += 1

    return tree


def gene_tree_cluster(tree,dissim_thres):
    '''
    cluster tree nodes based on sequence similarity
    '''
    if tree.is_leaf():
        tree.add_feature('centroid',tree.gene_set)
    else: 
        # get the children
        children = tree.get_children()

        # get child cluster size
        child_cluster_size = [len(child.centroid) for child in children]

        # skip if a child has multiple clusters
        if np.sum(np.asarray(child_cluster_size)>1) > 0:
            centroid = np.hstack(tuple(child.centroid for child in children))
            tree.add_feature('centroid',centroid)
            return None
        
        # create a zero distance matrix
        cluster_size = np.sum(child_cluster_size)
        dist = np.zeros((cluster_size,cluster_size))

        # fill out the up diagonal matrix
        for c1,c2 in itertools.combinations(np.arange(len(children)),2):
            d1 = int(np.sum(child_cluster_size[:c1]))  # position shift
            d2 = int(np.sum(child_cluster_size[:c2]))  # position shift
            k1 = 0
            k2 = 0
            # centroid in cluster k1 of node c1
            g1 = children[c1].centroid[k1]['gene']
            # centroid in cluster k2 of node c2
            g2 = children[c2].centroid[k2]['gene']
            # calculate the distance between cluster k1 and cluster k2
            dist_g1_g2 = tree.get_distance(str(g1),str(g2))
            # save to distance matrix
            dist[d1+k1][d2+k2] = dist_g1_g2

        # preform hierarchical clustering on the distance matrix
        Z = hierarchy.linkage(dist)
        C = hierarchy.fcluster(Z,t=dissim_thres,criterion='distance')-1  # make sure starting from 0

        # update tree data
        centroid = np.array([(None,None)]*len(np.unique(C)),dtype=gene_set_dt)

        for c in xrange(len(children)):
            d = int(np.sum(child_cluster_size[:c]))  # position shift
            for i in xrange(child_cluster_size[c]):
                k = C[d+i]
                if np.isnan(centroid[k]['abundance']):
                    centroid[k] = children[c].centroid[i]
                else:
                    if centroid[k]['abundance'] < children[c].centroid[i]['abundance']:
                        centroid[k]['gene'] = children[c].centroid[i]['gene']
                    centroid[k]['abundance'] += children[c].centroid[i]['abundance']

        tree.add_feature('centroid',centroid)


def load_gene_taxonomy(filename):
    '''load gene taxonomy annotation
    
    Parameters
    ----------
    filename : taxonomy annotation file
    
    Returns
    -------
    Dict
        a dictonary object storing gene taxonomy annotation
    '''
    gene_tax = defaultdict()
    with open(filename,'r') as f:
        for line in f:
            gene,tax = line.rstrip().split('\t')
            gene_tax[gene] = tax
    return gene_tax


def find_seed_otus():
    '''find seed otus based on phylogeny and gene abundance

    Parameters
    ----------

    Returns
    -------

    '''
    global opts

    # TODO load phylogeny tree
    if opts.verbose:
        logging.info('load phylogeny tree file')
    t0 = time.time()
    gene_tree = load_phylogeny_tree(os.path.abspath(opts.tree[0]))
    if opts.verbose:
        logging.info('elapsed time on loading phylogeny tree is %f minutes' % ((time.time()-t0)/60.))

    # TODO load gene abundance
    if opts.verbose:
        logging.info('load gene abundance file')
    t0 = time.time()
    gene_abun = load_gene_abundance(os.path.abspath(opts.abun[0]))
    if opts.verbose:
        logging.info('elapsed time on loading gene abundances is %f minutes' % ((time.time()-t0)/60.))

    # TODO load gene coverage
    gene_cover = load_gene_coverage(os.path.abspath(opts.abun[0]))
   
    # TODO load taxonomy annotation
    if opts.taxonomy is not None:
        if opts.verbose:
            logging.info('load gene taxonomy annotation')
        t0 = time.time()
        gene_tax = load_gene_taxonomy(os.path.abspath(opts.taxonomy))
        if opts.verbose:
            logging.info('elapsed time on loading gene taxonomy annotation is %f minutes' % ((time.time()-t0)/60.))

    # TODO add new attributes to tree nodes
    for node in gene_tree.iter_leaves():
        if node.name in gene_abun:
            node.add_feature('gene_set',np.array([(node.name,gene_abun[node.name])],dtype=gene_set_dt))
        else:
            node.add_feature('gene_set',np.array([(node.name,0)],dtype=gene_set_dt))

    # TODO propogate gene cluster from leaves to root
    if opts.verbose:
        logging.info('propogate gene cluster from bottom to top based on similarity')
    t0 = time.time()
    for t in gene_tree.traverse('postorder'):
        gene_tree_cluster(t,1.-opts.sim_thres[0])
    if opts.verbose:
        logging.info('elapsed time on propogating gene clustering is %f minutes' % ((time.time()-t0)/60.))
   
    # TODO compute total abundance
    total_abun = 0
    for node in gene_tree.iter_leaves():
        total_abun += gene_abun[node.name]
    abun_thres = opts.depth_thres[0]
    if opts.depth_ratio is not None:
        abun_thres = opts.depth_ratio*total_abun

    # TODO report seed OTUs
    if opts.verbose:
        logging.info('find gene cluster with depth threshold %f and coverage threshold %f' % (abun_thres,opts.gene_cover[0]))
    visited = Set()
    for t in gene_tree.traverse('preorder'):
        if t.id not in visited:
            visited.add(t.id)
            if len(t.centroid) == 1:
                if t.centroid[0]['abundance'] >= abun_thres:
                    gs = list()
                    if t.name in gene_cover:
                        gs.extend([(gene_cover[t.name],t.name)])
                    for n in t.get_descendants():
                        if n.name in gene_cover:
                            gs.extend([(gene_cover[n.name],n.name)])
                        visited.add(n.id)
                    gs = sorted(gs,reverse=True)
                    if gs[0][0] >= opts.gene_cover[0]:
                        out = {'gene':gs[0][1], \
                               'depth':gene_abun[gs[0][1]], \
                               'abundance':t.centroid[0]['abundance'], \
                               'ratio':gs[0][0]}
                        if opts.taxonomy is not None:
                            out['taxonomy'] = gene_tax[gs[0][1]]
                        else:
                            out['taxonomy'] = ''
                        print '%(gene)s\t%(abundance)f\t%(depth)f\t%(ratio)f\t%(taxonomy)s' % out


if __name__ == "__main__":
    try:
        # TODO set logging
        logging.basicConfig(format="[%(asctime)s] %(levelname)s : %(message)s", level=logging.INFO)

        # TODO set command-line parser
        parser = argparse.ArgumentParser(description=globals()['__doc__'],epilog=textwrap.dedent('''\
                                         Examples\n\
                                         --------\n\
                                         1. find_seed_otus 99_otus_unannotated.tree gg_99_otus_abun.txt\n\
                                         \n'''),formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument('tree',help='phylogeny tree',nargs=1,metavar="TREE")
        parser.add_argument('abun',help='abundance profile',nargs=1,metavar="ABUNDANCE")
        parser.add_argument('-s',help='sequence similarity threshold (default: 0.9)',nargs=1,\
                            default=[0.9],type=float,dest='sim_thres',metavar="SIMILARITY")
        parser.add_argument('-d',help='depth cutting threshold (default:10)',nargs=1,\
                            default=[10],type=float,dest='depth_thres',metavar='DEPTH')
        parser.add_argument('-c',help='gene coverage threshold (default:0.6)',nargs=1,\
                            default=[0.6],type=float,dest='gene_cover',metavar='COVERAGE')
        parser.add_argument('-r',help='use adaptive depth threshold (0<=r<=1 of total depth)',\
                            type=float,dest='depth_ratio',metavar='RATIO')
        parser.add_argument('-T',help='gene taxonomy annotation',default=None,type=str,\
                            dest='taxonomy',metavar="TAXONOMY")
        parser.add_argument('-v',help='verbose optput',action="store_true",default=False,dest='verbose')
        opts = parser.parse_args()
        
        # TODO main code
        t_start = time.time()
        find_seed_otus()   

        # TODO complete
        if opts.verbose:
            logging.info('total elapsed time is %f minutes' % ((time.time()-t_start)/60.))
        sys.exit(0)
    except KeyboardInterrupt,e:
        raise e
    except SystemExit,e:
        raise e
    except Exception,e:
        logging.exception('ERROR, UNEXCEPTED EXCEPTION')
        logging.exception(e)
        traceback.print_exc()
        os._exit(1)
