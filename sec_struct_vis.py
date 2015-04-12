#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Samuel Viana on 2015-04-08.
Copyright (c) 2015 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from Bio import SeqIO,ExPASy,SwissProt
from collections import OrderedDict
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pickle
import urllib2

def extract_accession(seq_record):
    id = seq_record.id
    tokens = id.split("|")
    return tokens[3].split('.')[0]

def get_SwissProt(dict,accession):
    try:
        handle = ExPASy.get_sprot_raw(accession)
        record = SwissProt.read(handle)
        dict[accession] = record
    except urllib2.HTTPError, error:
        print accession + ": protein not found on UniProt . "
    except ValueError as ve:
        print ve
    except Exception as ex:
        print ex

def dump_descriptions(dict):
    for key in dict.keys():
        record = dict [ key]
        print key + ' => ' +  record.description

def dump_features(dict):
    for key in dict.keys():
        record = dict [ key]
        print key, record.features

def get_ss_features(dict,accession):
    if accession in dict.keys():
        features = dict[accession].features
        features_collection = {}
        for feat in features:
            if feat[0] in ('STRAND','TURN','HELIX'):
                features_collection[ feat[1] ] = feat
        return features_collection
    else: return {}

def generate_dict_acc_struc(fasta_blast_results_file,uniprot_info,dict_acc_struc):
    no_sec_struc_counter = 0
    total_records_counter = 0

    for seq_record in SeqIO.parse ( fasta_blast_results_file, "fasta"):
        total_records_counter += 1
        accession = extract_accession ( seq_record )
        #print "Querying uniProt about %s" % (accession)
        get_SwissProt(uniprot_info,accession)
#        dump_features(uniprot_info)
        ss_feats = get_ss_features(uniprot_info,accession)
#        print ss_feats
        if len(ss_feats) > 0:
            print "Great ! Uniprot stores the secondary structure for " + accession
            seq_positions_sorted = sorted(ss_feats.keys())
#            print seq_positions_sorted
            ss_feats_sorted = OrderedDict()
            for pos in seq_positions_sorted:
                #print str(pos) + " => ", ss_feats[pos]
                ss_feats_sorted[pos] = ss_feats[pos]
            dict_acc_struc[ accession ] = ss_feats_sorted

        else: no_sec_struc_counter  += 1
    print  no_sec_struc_counter, " proteins with no secondary structure on Unirprot counted from a total of ", total_records_counter, 100.0 * ( float(no_sec_struc_counter) / float(total_records_counter)), " %"

# print the protein segments in one tow
def chart_protein_data(accession,protein_data,line):

    ypos = 100 - line

    plt.text(5, ypos + .2, accession, ha='left')

    for segment_start in protein_data.keys():
        segment_data = protein_data [ segment_start ]
        start = segment_data [1]
        end = segment_data [2]
        _type = segment_data [0]
        _color = ''
        if   _type == 'HELIX' : _color = 'red'
        elif _type == 'STRAND': _color = 'green'
        elif _type == 'TURN':   _color = 'blue'

        plt.plot([start,end], [ypos, ypos], color = _color, linewidth = 8)

def dump_dict_acc_struc(dict_acc_struc):
    _file = open("acc_struct.txt","w")
    for accession in dict_acc_struc.keys():
        segments = dict_acc_struc[accession]
        segment_string = []
        for segment_start in segments.keys():
            segment_data = segments [ segment_start ]
            start = segment_data [1]
            end = segment_data [2]
            _type = segment_data [0]
            segment_string.append( _type + ":(" + str(start) + ":" + str(end) + ")" )

        segments_info_string = ','.join(segment_string)
        _file.write(accession + "=>" + segments_info_string + "\n")
    _file.close()

def print_graph(dict_acc_struc):
    line = 0
    for accession in  dict_acc_struc.keys():
        a_protein = dict_acc_struc[accession]
        chart_protein_data(accession,a_protein, line)
        line += 1

    frame = plt.gca()

    #frame.axes.get_xaxis().set_visible(False)
    #frame.axes.get_yaxis().set_visible(False)
    frame.axes.get_yaxis().set_ticks([])

    plt.axis([0,1500,101 - len(dict_acc_struc) - 1, 101])
    plt.ylabel('PROTEINAS')
    plt.xlabel('COMPRIMENTO DAS SEQUENCIAS (numero de residuos)')
    #blue_patch = mpatches.Patch( color='blue', label='TURN')
    #plt.legend( handles = [blue_patch])
    plt.show()


def main():
    #os.chdir('D:/dropbox/code/python/HW0002')
    #fasta_blast_results_file = "HERG.ncbi-psiblast.500.swsissprot.fasta.txt"

    if ( len( sys.argv ) < 2 ):
        print "No fasta file with results provided!"
        exit()
    else:
        fasta_blast_results_file = sys.argv[1]

    uniprot_info = {}
    dict_acc_struc = OrderedDict()

    if os.path.isfile("swissprot.sec_struc.pickled"):
        f = open("swissprot.sec_struc.pickled","rb")
        dict_acc_struc = pickle.load(f)
        f.close()
    else:
        generate_dict_acc_struc(fasta_blast_results_file, uniprot_info, dict_acc_struc )
        pickle.dump(dict_acc_struc,open('swissprot.sec_struc.pickled','wb'))

    print_graph (dict_acc_struc)

    #dump_dict_acc_struc ( dict_acc_struc )


#    for accession in dict_acc_struc.keys():
#        print accession, dict_acc_struc[ accession ]
#        print "--"

    #print dict_acc_struc

if __name__ == '__main__':
    main()

