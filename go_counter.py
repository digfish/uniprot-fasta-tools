#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      sam
#
# Created:     06-04-2015
# Copyright:   (c) sam 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import os,sys
import urllib2
from Bio import SeqIO,ExPASy,SwissProt
import pickle
import operator



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

def dump_cs(dict):
    for key in dict.keys():
        record = dict [ key]
        print key, record.cross_references


def get_go_annot(dict,accession):
	cs = dict[accession].cross_references
	go_annot = {}
	for c in cs:
		if c[0] == 'GO':
			go_annot[ c[1] ] = c
	return go_annot

def list_go_acessions(uniprot_info):
    accession_and_gos = {}
    accessions = uniprot_info.keys()

    for an_accession in accessions:
    	accession_and_gos[ an_accession ] = get_go_annot(uniprot_info,an_accession)
    return accession_and_gos

def go_counter(_dict):
	go_counting = {}
	for accession in _dict.keys():
		go_record = _dict[accession]
		for go_id in go_record.keys() :
			cs = go_record[go_id][2]
			if cs in go_counting:
				go_counting [ cs ] += 1
			else:
				go_counting [ cs ] = 1
	return go_counting

def sort_dict(_dict):
	sorted_x = sorted(_dict.items(), key=operator.itemgetter(1))
	sorted_x.reverse()
	return sorted_x

def go_evidence_codes_counter(_dict):
	go_counting = {}
	for accession in _dict.keys():
		go_record = _dict[accession]
		for go_id in go_record.keys() :
			cs = go_record[go_id][3]
			if cs in go_counting:
				go_counting [ cs ] += 1
			else:
				go_counting [ cs ] = 1
	return go_counting


def main():
    #os.chdir('D:\Dropbox\code\python\HW0002')


    if ( len( sys.argv ) < 2 ):
        print "No fasta file with results provided!"
        exit()

    fasta_blast_results_file = sys.argv[1]
#    fasta_blast_results_file = "HERG.ncbi-psiblast.500.swsissprot.fasta.txt"

    uniprot_info = {}

    if os.path.isfile("swissprot.annot.pickled"):
        f = open("swissprot.annot.pickled","rb")
        uniprot_info = pickle.load(f)
        f.close()
    else:
        for seq_record in SeqIO.parse ( fasta_blast_results_file, "fasta"):
            accession = extract_accession ( seq_record )
            print "Querying uniProt about %s" % (accession)
            get_SwissProt(uniprot_info,accession)
        pickle.dump(uniprot_info,open('swissprot.annot.pickled','wb'))

    #dump_cs(uniprot_info)
##    for accession in uniprot_info.keys():
##        print get_go_annot(uniprot_info,accession)

    #print go_counter(uniprot_info)

    go_annot = list_go_acessions(uniprot_info)

    go_counts_dict = go_counter(go_annot)

    sorted_go_counts_dict = sort_dict(go_counts_dict)

    import csv

    writer = csv.writer(open('go_stats.csv','wb'), delimiter = '\t', quotechar = '"')

    writer.writerow(("Annotation","Ocurrences"))
    for row in sorted_go_counts_dict:
        writer.writerow(row)

    go_evidence_dict = go_evidence_codes_counter(go_annot)

    sorted_go_evidence_dict = sort_dict( go_evidence_dict )

    writer = csv.writer(open('evidence_go_stats.csv','wb'), delimiter = '\t', quotechar = '"')

    writer.writerow(("Evidence","Ocurrences"))
    for row in sorted_go_evidence_dict:
        writer.writerow(row)


    #print sorted_go_evidence_dict


if __name__ == '__main__':
    main()
