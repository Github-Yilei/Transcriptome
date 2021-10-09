#/usr/bin/python

#-*-coding: utf-8 -*-

import argparse
import sys
from pathlib import Path

# building a gene id dictionary
def gff2gene(gff3):
    gene_dict = dict()
    with open(gff3, 'r') as gff3:
        for line in gff3:
            if 'gene' in line:
                info = line.split()[8]
                gene = info.split('=')[2]
                gene_dict[gene] = [gene]

    return gene_dict

# parsing expression matrix files
def path2file(input_path):
    p = Path(input_path)
    file_list = list(p.glob("*/*gene_abundances.tsv"))

    return file_list

def TidyMatrix(args):
    gene_dict = gff2gene(args.gff3)
    file_list = path2file(args.input_path)
    matrix_header = ["gene_id"]
    result_dict = {}

    for f in file_list:
        sample = f.name.replace('_gene_abundances.tsv', '')
        matrix_header.append(sample)

        with open(f, 'r') as tsv:
            next(tsv)
            for line in tsv:

                line_spl = line.split()
                gene_id = line_spl[0]
                FPKM = line_spl[7]
                gene_dict[gene_id].append(FPKM)



    # saving
    out_put = args.input_path + '/tidymartix.tsv'
    tidymatrix = open(out_put, 'w')
    tidymatrix.write("\t".join(matrix_header) + "\n")

    for k,y in gene_dict.items():
       tidymatrix = open(out_put, 'a')
       tidymatrix.write("\t".join(y) + "\n")

if __name__=='__main__':

    parser = argparse.ArgumentParser(description =  "Parsing the josn files of fastp and building a ploting-ready data frame for R")
    parser.add_argument('--input',
                            dest = 'input_path',
                            help = 'The absolute path of the default results produced by StringTie')
    parser.add_argument('--gff3',
                            dest = 'gff3',
                            help = 'The gff3 file')

    if len(sys.argv) <= 1:
        parser.print_help()

        sys.exit()
    args = parser.parse_args()

    TidyMatrix(args)
