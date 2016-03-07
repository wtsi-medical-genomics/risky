#!/usr/bin/env python

"""
Combines all of the risk computation output form risky_chunk.py


AUTHORS
-------
    Daniel Rice
    Mari Niemi 
"""


import os
import numpy as np
import argparse
from collections import defaultdict


def split_strip(line, sep=None):
    if sep:
        return [el.strip() for el in line.split(sep)]
    else:
        return [el.strip() for el in line.split()]



def combine_results(results_dir, sample_file, n_chunks):
    
    f_total_risk = os.path.join(results_dir, 'total_risk')
    f_count_dict = os.path.join(results_dir, 'counts')
    f_snps_found = os.path.join(results_dir, 'snps_found')

    
    f_total_risk_combined = os.path.join(results_dir, 'combined_total_risk')
    f_count_dict_combined = os.path.join(results_dir, 'combined_counts')
    f_snps_found_combined = os.path.join(results_dir, 'combined_snps_found')


    # Calculate total risk and output to f_total_risk_combined
    chunks = []
    total_risk = []
    with open(f_total_risk) as f:
        for line in f:
            tokens = split_strip(line, sep='\t')
            chunks.append(tokens[0])
            total_risk.append(tokens[1:])

    assert len(chunks) == n_chunks

    total_risk = np.array(total_risk, dtype=float)
    print(total_risk)
    total_risk = np.sum(total_risk, axis=0)

    print(total_risk)

    if sample_file:
        sample_names = []
        with open(sample_file) as f:
            _, file_extension = os.path.splitext(sample_file)
            if file_extension.lower() == '.sample':
                f.readline()
                f.readline()
            for line in f:
                tokens = split_strip(line)
                sample_names.append(tokens[0])

        print(sample_names)
        print(total_risk)

        with open(f_total_risk_combined, 'w') as f:
            for sample_name, sample_risk in zip(sample_names, total_risk):
                f.write('{}\t{}\n'.format(sample_name, sample_risk))

    else:
        with open(f_total_risk_combined, 'w') as f:
            for sample_risk in total_risk:
                f.write('{}\n'.format(sample_risk))

    # Combine counts dict
    with open(f_count_dict) as f:
        header = f.readline()
        tokens = split_strip(header, '\t')
        keys = tokens[1:]
        this_n_chunks = 0
        count_dict = defaultdict(int)
        for line in f:
            tokens = line.split('\t')
            for k, v in zip(keys, tokens[1:]):
                count_dict[k] += int(v)
            this_n_chunks += 1

    assert this_n_chunks == n_chunks

    with open(f_count_dict_combined, 'w') as f:
        for k, v in count_dict.iteritems():
            f.write('{}\t{}\n'.format(k, v))

    # Combine SNPs found
    snps_found = []
    this_n_chunks = 0
    with open(f_snps_found) as f:
        for line in f:
            tokens = split_strip(line, '\t')
            snps_found.extend(tokens[1:])
            this_n_chunks += 1

    assert this_n_chunks == n_chunks
    print(snps_found)
    with open(f_snps_found_combined, 'w') as f:
        for snp_found in sorted(set(snps_found)):
            f.write(snp_found + '\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--results_dir')
    parser.add_argument('--sample_file', default=None)
    parser.add_argument('--n_chunks', type=int)
    
    args = vars(parser.parse_args())
    print(args)
    
    combine_results(**args)
    

if __name__ == '__main__':
    main()
