#!/usr/bin/env python
"""
This python script is intended to be called as part of a job array. Use risky.py to calculate risk score.

AUTHORS
-------
    Daniel Rice
    Mari Niemi 
"""

import os
import glob
import numpy as np
import re
import argparse
import time
import gzip
import random
from pprint import pprint

# RE_CHR = re.compile(r'/chr(\d+)/')

def get_weights_dict(weights_file, weights_p_value_max):
    weights_dict = {}
    dup_count = 0
    duplicates_out = weights_file +'.duplicates'

    # Multiple access to the weight file by all of the chunks can cause access errors so 
    # attempt and sleep upto 20 times before giving up.
    max_attempts = 100
    for attempt in range(max_attempts):
        try:
            with open(weights_file) as f:
                weights_file_contents = f.readlines()
            break
        except IOError:
            time.sleep(random.randint(10, 600))
    else:
        raise Exception('Attempted to read {} times but failed each time so exiting.'.format(attempt + 1))

    duplicates = []
    for i, line in enumerate(weights_file_contents[1:]):
        line_split = [el.strip() for el in line.split()]
        chrom, pos, effect_allele, alt_allele, reference_allele, beta, pvalue, rsid = line_split
        pvalue = float(pvalue)
        beta = float(beta)

        if pvalue > weights_p_value_max:
            continue

        assert effect_allele == reference_allele

        k = (chrom, pos)
        if k in weights_dict:
            duplicates.append(line)
            weights_dict.pop(k, None)
        else: 
            weights_dict[k] = {
                    'effect_allele' : effect_allele,
                    'beta' : beta,
                    'rsid' : rsid,
                    'pvalue' : pvalue,
                    'alt_allele' : alt_allele,
                    }
        if not i % 1e6:
            print 'Processed {} entries from weights-file. Latest variant is: {}'.format(i, k)
            print weights_dict[k]

    print 'Total variants processed from weights_file: ', i
    print 'Number of chr:pos duplicates found in weights-file: ', len(duplicates)
    print 'Here are the duplicates:\n'
    
    for el in duplicates:
        print(el)
    
    return weights_dict


def calculate_risk_score(chunk_path, weights_dict, n_samples):
   
    # GWAS3 files
    RE_CHR = re.compile(r'/(\d+).gen.gz')
    
    # Mari's files
    # RE_CHR = re.compile(r'/chr(\d+)/') 
    
    complement_dict = {
        'A' : 'T',
        'T' : 'A',
        'G' : 'C',
        'C' : 'G',
        }

    count_dict = {
        'Total variants matched to weights file by chr, pos' : 0,
        'Alleles_discordant': 0,
        'Alleles_concordant': 0,
        'Alleles_swapped': 0,
        'Indels discarded' : 0,
        }

    weights_array = np.tile([2,1,0], n_samples)
    weights_array_flipped = np.tile([0,1,2], n_samples)
    total_risk = np.zeros(n_samples*3)

    print(chunk_path)
    
    m = RE_CHR.search(chunk_path)
    g = m.groups()
    chromosome = g[0]
    assert chromosome in [str(c) for c in range(1, 24)]
    snps_found = []
    
    _, file_extension = os.path.splitext(chunk_path)
    if file_extension.lower() == '.gz':
        f = gzip.open(chunk_path)
    else:
        f = open(chunk_path)       

    for line in f:
        # print '\nline\n', line
        # I don't think we need the strip here too, split should remove all whitespace
        # It can't hurt to keep it for the time being.
        line_split = line.strip().split()
        pos, ref_allele, alt_allele = line_split[2:5]
        if (len(ref_allele) > 1) or \
           (len(alt_allele) > 1) or \
           (ref_allele not in complement_dict.keys()) or \
           (alt_allele not in complement_dict.keys()):
            count_dict['Indels discarded'] += 1
            continue
        k = (chromosome, pos)
        print k
        if k in weights_dict:
            v = weights_dict[k]
            pprint(v)
            genotypes = np.array(line_split[5:], dtype= float)
            count_dict['Total variants matched to weights file by chr, pos'] += 1

            c_ref_allele = complement_dict[ref_allele]
            c_alt_allele = complement_dict[alt_allele]

            if ((ref_allele == v['effect_allele']) and (alt_allele == v['alt_allele'])) or \
               ((c_ref_allele == v['effect_allele']) and (c_alt_allele == v['alt_allele'])):
                # 2, 1, 0
                snps_found.append(k)
                total_risk += genotypes * weights_array * v['beta']
                count_dict['Alleles_concordant'] += 1

            elif ((ref_allele == v['alt_allele']) and (alt_allele == v['effect_allele'])) or \
                 ((c_ref_allele == v['alt_allele']) and (c_alt_allele == v['effect_allele'])):
                # 0, 1, 2
                snps_found.append(k)
                total_risk += genotypes * weights_array_flipped * v['beta'] 
                count_dict['Alleles_swapped'] += 1

            else:
                count_dict['Alleles_discordant'] += 1

    f.close()
    total_risk = total_risk.reshape((n_samples, 3)).sum(axis=1)

    return total_risk, count_dict, snps_found


def get_chunks(root_path):
    chunks = glob.glob(os.path.join(root_path, 'chr*/coreex.*.imputed'))
    return chunks


def get_n_samples(fam_path):
    with open(fam_path) as s:
        n_samples = sum(1 for line in s)
    print 'Number of samples in fam: ', n_samples
    return n_samples


def combine_results(results):

    total_risk = 0
    count_dict = Counter()
    snps_found = []

    for result in results:
        chunk_total_risk, chunk_count_dict, chunk_snps_found = result

        total_risk += chunk_total_risk        
        count_dict += chunk_count_dict
        snps_found += chunk_snps_found     
        
    return total_risk, count_dict, snps_found


def save_results(results_dir, chunk, total_risk, count_dict, snps_found):
    
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    f_total_risk = os.path.join(results_dir, 'total_risk')
    f_count_dict = os.path.join(results_dir, 'counts')
    f_snps_found = os.path.join(results_dir, 'snps_found')

    
    # Print total_risk to file
    with open(f_total_risk, 'a') as f:
        f.write('{}\t{}\n'.format(chunk, '\t'.join([str(el) for el in total_risk])))

    # Print count_dict to file
    keys = sorted(count_dict)
    line = '{}\t{}\n'.format(chunk, '\t'.join([str(count_dict[k]) for k in keys]))
    if not os.path.exists(f_count_dict):
        with open(f_count_dict, 'w') as f:
            header = 'chunk\t{}\n'.format('\t'.join(keys))
            f.write(header)
            f.write(line)
    else:
        with open(f_count_dict, 'a') as f:
            f.write(line)

    # Print snps_found to file
    with open(f_snps_found, 'a') as f:
        f.write('{}\t{}\n'.format(chunk, '\t'.join([str(el) for el in snps_found])))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--weights_path')
    parser.add_argument('--results_dir')
    parser.add_argument('--n_samples', type=int)
    parser.add_argument('--weights_p_value_max', type=float)
    parser.add_argument('--chunk')
    
    args = parser.parse_args()
   
    print('\nArguments:')
    print(vars(args))

    weights_path = args.weights_path
    results_dir = args.results_dir
    n_samples = args.n_samples
    weights_p_value_max = args.weights_p_value_max
    chunk = args.chunk
    
    weights_dict = get_weights_dict(weights_path, weights_p_value_max)
    total_risk, count_dict, snps_found = calculate_risk_score(chunk, weights_dict, n_samples)
    save_results(results_dir, chunk, total_risk, count_dict, snps_found)


if __name__== '__main__': 
    main()

    





