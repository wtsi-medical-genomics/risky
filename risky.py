#!/usr/bin/env python

import gzip
import os
import glob
import sys
import shutil
import subprocess
import re
import pprint
import yaml
import argparse

"""
This python script will:
    1. Submit a job array over all chunks.
    2. Submit a final job that will run when the chunks job array finshes.

The parameters should be specified in the p dictionary below, nothing else should need to be changed.


PARAMS
------

    chunks
        A directory of the chunks to be computed. Use the wildcard character * to grab all of your files.
    
    chr-regex
        The regular expression to find the chromosome.

    memory
        How many Mb you want to use for the job array risk calculation.
    
    weights
        The path to the weights file.
    
    working-dir
        Where all of the params, scripts, o, e, and results are saved. See OUTPUT below.
    
    sample-file
        The path to a sample file.

    weights-p-value max
        The maximum p-value of a weight to include in the risk calculation.

    farm-group (OPTIONAL)
        The farm-group to run the jobs against.

    queue
        Which farm queue to send the jobs to.

OUTPUT
------

    working dir/
        name/
            params
            o/
            e/
            scripts/
            results/


AUTHORS
-------
    Daniel Rice
    Mari Niemi 

"""

def get_choice():
    yes = set(['yes','y', 'ye', ''])
    no = set(['no','n'])
    choice = raw_input('([y]/n)? ').lower()
    if choice in yes:
       return True
    else:
       return False


def get_script_path():
    return os.path.realpath(__file__)


def get_script_directory():
    return os.path.dirname(os.path.realpath(__file__))


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('params')
    args = parser.parse_args()

    with open(args.params) as f:
        p = yaml.load(f)
    
    observed = set(p.keys())
    required = set(['working-dir', 'weights', 'memory', 'chunks', 'weights-p-value-max'])
    short = required - observed
    if short:
        print('Error, the params file does not contain all of the paramters! Missing the following:\n')
        for el in short:
            print('\t{}'.format(el))
        sys.exit()

    print('\n\nUsing the following parameters from file {}:\n'.format(args.params))
    pprint.pprint(p)

    p['job-name'] = os.path.basename(p['working-dir'])

    if os.path.exists(p['working-dir']):
        print('\n\n{working-dir} already exists. Do you want me to delete it and proceed?'.format(**p))
        choice = get_choice()
        if not choice:
            sys.exit()
        else:
            shutil.rmtree(p['working-dir'])

    os.makedirs(p['working-dir'])

    params_file = os.path.join(p['working-dir'], 'params')
    with open(params_file, 'w') as f:
        f.write('Running python script:\n')
        f.write(get_script_path())
        f.write('\n\nWith paramsters:\n')
        pprint.pprint(p, f)


    for directory in ['o', 'e', 'scripts', 'results']:
        path = os.path.join(p['working-dir'], directory)
        os.makedirs(path)
        p[directory] = path
   
    chunks = glob.glob(p['chunks'])
    p['n-chunks'] = len(chunks)

    if 'sample-file' in p:
        with open(p['sample-file']) as f:
            n_samples = sum(1 for line in f)
            _, file_extension = os.path.splitext(p['sample-file'])
            if file_extension == '.sample':
                n_samples -= 2
            p['n-samples'] = n_samples 
            p['sample-file'] = '--sample_file={sample-file}'.format(**p)
    else:
        _, file_extension = os.path.splitext(chunks[0])
        if file_extension.lower() == '.gz':
            f = gzip.open(chunks[0])
        else:
            f = open(chunk_path)       
        header = f.readline()
        f.close()
        header = header.split()
        genotypes = header[5:]
        n_samples = len(genotypes)/3.
        assert n_samples == int(n_samples)
        p['n-samples'] = int(n_samples)
        p['sample-file'] = ''
    p['risky-chunk'] = os.path.join(get_script_directory(), 'risky_chunk.py')
    
    array_job_element = """
    #!/usr/bin/env bash

    CHUNKS=({chunks})

    {risky-chunk} \\
    --weights_path={weights} \\
    --results_dir={results} \\
    --n_samples={n-samples} \\
    --weights_p_value_max={weights-p-value-max} \\
    --chunk=${{CHUNKS[$LSB_JOBINDEX - 1]}}
    """.format(**p)

    p['array-job-element'] = os.path.join(p['scripts'], 'array_job_element.sh')
    with open (p['array-job-element'], 'w') as f:
        f.write(array_job_element)
    os.chmod(p['array-job-element'], 0700)


    if 'farm-group' in p:
        p['farm-group'] = '-G {farm-group}'.format(**p)
    else:
        p['farm-group'] = ''


    if 'queue' in p:
        p['queue'] = '-q {queue}'.format(**p)
    else:
        p['queue'] = ''
    
    run_array = """
    #!/usr/bin/env bash

    bsub \\
    -J"{job-name}[1-{n-chunks}]" \\
    {farm-group} \\
    {queue} \\
    -o {o}/chunk-%J-%I \\
    -e {e}/chunk-%J-%I \\
    -M{memory} \\
    -R"select[mem>{memory}] rusage[mem={memory}]" \\
    {array-job-element}
    """.format(**p)
    fpath = os.path.join(p['scripts'], 'run_array.sh')
    with open (fpath, 'w') as f:
        f.write(run_array)
    os.chmod(fpath, 0700)

    print '\nRunning chunks job array from script at: ' + fpath
    job_submit = subprocess.check_output(fpath, shell=True)
    print(job_submit)

    RE_JOBID = re.compile(r'Job \<(\d+)\>')
    m = RE_JOBID.search(job_submit)
    if not m:
        print('Could not find the Job ID for the chunks job array!')
        sys.exit()

    g = m.groups()
    p['job-id'] = g[0]
    p['risky-combine'] = os.path.join(get_script_directory(), 'risky_combine.py')

    combine = """
    #!/usr/bin/env bash

    bsub \\
    -w "done({job-id})" \\
    {farm-group} \\
    {queue} \\
    -o {o}/combine-%J \\
    -e {e}/combine-%J \\
    -M{memory} \\
    -R"select[mem>{memory}] rusage[mem={memory}]" \\
    {risky-combine} \\
    --results_dir={results} \\
    {sample-file} \\
    --n_chunks={n-chunks} \\
    """.format(**p)
    
    print(combine)

    fpath = os.path.join(p['scripts'], 'combine.sh')
    with open (fpath, 'w') as f:
        f.write(combine)
    os.chmod(fpath, 0700)

    print 'Queueing combine results job to run when {job-id} finishes'.format(**p)
    job_submit = subprocess.check_output(fpath, shell=True)
    print(job_submit)


if __name__ == '__main__':
    main()

