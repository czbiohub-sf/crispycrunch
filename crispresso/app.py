"""
Minimal Flask app that does two things:
1) downloads fastq files from s3 bucket
2) executes crispresso
"""

import os
import shutil
import sys

from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count

from flask import Flask, jsonify, request

import pandas
import s3

from CRISPResso.CRISPRessoCORE import main as CRISPResso_main  # noqa

OUTPUT_DIR = 'output'

app = Flask(__name__, static_folder=OUTPUT_DIR)


@app.route('/')
def hello_world():
    return 'Hello world'


@app.route('/analyze', methods=['POST'])
def crispresso():
    """
    `s3_bucket` and `s3_prefix` make up the location of all of the fastq files
    produced by one illumina reading of a sample plate.
    """
    post_data = request.get_json()
    sheet = pandas.DataFrame(post_data['sheet'])

    # TODO (gdingle): make fastqs custom named such as
    # CC8763_chr3:1000-2000_chr3:1200-1220__S1_L001_R1_001.fastq.gz
    fastqs = s3.download_fastqs(post_data['s3_bucket'], post_data['s3_prefix'])

    # how to get hdr amplicons here?

    # TODO (gdingle): is this correct? what about rev direction?
    # assert all(pair[1] in pair[0] for pair in amplicon_seqs), \
    #     'The guide sequence must be present in the amplicon sequence'

    # Although threads would be more efficient, CRISPResso is not thead-safe.
    # # TODO (gdingle): what is the optimal number of processes? CRISPResso
    # has its own pool internally. See n_processes.
    pool = ProcessPoolExecutor(max_workers=cpu_count())
    futures = _start_all_analyses(pool, sheet, fastqs, post_data.get('dryrun'))

    # In async mode, return the paths where results are expected to be written.
    if post_data.get('async'):
        results = [_results_path(fastqs[i])
                   for i in range(0, len(sheet), 2)]
    else:
        results = [f.result() for f in as_completed(futures)]

    return jsonify({
        'fastqs': fastqs,
        'results': results,
    })


def _start_all_analyses(pool, sheet, fastqs, dryrun):
    futures = []
    for i, row in enumerate(sheet.to_records()):
        fwd = fastqs[i * 2]
        rev = fastqs[i * 2 + 1]
        assert'_R1_' in fwd and '_R2_' in rev, 'Fastq files must be paired and sorted'
        # TODO (gdingle): how to handle missing or extra fastq files?
        # TODO (gdingle): match fastq files and selected_guides on sample names
        futures.append(
            pool.submit(
                _analyze_fastq_pair,
                fwd,
                rev,
                row['target_seq'],
                row['guide_seq'],
                None,  # TODO (gdingle): donor_guides
                dryrun,
            ))
    return futures

# Eliminate illumina boilerplate. See https://goo.gl/fvPLMa.


def _results_name(fwd):
    return fwd.split('/')[-1].split('_')[0]


def _crispresso_results_path(fwd):
    return OUTPUT_DIR + '/CRISPResso_on_' + _results_name(fwd)


def _results_path(fwd):
    return _crispresso_results_path(fwd).replace('CRISPResso_on_', '')


def _analyze_fastq_pair(fwd,
                        rev,
                        amplicon_seq,
                        guide_seq,
                        expected_hdr_amplicon_seq,
                        dryrun):

    crispresso_args = [
        '/opt/conda/bin/CRISPResso',
        '--fastq_r1', fwd,
        '--fastq_r2', rev,
        '--amplicon_seq', amplicon_seq,
        '--guide_seq', guide_seq,
        '--expected_hdr_amplicon_seq', expected_hdr_amplicon_seq,
        '-o', OUTPUT_DIR,
        '--save_also_png',
        '--trim_sequences',
        '--trimmomatic_options_string', _get_trim_opt(),
        '--n_processes', str(cpu_count()),
        '--name', _results_name(fwd),
    ]

    if not dryrun:
        _import_and_execute(crispresso_args)

    # Remove dir prefix. See https://goo.gl/s7dzEK .
    results_path = _results_path(fwd)
    if os.path.exists(_crispresso_results_path(fwd)):
        if os.path.exists(results_path):
            shutil.rmtree(results_path)
        os.rename(_crispresso_results_path(fwd), results_path)

    success = os.path.exists(
        results_path + '/Quantification_of_editing_frequency.txt')

    return success, results_path


def _get_trim_opt(adapterloc='fastqs/TruSeq3-PE-2.fa'):
    # TODO (gdingle): where to put adapterloc permanently?
    # TODO (gdingle): these are all same as defaults... why?
    leadingqual = 3
    trailingqual = 3
    slidewindsiz = 4
    # TODO (gdingle): these are NOT same as defaults... why?
    slidewindqual = 20
    minlength = 50

    # TODO (gdingle): clean up
    return 'ILLUMINACLIP:' + adapterloc + ':2:30:10 LEADING:' + str(leadingqual) \
        + ' TRAILING:' + str(trailingqual) + ' SLIDINGWINDOW:' + str(slidewindsiz) + \
        ':' + str(slidewindqual) + ' MINLEN:' + str(minlength)


def _import_and_execute(crispresso_args):
    # HACK ALERT! Here we override sys.argv to run the CRISPResso script
    # in the same process as the server, so we can get better tracebacks
    # and more python goodness.
    sys.argv = crispresso_args
    try:
        CRISPResso_main()
    except SystemExit as e:
        if e.code != 0:
            raise e


if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5000, debug=True, threaded=True)
