"""
Minimal Flask app that does two things:
1) downloads fastq files from s3 bucket
2) executes crispresso
"""

import multiprocessing
import os
import shutil
import sys

from concurrent.futures import ProcessPoolExecutor

from flask import Flask, request, jsonify

import s3

OUTPUT_DIR = 'output'

app = Flask(__name__, static_folder=OUTPUT_DIR)


@app.route('/')
def hello_world():
    return 'HEllo world'


@app.route('/crispresso')
def crispresso():

    fastqs = s3.download_fastqs(
        request.args['s3_bucket'],
        request.args['s3_prefix'],
        overwrite=not request.args.get('dryrun'))

    with ProcessPoolExecutor() as pool:
        futures = []
        for i in range(0, len(fastqs), 2):
            fwd = fastqs[i]
            rev = fastqs[i + 1]
            assert 'R1' in fwd and 'R2' in rev, 'Fastq files must be paired and sorted'
            futures.append(pool.submit(_analyze_fastq_pair, fwd, rev))

    urls = [f.result() for f in futures]
    return jsonify(urls)


def _analyze_fastq_pair(fwd, rev):

    # Eliminate illumina boilerplate. See https://goo.gl/fvPLMa.
    results_name = fwd.split('/')[-1].split('_')[0]

    # TODO (gdingle): make it POST only
    crispresso_args = [
        '/opt/conda/bin/CRISPResso',
        # TODO (gdingle): group by fwd and reverse
        '--fastq_r1', fwd,
        '--fastq_r2', rev,
        '--amplicon_seq', request.args['amplicon_seq'],
        '--guide_seq', request.args['guide_seq'],
        '--expected_hdr_amplicon_seq', request.args.get('expected_hdr_amplicon_seq'),
        '-o', OUTPUT_DIR,
        '--save_also_png',
        '--trim_sequences',
        # TODO (gdingle): what is fasta adapter?
        '--trimmomatic_options_string', _get_trim_opt(),
        '--n_processes', str(multiprocessing.cpu_count()),
        '--name', results_name,
    ]

    if not request.args.get('dryrun'):
        _import_and_execute(crispresso_args)

    # Remove dir prefix. See https://goo.gl/s7dzEK .
    crispresso_results_path = OUTPUT_DIR + '/CRISPResso_on_' + results_name
    results_path = crispresso_results_path.replace('CRISPResso_on_', '')
    if os.path.exists(crispresso_results_path):
        shutil.rmtree(results_path)
        os.rename(crispresso_results_path, results_path)

    return results_path


def _get_trim_opt(adapterloc='fastqs/TruSeq3-PE-2.fa'):
    # TODO (gdingle): where to put adapterloc permanently? does it ever change?
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
    from CRISPResso.CRISPRessoCORE import main  # noqa
    try:
        main()
    except SystemExit as e:
        if e.code != 0:
            raise e


if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5000, debug=True, threaded=True)
