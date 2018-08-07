"""
Minimal Flask app that does two things:
1) downloads fastq files from s3 bucket
2) executes crispresso
"""

import multiprocessing
import sys

from flask import Flask, request

import s3

app = Flask(__name__)


@app.route('/')
def hello_world():
    return 'HEllo world'


@app.route('/crispresso')
def crispresso():

    fastqs = s3.download_fastqs(
        request.args['s3_bucket'],
        request.args['s3_prefix'],
        overwrite=not request.args.get('dryrun'))

    # TODO (gdingle): make it POST only
    crispresso_args = [
        '/opt/conda/bin/CRISPResso',
        # TODO (gdingle): group by fwd and reverse
        '--fastq_r1', fastqs[0],
        '--fastq_r2', fastqs[1],
        '--amplicon_seq', request.args['amplicon_seq'],
        '--guide_seq', request.args['guide_seq'],
        '--expected_hdr_amplicon_seq', request.args['expected_hdr_amplicon_seq'],
        '-o', 'output',
        '--save_also_png',
        '--trim_sequences',
        # TODO (gdingle): what is fasta adapter?
        '--trimmomatic_options_string', _get_trim_opt(),
        '--n_processes', str(multiprocessing.cpu_count()),
    ]

    if not request.args.get('dryrun'):
        _import_and_execute(crispresso_args)

    # TODO (gdingle): return something useful
    return ' '.join(crispresso_args)


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
