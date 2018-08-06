import os
import subprocess

from flask import Flask, request

app = Flask(__name__)


@app.route('/')
def hello_world():
    return 'HEllo world'


def _get_trim_opt(adapterloc):
    # adapterloc = args.input_dir + "/TruSeq3-PE-2.fa"
    leadingqual = 3
    trailingqual = 3
    slidewindsiz = 4
    slidewindqual = 20
    minlength = 50

    # TODO (gdingle): clean up
    return 'ILLUMINACLIP:' + adapterloc + ':2:30:10 LEADING:' + str(leadingqual) \
        + ' TRAILING:' + str(trailingqual) + ' SLIDINGWINDOW:' + str(slidewindsiz) + \
        ':' + str(slidewindqual) + ' MINLEN:' + str(minlength)


@app.route('/crispresso')
def crispresso():
    # TODO (gdingle): make it POST only
    crispresso_cmd = [
        '/opt/conda/bin/CRISPResso',
        '--fastq_r1', request.args['fastq_r1'],
        '--fastq_r2', request.args['fastq_r2'],
        '--amplicon_seq', request.args['amplicon_seq'],
        '--guide_seq', request.args['guide_seq'],
        '--expected_hdr_amplicon_seq', request.args['expected_hdr_amplicon_seq'],
        '-o', 'output',
        '--save_also_png',
        '--trim_sequences',
        # TODO (gdingle): how to verify this is working?
        '--trimmomatic_options_string',
        "'" + _get_trim_opt(request.args['adapterloc']) + "'",
    ]

    # TODO (gdingle): why doesn't subprocess.check_output work here?
    os.system(' '.join(crispresso_cmd))

    return ' '.join(crispresso_cmd)


if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5000, debug=True)
