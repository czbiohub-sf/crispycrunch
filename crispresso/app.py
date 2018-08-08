"""
Minimal Flask app that does two things:
1) downloads fastq files from s3 bucket
2) executes crispresso
"""

import os
import shutil
import sys

from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count

from flask import Flask, jsonify, request

import s3
import seqs

from CRISPResso.CRISPRessoCORE import main as CRISPResso_main  # noqa

OUTPUT_DIR = 'output'

app = Flask(__name__, static_folder=OUTPUT_DIR)


@app.route('/')
def hello_world():
    return 'HEllo world'


# TODO (gdingle): make it POST only
@app.route('/crispresso', methods=['POST'])
def crispresso():
    """
    `s3_bucket` and `s3_prefix` make up the location of all of the fastq files
    produced by one illumina reading of a sample plate.

    `selected_guides` are expected to be a JSON map of target region to gRNAs.
    {
      "chr1:11130541-11130751": {
        "s180+": "TTTCTACTACTACAGTTGAG TGG",
        "s185+": "ACTACTACAGTTGAGTGGTC TGG",
        "s76+": "AAAAGAGTTCCAGAGTGCTC TGG"
      },
      "chr7:5569177-5569415": {
        "s191+": "TTCCGGCGCGCCGAGTCCTT AGG",
        "s207-": "TCCCCAATCTGGGCGCGCGC CGG",
        "s225+": "CGCCGGCGCGCGCCCAGATT GGG"
      }
    }

    Likewise, `selected_donors` are a map of target region to ssDNAs.
    """
    posted_data = request.get_json()

    fastqs = s3.download_fastqs(
        posted_data['s3_bucket'],
        posted_data['s3_prefix'],
        # overwrite=not posted_data.get('dryrun'))
        # TODO (gdingle):
        overwrite=False)

    # TODO (gdingle): should we do this transform on the crispycrunch side?
    # TODO (gdingle): cache get_reference_amplicon somehow... unfortunate lru_cache not python2
    amplicon_seqs = [
        (seqs.get_reference_amplicon(chr_loc), guide_seq.split(' ')[0])
        for chr_loc, guide_seqs in posted_data['selected_guides'].items()
        for guide_seq in guide_seqs.values()
    ]
    # how to get hdr amplicons here?
    # donor_guides = posted_data['donor_guides']

    # Although threads would be more efficient, CRISPResso is not thead-safe.
    # # TODO (gdingle): what is the optimal number of processes? CRISPResso
    # has its own pool internally. See n_processes.
    with ProcessPoolExecutor(max_workers=cpu_count()) as pool:
        futures = []
        for i, amplicon_seq in enumerate(amplicon_seqs):
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
                    amplicon_seq[0],
                    amplicon_seq[1],
                    None,  # TODO (gdingle): donor_guides
                    posted_data.get('dryrun'),
                ))

    urls = [f.result() for f in futures]
    return jsonify(urls)


def _analyze_fastq_pair(fwd,
                        rev,
                        amplicon_seq,
                        guide_seq,
                        expected_hdr_amplicon_seq,
                        dryrun):

    # Eliminate illumina boilerplate. See https://goo.gl/fvPLMa.
    results_name = fwd.split('/')[-1].split('_')[0]

    crispresso_args = [
        '/opt/conda/bin/CRISPResso',
        # TODO (gdingle): group by fwd and reverse
        '--fastq_r1', fwd,
        '--fastq_r2', rev,
        '--amplicon_seq', amplicon_seq,
        '--guide_seq', guide_seq,
        '--expected_hdr_amplicon_seq', expected_hdr_amplicon_seq,
        '-o', OUTPUT_DIR,
        '--save_also_png',
        '--trim_sequences',
        # TODO (gdingle): what is fasta adapter?
        '--trimmomatic_options_string', _get_trim_opt(),
        '--n_processes', str(cpu_count()),
        '--name', results_name,
    ]

    if not dryrun:
        _import_and_execute(crispresso_args)

    # Remove dir prefix. See https://goo.gl/s7dzEK .
    crispresso_results_path = OUTPUT_DIR + '/CRISPResso_on_' + results_name
    results_path = crispresso_results_path.replace('CRISPResso_on_', '')
    if os.path.exists(crispresso_results_path):
        if os.path.exists(results_path):
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
    try:
        CRISPResso_main()
    except SystemExit as e:
        if e.code != 0:
            raise e


if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5000, debug=True, threaded=True)
