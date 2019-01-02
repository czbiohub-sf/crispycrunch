# crispresso-test.sh
# This is for testing on remote machine inside docker container.
# docker exec -it mystifying_lovelace bash
cd /CRISPResso
curl -O https://s3-us-west-2.amazonaws.com/jasonli-bucket/CrispyCrunch/mNGplate3_unsorted_A1_POLR1A-C_S1_R1_001.fastq.gz
curl -O https://s3-us-west-2.amazonaws.com/jasonli-bucket/CrispyCrunch/mNGplate3_unsorted_A1_POLR1A-C_S1_R2_001.fastq.gz
python CRISPResso.py -r1 mNGplate3_unsorted_A1_POLR1A-C_S1_R1_001.fastq.gz -r2 mNGplate3_unsorted_A1_POLR1A-C_S1_R2_001.fastq.gz -g CTCTGAGATAGCAGCTACCC -a TGTACTGTCACTTGGAACTGCCCTGGATTCCAGTCTCCTGGTCCTCTCATGCAGAAGGCAGGCTGGGCCACGCCCTCACCAAGGGTCCTTGGAGCTGGGCAGATGGTGCCGGGGTAGCTGCTATCTCAGAGGCTGCTTGAGCTCGAACAGGCCTGTCCCGCCCCTGACGACCTTCCCGACCACAAGGCAGGCAGAAGGAGACCTCAGCTCATCGTGGGATCCTGACAGAGACACAAAAACATGTGTCAGGGTGTTGAG -e TGTACTGTCACTTGGAACTGCCCTGGATTCCAGTCTCCTGGTCCTCTCATGCAGAAGGCAGGCTGGGCCACGCCCTCACCAAGGGTCCTTGGAGCTGGGCAGATGGTGCCGGGGTAGCTCTTACATCATATCGGTAAAGGCCTTTTGCCACTCCTTGAAGTTGAGCTCGGTACCACTTCCTGGACCTTGAAACAAAACTTCCAATCCGCCACCTCTCAGAGGCTGCTTGAGCTCGAACAGGCCTGTCCCGCCCCTGACGACCTTCCCGACCACAAGGCAGGCAGAAGGAGACCTCAGCTCATCGTGGGATCCTGACAGAGACACAAAAACATGTGTCAGGGTGTTGAG \
  --debug \
  --suppress_report \
  --trim_sequences \
  # Number of seeds to test whether read is forward or reverse',default=5
  --aln_seed_count=2 \
  # Length of seeds to test whether read is forward or reverse',default=10
  --aln_seed_len=10 \
  # of seeds that must match to call the read forward/reverse',default=2
  --aln_seed_min=2

# cProfile.run("""main()""", sort='tottime')
#


