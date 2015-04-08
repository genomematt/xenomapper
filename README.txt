xenomappability
===============

xenomappability --fasta tests/data/test_from_EcoliK12DH10B.fasta --readlength 10 > tests/data/test_from_EcoliK12DH10B_10reads.fasta

bowtie2-build tests/data/test_from_EcoliK12DH10B.fasta tests/data/test_from_EcoliK12DH10B
bowtie2 -x tests/data/test_from_EcoliK12DH10B -f -U tests/data/test_from_EcoliK12DH10B_10reads.fasta -S tests/data/test_from_EcoliK12DH10B_10reads.sam

xenomappability --mapped_test_data tests/data/test_from_EcoliK12DH10B_10reads.sam > tests/data/test_from_EcoliK12DH10B_10reads.wig
xenomappability --single_end_wiggle tests/data/test_from_EcoliK12DH10B_10reads.wig --sam_for_sizes tests/data/paired_end_testdata_human.sam
