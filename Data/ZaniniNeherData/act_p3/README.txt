The data is stored in the following format. Each file is in Tab Separated
Value (TSV) format. Each file contains a separate time point. The estimated
number of days from infection is in the filename and again in the first line
of the file.

Each row is a position in the genome, in nucleotides. The coordinates refer to
the reference sequence for each patient, which can be downloaded in FASTA or
Genbank (GB) format from the website.

Each column is a possible nucleotide, according to the alphabet ACGT-N. The -
(dash) indicates deletions (with respect to the patient reference sequence),
N an ambiguous nucleotide. The latter indicates ambiguous nucleotides, should
be rare and can be ignored for most purposes.

Insertions are not included.
