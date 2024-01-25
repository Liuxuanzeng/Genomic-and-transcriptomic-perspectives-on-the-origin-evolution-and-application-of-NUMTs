# Genomic and transcriptomic perspectives on the origin, evolution and application of NUMTs sequences in Orthoptera 
The major codes in the analysis of Nuclear mitochondrial DNA fragments (NUMTs).

NUMT_blast_allmaker.py

To extract NUMT fragments, flanking 25bp, flanking 100bp, and flanking 2000bp sequences and calculate the GC content. Because mtDNA is a circular molecule, this code merged the NUMT fragments in the BLASTn result due to linear alignment.NUMT fragments, flanking 25bp, flanking 100bp, and flanking 2000bp sequences.

usage： python NUMT_blast_allmaker.py [genome.fasta] [mt.fasta] [mt.gff] [species name]



same_NUMT.py

Analyze NUMTs replication events after NUMTs insertion. Look for identical NUMTs and their flankers within species.The output result is same_NUMT_link.txt.

usage： same_NUMT.py	[species name]

NUMT_map_to_MT_isoseq_trans.py

Find mitochondrial transcripts consistent with NUMTs based on Pacbio iso-seq CCS data. And extract these matched mitochondrial transcripts in fasta format.

usage： python3 NUMT_map_to_MT_isoseq_trans.py [minimap.out] [spname] [0.8] [20] [isoseq.fa]

[0.8] Filter alignment results with sequence similarity below 80%. [float] can be modified as needed.

[20] The length difference between mitochondrial transcripts and NUMTs is less than 20bp. [int] can be modified as needed.
