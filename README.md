# Transcription and Translation Review (Python & R)
I recently mixed up Exons with Introns :S ... here's a review of Exons, Introns, Promoters, and Open Reading Frames with examples done in R Studio & Python:

#### Exons and Introns (Python)
These programs can be found in the Exons_and_Introns folder:
* DNA_to_RNA_Simple_String.py takes a simple DNA string as input and returns an RNA string.
* RNA_to_AA_Simple_String.py takes a simple RNA string as input and returns a protein string. It *does* take stop codons into account, but this version can only do one reading frame. 
* Splicing_Out_Introns_1.py returns a sequence of ONLY the coding regions of the DNA string.
* Splicing_Out_Introns_2.py returns a sequence of ONLY the coding regions of the DNA and it returns the percentage of coding DNA.
* Splicing_Out_Introns_3.py returns a sequence of coding bases as UPPER CASE and non-coding bases as lowercase.
* RNA_Splicing.py is a simple example of removing the introns from a small RNA string.
* NOTE: These programs assume only one reading frame and one DNA string. Computing all six reading frames is *way* beyond my skills! I SHOULD DO AN EXAMPLE WITH BIOPYTHON: https://www.biostars.org/p/97746/

![PLACEHOLDER](/Assets/rna-protein-dictionary.jpg "PLACEHOLDER")
![PLACEHOLDER](/Assets/stop_codons_trinucleotideCURRENT.jpg "PLACEHOLDER")

#### Promoter Regions (R Studio)
These files can be found in the Promoters folder:
* how_to_extract_promoters_positions.Rmd is an R tutorial that I followed on determinging promoter regions from TxDb.Hsapiens.UCSC.hg19.knownGene.
* TF_Binding_to_DNA_Promoter_Regions.Rmd is an R tutorial that I followed for searching out potential transcription factor / DNA promoter region interactions.  

#### Detecting Open Reading Frames (Python, R Studio)
These files can be found in the Open_Reading_Frames folder:
* identifying_open_reading_frames.py is a BioPython tutorial for finding ORFs.
* identifying_open_reading_frame_positions.py is a BioPython tutorial for finding the positions of ORFs. 
* ORFik_Overview.Rmd is an R Studio tutorial that I completed on ORFik, which is a package for exploring open reading frames. 
* FindingGenesWithORFs.py this is an incorrect solution to a Rosalind.info ORF problem ... trying to figure things out.

![PLACEHOLDER](/Assets/identifying-open-reading-frames.jpg "PLACEHOLDER")

