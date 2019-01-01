# Goals
This is what I want to present on Wednesday, January 2nd @7am:
* Coding
	* Wrangling_CIRM_Data: get downloadling shell script working for logging failures
	* Breaking_Ontogeny: git Clays' bash script working, describe missing files, add example data to play w/ and mini-tutorial
	* Review Introns, Exons, and Promoters
* Scientific Writing
	* ALT
		* Is TMM inhibition a universal cancer treatment?
		* Markers of ALT Activity
		* Promoter Compaction and Exon Deletion can Initiate ALT
		* ALT remodels telomere architecture
	* Stem Cells
		* MSC promoter and ATRX => ALT
		* HSC dysfunction may lead to AD (microglia)
		* iPSC seem to use ALT
	* Sequencing	
		* Review of DNA Sequencing methods from CIRM
		* Exon sequencing ALT
		* Sequencing Identifying Gene Signature in ALT
==========================================================================
Repo: ALT_Introns_Exons_and_Promoters
Description: I used Python and R to explore the Alternative Lengthening of Telomeres literature in humans, mice, c. elegans, and yeast.
Website: davehalvorsen.github.oi/ALT

# Introduction
This is a programmatic exploration of the Alternative Lengthening of Telomeres (ALT) literature in humans, mice, c. elegans, and yeast. 

# POT-1 Deficiency Creates ALT+ C. Elegans Strains
ALT can happen in C elegans! Most human cancers have long/heterogenous telomeres 
telomeres cap linear chromosomes cause DNA polymerases can't comleteely copy chromosomes
telomerase adds telomeric repeats w/ reverse trascription
shortening leads to senesnce AND is part of stop cancer
senscence => proliferation invovles teloemre loss and crisis => cell death and crhomsome isntability
10-15% cancers use ALTspontaneous ALT have long/heterogenous OR normal OR both. 
POT2 represses normal telomere length ALT. mamallian POT1 has homologs in C elegans
pot-1 CeOB2 and pot-2 CeOB1 
pot1 mutants reported to have wild telomere lenghts vs. pot2 normal lengths. 
authors created a variety of double mutants
authors found ALT C elegans that had long or normal length
FIG3D trt-1 & POT-2 absence lead to ALT Caenorhabditis elegans with expected telomere length
FIG 3D trt-1 & pot-1 mutants reported to have heterogenous telomeres like human ALT
CHENG 2012 CelegansTELOMEREandSURVIVAL.png
#### Reading and Writing pot-1 FASTA Files

#### Multiple Sequence Alignment of pot-1 Genes

#### Multiple Sequence Alignment of pot-1 Proteins

#### Identifying pot-1 Introns

#### Identifying pot-1 Exons

#### Displaying pot-1 Open Reading Frames

##### pot-1 Alternative Splicing

# TERT Promoter Compaction Causes ALT in Human and Mouse Cells
Lack of Telomerase Gene Expression in Alternative Lengthening of Telomere Cells Is Associated with chromatin Remodeling of the hTR and hTERT Gene Promoters
Reversible Conversion of Immortal Human Cells from Telomerase-Positive to Telomerase-Negative Cells
DNA methyltransferases control telomere length and telomere recombination in mammalian cells


### Identifying The hTERT Promoter Sequence
Cong 1999 The core promoter region (of hTERT) identified in this
study extends from –330 to +361 bp

#### Identifying the hTERT CpG Island Region
Kumakura 2005
the hTERT CpG island which is from 654 bp upstream of the putative
transcription start site to 510 bp downstream of the transcription start site
were the same as those described by Dessain et al. (8) 
8 is Dessain SK, Yu H, Reddel RR, Beijersbergen RL,
Weinberg RA. Methylation of the human telomerase
gene CpG island. Cancer Res 2000;60:537–41.

Dessain 2000
On the basis of the quantitative
criteria proposed by Antequera and Bird (18), this CpG island is from
654 bp upstream of the putative transcription start site (6) to 510 bp
downstream of the transcription start site, ending 56 bp after the start
of the first intron (6 – 8, 18). 5
Within this region, the DNA has a GC
content of 74% and a CG:GC ratio of 0.87
#### Obtaining hTERT WITH the CpG Island Region
IMAGE TERT REVERSE ARROW INDICATES REVERSE STRAND
https://www.ncbi.nlm.nih.gov/gene/7015
INITIAL SETUP
from:  
1253167
 to:  
1295047
Cong 1999 start codon of ATG seems to be on line 1 of the hTERT
"ATGCCGCGCGCT" is end of the line, which is 59 in from the left (59 is A of atg)
CpG is 654 bp upstream of transcriptoin start site, SO 595 back from current start
SO 1253167-595 = 1252572 should be start site now
I might be off by 1 or so, BUT I nkow that its 654 upstream of start to 510 bp downstream of start site (ending 56 bp after start of first intron). I can check the m ath on that location to be sure :)
I THINK I WAS BACKWARDS?!?
What about 1295047 + 595 = 1295642 YESSSSSSS, that's right :)
now i have more at the beginniig, so "atgccgcgcgctccccgct" is fully searchable!
https://www.ncbi.nlm.nih.gov/nuccore/NC_000005.10?report=fasta&from=1253167&to=1295642&strand=true


#### GC content of CpG Promoter R

#### Identifying Exons and Introns of hTERT

#### hTERT Alternative Splicing

# ATRX Exon Deletion is Required for Human ALT Activity
BLASTAlignRetrieve/ID mappingPeptide searchContactHelp
UniProtKB - Q9H668 (STN1_HUMAN)
PART OF CST COMPLEX AND ALT-INVOLVED!!!!!!
https://www.uniprot.org/uniprot/Q9H668
#### Alternative Splicing in ATRX


# STN1 Mutation Triggers ALT in Yeast
A Mutation in the STN1 Gene Triggers an Alternative Lengthening of Telomere-Like Runaway Recombinational Telomere Elongation and Rapid Deletion in Yeast
https://mcb.asm.org/content/25/18/8064

# Data Sources

# Citations
Cheng 2012 Caenorhabditis elegans POT-2 telomere protein represses a mode of alternative lengthening of telomeres with normal telomere lengths

Cong 1999 The human telomerase catalytic subunit hTERT: organization of the gene and characterization of the promoter














# Reading and Writing FASTA
* Print_hATRX_Gene_Fasta.py
* Print_hATRX_Protein_Fasta.py
* Create_New_ATRX_FASTA_Files.py

# Simplified (Incorrect) DNA => Protein
* Print_Standard_Codon_Table.py
* 

# Introns
a segment of a DNA or RNA molecule that does not code for proteins and interrupts the sequence of genes.


# Exons
a segment of a DNA or RNA molecule containing information coding for a protein or peptide sequence.


# Promoters
n genetics, a promoter is a region of DNA that initiates transcription of a particular gene. Promoters are located near the transcription start sites of genes, on the same strand and upstream on the DNA (towards the 3' region of the anti-sense strand).


# Open Reading Frames
In molecular genetics, an open reading frame (ORF) is the part of a reading frame that has the ability to be translated. An ORF is a continuous stretch of codons that begins with a start codon (usually AUG) and ends at a stop codon (usually UAA, UAG or UGA).
























========================
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

