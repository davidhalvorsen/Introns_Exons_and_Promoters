Repo: ALT_Introns_Exons_and_Promoters
Description: I used Python and R to explore the Alternative Lengthening of Telomeres literature in Homo sapiens, Caenorhabditis elegans, and yeast.
Website: davehalvorsen.github.oi/ALT

# Telomere Maintenance Mechanisms
Alternative Lengthening of Telomeres (ALT)
Telomerase (TEL)
Shay FIGURE DUAL INHIBITION
REDDEL EST
these therapies complicated cause might need to shorten telomere lengths too !
ALT+ cells using ALT
TEL+ cells using TEL

# ATRX Exon Deletion is Common in ALT
This project can be found in the Human_ATRX_ALT folder. ATRX gene mutations are found in a range of cancers. 10-15% of cancers are estimated to use ALT. ALT involves homologous recombination-based telomere elongation. Inactivating mutations in either ATRX or DAXX are found in many cancers. Depletion of ATRX seems insufficient to trigger ALT, but it does seem to play a key role in the ALT pathway. The absence of ATRX might lead to the failure of stalled replication forks to get resolved. The required fork restart would require homologus recombination and could jumpstart the ALT pathway (Clynes 2013). ALT involves a template-based lengthening of telomeres with homologous recombination. The genetic and epigenetic changes are not full understood. Lovejoy 2012 reported that ATRX gene mutations are a common feature of ALT. Specifically 19/22 ALT+ cell lines had an issue with the expression of ATRX or DAXX (Lovejoy 2012). See the Lovejoy 2012 supplementary information for the Excel table of Exon deletions in ALT cell lines. 
![ATRX_Prevents_Fork_Collapse](/Assets/ATRX_Prevents_Fork_Collapse.jpg "ATRX_Prevents_Fork_Collapse")

(Clynes 2013)

#### Getting ATRX DNA
Searching Ensembl for human ATRX yielded ATRX-201 and ATRX-202. I picked ATRX-201 cause it has 35 exons (which matches the Lovejoy 2012 paper). It was Ensembly ENST00000373344.10. Ensembl refseq switch to NCBI Reference Sequence yielded NM_000489.5 for the gene. I saved it as NM_000489.5_homo_sapiens_ATRX_Gene.fasta.

#### Removing ATRX Exons 2-29 
See the Lovejoy 2012 supplementary Excel table for a list of commonly missing ATRX Exons. I decided to play with the U2OS variant because that is a cell line that I used to grow :) U2OS is missing ATRX exons 2-29. NCBI says exon 2 is [236:348] and exon 29 is 6542..6719 https://www.ncbi.nlm.nih.gov/nuccore/NM_000489. I used R to remove those exons.

```{r}
library(seqinr)
WT_hATRX_Gene <- read.fasta("NM_000489.5_homo_sapiens_ATRX_Gene.fasta")
WT_hATRX_Gene_Nucleotides <- WT_hATRX_Gene[[1]]
length(WT_hATRX_Gene_Nucleotides)
typeof(WT_hATRX_Gene_Nucleotides)
WT_hATRX_Gene_Nucleotides
U2OS_hATRX_Gene_Nucleotide_FIRST <- WT_hATRX_Gene_Nucleotides[1:235]
U2OS_hATRX_Gene_Nucleotide_SECOND <- WT_hATRX_Gene_Nucleotides[6720:length(WT_hATRX_Gene_Nucleotides)]
U2OS_ATRX_Characters <- c(U2OS_hATRX_Gene_Nucleotide_FIRST, U2OS_hATRX_Gene_Nucleotide_SECOND)
U2OS_ATRX_DNAstring <- DNAString(paste(toupper(U2OS_ATRX_Characters), collapse = ""))
```

#### Sequence Alignment of WT ATRX to Mutant ATRX
The R msa package can't handle the full length of the ATRX gene, so I shortened it down to 400 nucleotides.
```{r}
# limit it to 400 ... that's more than enough to see exon absence
U2OS_ATRX_DNA_Short <- U2OS_ATRX_DNAstring[1:400]
WT_hATRX_DNA_Short <- toupper(WT_hATRX_Gene_Nucleotides[1:400])
write.fasta(sequences = U2OS_ATRX_DNA_Short, names = "U2OS_ATRX_DNA_Short", file.out = "U2OS_ATRX_DNA_Short.fasta", open = "w", nbchar = 70, as.string = FALSE)
write.fasta(sequences = WT_hATRX_DNA_Short, names = "WT_hATRX_DNA_Short", file.out = "WT_hATRX_DNA_Short.fasta", open = "w", nbchar = 70, as.string = FALSE)
library(msa)
# both_ATRX_Sequences <- read.fasta("WT_and_U2OS_hATRX.fasta")
both_ATRX_Sequences_SHORT <- "both_ATRX_Sequences_SHORT.fasta"
# typeof(both_ATRX_Sequences)
#both_ATRX_DNAStringSet <- readDNAStringSet(both_ATRX_Sequences)
#both_ATRX_Sequences_Alignment <- msa(both_ATRX_DNAStringSet)
both_ATRX_Sequences_SHORT_StringSet <- readDNAStringSet(both_ATRX_Sequences_SHORT)
both_ATRX_Sequences_Alignment_SHORT <- msa(both_ATRX_Sequences_SHORT_StringSet)
both_ATRX_Sequences_Alignment_SHORT
msaPrettyPrint(both_ATRX_Sequences_Alignment_SHORT, output="pdf", showNames="none",
showLogo="none", askForOverwrite=FALSE, verbose=TRUE)
#texi2pdf("both_ATRX_Sequences_Alignment.tex", clean=TRUE)
```
You can see that the sequences are the same until postion 236. That is where the Exon deletion for U2OS starts! 
![ATRX_Exon_Deletion_Alignment](/Assets/ATRX_Exon_Deletion_Alignment.jpg "ATRX_Exon_Deletion_Alignment")

# TERT Promoter Compaction is Found in ALT
ALT cells commonly have long, heterogenous telomere lengths (Kumakura 2005). Mouse embryonic stem cells deficient for DNMT have HUGE telomeres. Under normal conditions, mouse subtelomeres are heavily methylated, BUT that is not the case in mESC deficient for DNMT. The lack of DNMT increased the rate of telomeric sister chromatid exchanges (T-SCE), and ALT-associated Promyelocytic Nuclear Bodies (APBs). T-SCE and APBs are both common features of ALT activity. The authors concluded that the increased telomeric recombination MIGHT lead to telomere length changes, BUT they do not exclude the involvement of telomerase in the weirdly long telomeres that were seen (Gonzalo 2006). Luckily, I found these two other papers that go into more detail about TERT chromatin compaction in ALT!

Atkinson 2005 found that chromatin modifications of hTR and hTERT promoters were commonly found in ALT activity. Treatment of ALT+ cells with 5-AZC or Trichostatin A lead to chromatin remodeling of hTR and hTERT. This induced telomerase expression. Interestingly enough they found that mehtylated Lys20 Histon H4 was not associated with gene expression, BUT does seem to be ALT specific (Atkinson 2005). This might be a new marker of ALT activity! Acetylation of H3K9 and methylation of H3K4 is known to be associated with an open chromatin conformation. In Kumakura 2005, the authors found that ALT+ cells had H3K9 methylation and low levels of H3K4 methylation and H3K9+H3K14 acetylations. The ratio of H3K9 methylation / H3K4 methylation was different across ALT+ and TEL+ cell lines. They found that treating ALT+ cells with TSC or 5-AZC caused a reversion from complete to partial methylation of the CpG islands on the hTERT promoter. They switched an E6CL TEL+ line to TEL- and it was able to grow for well over 240 population doublings (Kumakura 2005). That's some ALT activity right there!
![H3K9_H3K4_Methylation_Ratio](/Assets/H3K9_H3K4_Methylation_Ratio.jpg "H3K9_H3K4_Methylation_Ratio")

(Kumakura 2005)

#### Getting The hTERT Sequence
The UniProt entry O14746 is for hTERT https://www.uniprot.org/uniprot/O14746. Following GeneID 7015 gets TERT telomerase reverse transcriptase for humans https://www.ncbi.nlm.nih.gov/gene/7015. Note that the reverse arrow on TERT indicates that the sequence is on the reverse strand (this will become important later). I downloaded the FASTA as "NC_000005.10_hChrom5_TERT_CpG_Start.fasta". A quick text search shows that the start codon is at position 59. Take care with this sequence cause it's the reverse complement of the actual sequence!

![hTERT_NCBI_Reverse_Strand](/Assets/hTERT_NCBI_Reverse_Strand.jpg "hTERT_NCBI_Reverse_Strand")

#### Obtaining hTERT WITH the CpG Island Region
Stay with me ... we're about to dig a bit into the literature! The hTERT sequence that I grabbed from NCBI DOES NOT contain the CpG island for hTERT. It doesn't even contain the normal promoter region for hTERT! Cong 1999 reports that the core hTERT promoter region is from -330 to +361 bp of the ATG start codon. HOWEVER, Kumakura 2005 found the hTERT CpG island to be from 654 bp upstream to 510 bp downstream of the ATG start codon, so this is the actual region that I need to grab! 

Grabbing the hTERT FASTA sequence from https://www.ncbi.nlm.nih.gov/gene/7015 INITIALLY is from: 1253167 to: 1295047. Checking Cong 1999 and The FASTA file, I can see that the hTERT start codon, AND a bit more of that region, of "ATGCCGCGCGCT" is at the end of the first FASTA line, which is 59 in from the left (59 is A of ATG). CpG is 654 bp upstream of the transcriptoin start site, SO going 595 back from current start site 1253167-595 = 1252572. 1252572 should be the start of the CpG island, RIGHT?!? WRONG!!! ... Why aren't I getting any more nucleotides before the current start read?!?!? OH!!! I'm looking at the reverse complement, lol! ;) 

It should be  1295047 + 595 = 1295642 YESSSSSSS, that's right :) Now I have more at the beginniig, so "atgccgcgcgctccccgct" is fully searchable! This is the new range:
https://www.ncbi.nlm.nih.gov/nuccore/NC_000005.10?report=fasta&from=1253167&to=1295642&strand=true. I saved the FASTA as NC_000005.10_hChrom5_TERT_CpG_Start.fasta. 

#### Analyzing the Alleged CpG Promoter Region
Dessain 2000 reported that the hTERT CpG island has a GC content of 74% and a CG:GC ratio of 0.87. Is that what I get for the same region?!? I wrote code in R to get the GC content and CG:GC ratio of the hTERT CpG promoter region that I identified. I didn't comment my code ... I am sorry. Note that I picked i=1164 cause Kumakura 2005 says that the hTERT CpG island to be from 654 bp upstream to 510 bp downstream of the ATG start codon, which is 654 + 510 = 1164 :) Here's R code that I didn't bother commenting :( 

```{r}
dna_file <- read.fasta("/media/david/Linux/Introns_Exons_and_Promoters/hTERT_CpG_Island/NC_000005.10_hChrom5_TERT_CpG_Start.fasta")
individual_characters <- dna_file[[1]]
i <- 1
CG_count <- 0
GC_count <- 0
g_count <- 0
c_count <- 0
for (letter in individual_characters) {
  #print(letter)
  i <- i + 1
  if (letter == "g") {
    g_count <- g_count + 1
  }
  if (letter == "c") {
    c_count <- c_count + 1
  }
  if ((letter == "c") & (last_letter == "g")) {
    GC_count <- GC_count + 1
  }
  if ((letter == "g") & (last_letter == "c")) {
    CG_count <- CG_count + 1
  }
  last_letter <- letter
  if (i >= 1164) {
    break
  }
}
print(100*(c_count+g_count)/1164)
print(GC_count)
print(CG_count)
print(CG_count/GC_count)
```
The lazily unlabeled output is:

[1] 76.03093
[1] 167
[1] 141
[1] 0.8443114

Cool-ness! My GC content is 76 % and the ratio of CpG/GpC is 0.84. Recall that Dessain 2000 reported that the hTERT CpG island has a GC content of 74% and a CG:GC ratio of 0.87. I'm 2 % off of the GC content and 0.03 off of the CpG/GpC. THAT'S PRETTY GOOD FOR REPLICATING DATA FROM A PAPER THAT IS ALMOST TWO DECADES OLD :D But, what if that was just random luck? I re-ran that same R code on the region AFTER the CpG island and I got wildly different data. Here's that code (yes, I should've used a function; I am lazy, lol):

```{r}
dna_file <- read.fasta("/media/david/Linux/Introns_Exons_and_Promoters/hTERT_CpG_Island/NC_000005.10_hChrom5_TERT_CpG_Start.fasta")
individual_characters <- dna_file[[1]]
#individual_characters[5]
i <- 1164
CG_count <- 0
GC_count <- 0
g_count <- 0
c_count <- 0
for (letter in individual_characters) {
  if (i <1164) {
    next
  }
  #print(letter)
  i <- i + 1
  if (letter == "g") {
    g_count <- g_count + 1
  }
  if (letter == "c") {
    c_count <- c_count + 1
  }
  if ((letter == "c") & (last_letter == "g")) {
    GC_count <- GC_count + 1
  }
  if ((letter == "g") & (last_letter == "c")) {
    CG_count <- CG_count + 1
  }
  last_letter <- letter
  
  
  if (i >= 42476) {
    break
  }
}
print(100*(c_count+g_count)/41312)
print(GC_count)
print(CG_count)
print(CG_count/GC_count)
```

The lazily unlabeled output is:

[1] 58.55926
[1] 3047
[1] 1654
[1] 0.542829

My GC content is 58.5 % and the ratio of CpG/GpC is 0.54. Recall that Dessain 2000 reported that the hTERT CpG island has a GC content of 74% and a CG:GC ratio of 0.87. THE REGION THAT IS NOT A CpG ISLAND IS 15.5% off of the GC content and 0.33 off of the CpG/GpC. I could dig into this with more statistical rigor, but I think you get the idea. I'M EXCITED!!! This was a really cool biological programming exercise!!! :D

 



# POT-1 Deficiency Creates ALT+ C. Elegans Strains
Telomeres cap linear chromosomes because DNA polymerase can't completely copy chromosomes. Telomerase adds telomeric repeats to the ends of linear chromosomes with reverse transcription. Most human cancers have long, heterogenous telomeres. Telomere shortening leads to senescence and potentially crisis. Cancer emerges as part of massive cell death and genomic rearrangements after crisis. 10-15% of cancers are estimated to use ALT (Cheng 2012). 

ALT can happen in Caenorhabditis elegans! Mammalian POT1 has homologs in C. elegans as pot-1 (CeOB2) and pot-2 (CeOB1). What's the deal with the reversing of 1 and 2? That's how it's reported in the paper ... it's odd. pot-1 mutant C. elegans have HUGE telomere lengths while pot-2 mutants have normal telomere lengths. The authors of Cheng 2012 created a variety of mutants in C. elegans.  The trt-1 C. elegans mutant has a deletion in telomerase reverse transcriptase. trt-1 & pot-2 absence led to ALT+ Caenorhabditis elegans with normal telomere lengths. trt-1 and pot-1 mutants were found to have long, heterogenous telomere lengths like those seen in human ALT. Here is the survival figure showing that C. elegans can survive in the absence of telomerase reverse transcriptase.

![Celegans_ALT_Generation_Survival](/Assets/Celegans_ALT_Generation_Survival.jpg "Celegans_ALT_Generation_Survival")

(Cheng 2012)

#### Multiple Sequence Alignment of pot-1 Genes
YES, pot-2 was the central point of the paper, but it won't be as fun to play with because it only has one isoform. I picked pot-1 cause there is a lot of cool stuff to play with. There were a lot of workup steps to get all of the sequences ... It would take a long while to review them. Essentially, I looked up the proteins on UniProt and then grabbed the DNA files from NCBI GenBank and WormBase. Check out the Celegans_POT1_ALT folder for the file names of everything. The file containing all the C. elegans genes is Celegans_POT1_genes.fasta. I used the R package "msa" for multiple sequence alignment with this code:

```{r}
library(msa)
Celegans_POT1_genes <- "/media/david/Linux/Introns_Exons_and_Promoters/Celegans_POT1_ALT/DNA/Celegans_POT1_genes.fasta"
Celegans_POT1_genes_DNA <- readDNAStringSet(Celegans_POT1_genes)
Celegans_POT1_gene_alignment <- msa(Celegans_POT1_genes_DNA)
msaPrettyPrint(Celegans_POT1_gene_alignment, output="pdf", showNames="none",
showLogo="none", askForOverwrite=FALSE, verbose=FALSE)
```

The aligned sequences aren't very pretty ... I decided not to include sequence labels cause it shortened the available nucleotide space for each new line. Here's part of the output for you to get the idea of the work:

![Celegans_POT1_gene_alignment](/Assets/Celegans_POT1_gene_alignment.jpg "Celegans_POT1_gene_alignment")

#### Multiple Sequence Alignment of pot-1 Proteins

#### Displaying pot-1 Open Reading Frames

# STN1 Mutation Triggers ALT in Yeast
BLASTAlignRetrieve/ID mappingPeptide searchContactHelp
UniProtKB - Q9H668 (STN1_HUMAN)
PART OF CST COMPLEX AND ALT-INVOLVED!!!!!!
https://www.uniprot.org/uniprot/Q9H668

FIG4 of paper shows sequence analysis of K. lactis STN1 gene and homologues from S. cerevisiae (Sc) and Candida glabrata (Cgl) (GenBank accession numbers P_38960 and XP_448655, respectively).
BUT NCBI says "The following term was not found in Nucleotide: P_38960."for P_38960. XP_448655 is here: https://www.ncbi.nlm.nih.gov/protein/XP_448655. 
I couldn't find the sequence for K. lactis they were talking about ... I think it's NCBI Reference Sequence: XM_452728.1 Kluyveromyces lactis uncharacterized protein (KLLA0_C11825g), partial mRNA BECAUSE /note="weakly similar to uniprot|P38960 Saccharomyces
                     cerevisiae YDR082W STN1 Protein involved in telomere
                     length regulation functions in telomere metabolism during
                     late S phase"
and it's the only result (other than whole chromosome chunks) w/ this crazy search term I made:
(((stn1) NOT "Pyrenophora tritici-repentis"[porgn:__txid45151] NOT "Fusarium fujikuroi"[porgn:__txid5127] NOT "[Candida] glabrata"[porgn:__txid5478] NOT "Hortaea werneckii"[porgn:__txid91943] NOT "Saccharomyces cerevisiae"[porgn:__txid4932]) NOT "Metarhizium robertsii"[porgn:__txid568076] NOT "Fusarium sp. FIESC_5 CS3069"[porgn:__txid1318460] NOT "Fusarium pseudograminearum CS3487"[porgn:__txid1318458] NOT "Fusarium pseudograminearum CS3427"[porgn:__txid1318457] NOT "Fusarium pseudograminearum CS3220"[porgn:__txid1318456] NOT "Fonsecaea multimorphosa"[porgn:__txid979981] NOT "Cladophialophora immunda"[porgn:__txid569365] NOT "Aspergillus nidulans FGSC A4"[porgn:__txid227321] NOT "Candida viswanathii"[porgn:__txid5486] NOT "Zygosaccharomyces bailii"[porgn:__txid4954] NOT "Metarhizium anisopliae"[porgn:__txid5530] NOT "Aspergillus flavus"[porgn:__txid5059] NOT "Talaromyces atroroseus"[porgn:__txid1441469] NOT "[Candida] auris"[porgn:__txid498019] NOT "Zygosaccharomyces rouxii"[porgn:__txid4956] NOT "[Candida] boidinii"[porgn:__txid5477] NOT "Komagataella phaffii"[porgn:__txid460519] NOT "Aspergillus fumigatus"[porgn:__txid746128] NOT "Candida albicans SC5314"[porgn:__txid237561] NOT "Yarrowia lipolytica"[porgn:__txid4952]) AND "Kluyveromyces lactis"[porgn:__txid28985] 

sceriviasa yeast was easier to find cause it's got stn1p in the title
https://www.ncbi.nlm.nih.gov/nuccore/NM_001180390.1


Lyer 2005 A Mutation in the STN1 Gene Triggers an Alternative Lengthening of Telomere-Like Runaway Recombinational Telomere Elongation and Rapid Deletion in Yeast
https://mcb.asm.org/content/25/18/8064

#### getting k lactis

#### replicating attempt of lyer 2005

#### multiple sequence3 alignment of human STN1 protein vs. k lactis



# Data Sources

# Citations
* Cheng 2012 Caenorhabditis elegans POT-2 telomere protein represses a mode of alternative lengthening of telomeres with normal telomere lengths
* Cong 1999 The human telomerase catalytic subunit hTERT: organization of the gene and characterization of the promoter
* Lovejoy 2012 PLoS Genet Loss of ATRX, genome instability, altered DNA damage response hallmarks of ALT pathway
* Clynes 2013 Curr Opin Genet Dev ATRX and the replication of structured DNA
* Gonzalo 2006 DNA methyltransferases control telomere length and telomere recombination in mammalian cells.pdf
* Atkinson 2005 Lack of Telomerase Gene Expression in Alternative Lengthening of Telomere Cells Is Associated with Chromatin Remodeling of the hTR and hTERT Gene Promoters
* Kumakura 2005 Reversible Conversion of Immortal Human Cells from Telomerase-Positive to Telomerase-Negative Cells
* Cong 1999 The human telomerase catalytic subunit hTERT: organization of the gene and characterization of the promoter
* Dessain 2000 Methylation of the Human Telomerase Gene CpG Island
* Cheng 2012 Caenorhabditis elegans POT-2 telomere protein represses a mode of alternative lengthening of telomeres with normal telomere lengths

NOTES NOTES NOTES NOTES NOTES NOTES NOTES NOTES NOTES NOTES NOTES NOTES NOTES NOTES NOTES NOTES NOTES NOTES NOTES NOTES 
* Lovejoy 2012 PLoS Genet Loss of ATRX, genome instability, altered DNA damage response hallmarks of ALT pathway
	* ALT proposedtemplate extension of teloemres w/ HR, BUT getic / epigenetic changes aren't fully known. this paper used genomic, molecular  biological and cell biological sciencyi9ng on panel of 22 ALT cell lines (haz lines derived in vitro). loss of ATRX protein + ATRX gene mutations are feature of ALT lines. ALT linked w/ genome rearrangements, micronucleation and wacky DSB repair ... these diagnositc features may help to identify ALT.
	* table of exon deletions from supplementary excel
* Clynes 2013 Curr Opin Genet Dev ATRX and the replication of structured DNA(1)
	* ATRX mutations found in vareity of cancer types. 10-15% cancers use ALT. HR w/ telomeres OR ectr. inactivating mutations in ATRX OR DAXX found in many ALT cancers. loss of ATRX found in 90% of ALT cell lines. Depletion of ALT seems insufficient to cause ALT ... + telomerase inhibition? 
	* absence of ATRX might lead to stalled replication forks are not processed well. fork restart is depedent on HR and could trigger the ALT pathway FIGURE OF ALT PATHWAY FROM CLYNES 2013



* Gonzalo 2006 DNA methyltransferases control telomere length and telomere recombination in mammalian cells.pdf
	* mESC deficient for DNMT have HUGE telomeres. Mouse subtelomeres HEAVILY methylated, BT NOT in DNMT cells?!?	
	* Histone 3 Lys9 H3K9 and histone 4 Lys 20 H4K20 trimethylation REMAIN.
	* lack of DNMT INCREASED telomeric recombination (T-SCE) and APBs. authors conclude that increased telomeric recombination MIGHT lead to telomere length changes, BUT they do not exclude the invovlement of telomerase in the aberrant telomere elongation, BUT here are two papers that better describe the link!!!!
* Kumakura 2005 Reversible Conversion of Immortal Human Cells from Telomerase-Positive to Telomerase-Negative Cells
	* Triochostatin A and/or 5 AZC reversion from complete to partial methylation of CpG islands of hTERT promoters.
		* ALT+ cells had long, heterogenous telomere lengths
		* acetylation H3K9 and methylation of H3K4 OPEN chromatin KNOWN
		* parent WHE NO TMM had H3K9 hypo methylation, H3K4 hypermethylation and acetylations of H3K9 and H3K14
		* ALT H3K9 methylation and low levels of H3K4 methylation and H3K9+H3K14 acetylations
		* treating with TSC or 5AZC decreased ratio of H3K9 methylation / H3K4 methylation of ALT into TEL- lines
		* authors made E6CL from TEL+ to TEL- w/ >240 PDs!!!
		* evidence that TMM reverseible based on chromatin structure!!!
		* FIGURE OF RATIOS!
* Atkinson 2005 found chromatin modifications of hTR and hTERT promoters linked w/ ALT
	* hTR and hTERT lower expression in ALT is linked w/ histone H3 and H4 hypoacetylation amd methylation of Lys9 of Histone H3
	* TEL+ cells had hyperacetylation of H3 and H4 and methylation of Lys4H3. 
	* Treatment w/ 5 AZC and Trichostatin A => chromatin remodeling hTR, hTERT and therefore expression. 
	* ALT+ methylated Lys20 Histone H4 not associated w/ gene expression, but does seem ALT-specific ... new marker?
	* authors propose ALT may arise from tight repression of hTR and hTERT promoters ... possibly why MSC go ALT?!?




#### 
Cong 1999 The human telomerase catalytic subunit hTERT: organization of the gene and characterization of the promoter
Dessain 2000 Methylation of the Human Telomerase Gene CpG Island

#### Identifying the hTERT CpG Island Region
	Kumakura 2005 found the hTERT CpG island to be from 654 bp upstream to 510 bp downstream of the ATG start codon. 

	Cong 1999 The human telomerase catalytic subunit hTERT: organization of the gene and characterization of the promoter
	Kumakura 2005 Reversible Conversion of Immortal Human Cells from Telomerase-Positive to Telomerase-Negative Cells
	I want to explore hTERT promoter CpG island methylation, BUT I'll need to do some digging to get the sequence and identify the CpG island region. Cong 1999 reports that the core hTERT promoter region is from -330 to +361 bp of the ATG start codon. HOWEVER, that doesn't necessarily imply the CpG region to be ONLY from -330 to +361bp (see next section).

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

Cheng 2012 Caenorhabditis elegans POT-2 telomere protein represses a mode of alternative lengthening of telomeres with normal telomere lengths
	Most human cancers have long/heterogenous telomeres 
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

