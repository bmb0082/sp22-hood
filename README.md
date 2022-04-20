# Using bioinformatics to detect SNPs and predict functional changes in COX1, CytB, and the mitochondrial control region

## Project Summary

This project’s central aim is to identify and characterize mutations from reproductive female (n=3) and 14-day old (n=3) House mice maintained in seminatural enclosures. Given mutations in the mitochondrial genome can arise from replication error in the form of transitions or from damaging reactive oxygen species (ROS)as a unique transversion type, we further investigate the origin of the mutations detected and use protein modelling to predict which type of mutation could induce functional changes. The project begins with a proof of concept providing a background in single nucleotide polymorphisms (SNP), theory, and the central dogma to show how insertions and deletions might affect translated products while showcasing functional changes using a predicted 3D model in COX1. While methods for the proof of concept used reference genes from the NCBI database, we used bioinformatics to identify mutations in COX1, CytB, and the control region of the mitochondrial genome from raw read fastq-formatted data. Of the mutations detected in COX1, CytB, and the control region, 24.0%, 24.4%, and 25.0% were ROS-mediated respectively, while 33.0%, 30.8%, and 33.7% of mutations detected were due to errors in replication. Although none of the mutations detected within COX1 and CytB were predicted to result in functional changes, we are left to ask if the mutational load detected in the control region might have any effect on mtDNA replication frequency.

# Table of Contents

<details>
<summary>Sections</summary>

 1. [Modeling transcript directionality in *Mus musculus* mitochondrial cytochrome c oxidase subunit I](https://github.com/bmb0082/lab-project/blob/main/README.md#modeling-transcript-directionality-in-mus-musculus-mitochondrial-cytochrome-c-oxidase-subunit-i)
 2. [Modeling insertions, deletions, and INDEL mutations in *Mus musculus* mitochondrial cytochrome c oxidase subunit I](https://github.com/bmb0082/lab-project/blob/main/README.md#modeling-insertions-deletions-and-indel-mutations-in-mus-musculus-mitochondrial-cytochrome-c-oxidase-subunit-i)
 3. [Analyzing SNP transition and transversion ratios in three *Mus musculus* mitochondrial genes: cytochrome c oxidase subunit I, cytochrome b subunit, and the D-loop control region](https://github.com/bmb0082/lab-project/blob/main/README.md#analyzing-snp-transition-and-transversion-ratios-in-three-mus-musculus-mitochondrial-genes-cytochrome-c-oxidase-subunit-i-cytochrome-b-subunit-and-the-d-loop-control-region)
 4. [Investigating the origin of SNPs by type in three *Mus musculus* mitochondrial genes: cytochrome c oxidase subunit I, cytochrome b subunit, and the D-loop control region](https://github.com/bmb0082/lab-project/blob/main/README.md#investigating-the-origins-of-snps-by-type-in-three-mus-musculus-mitochondrial-genes-cytochrome-c-oxidase-subunit-i-cytochrome-b-subunit-and-the-d-loop-control-region)
 5. [Modeling the effects of SNPs from *Mus musculus* consensus read data in three mitochondrial genes: cytochrome c oxidase subunit I, cytochrome b subunit, and the D-loop control region](https://github.com/bmb0082/lab-project/blob/main/README.md#modeling-the-effects-of-snps-from-mus-musculus-consensus-read-data-in-three-mitochondrial-genes-cytochrome-c-oxidase-subunit-i-cytochrome-b-subunit-and-the-d-loop-control-region)
 
</details>

# Proof of Concept: Modeling transcript directionality in *Mus musculus* mitochondrial cytochrome c oxidase subunit I
Mitochondrial cytochrome c oxidase I (COX1), encoded by the mt-Co1 gene, is a subunit of cytochrome c oxidase (RC-IV), the last enzyme in the electron transport chain. Correct native protein structure is essential for safe acceptance of terminal electons, step-wise reduction of oxygen to peroxide to water, and overall mitochondiral efficiency for adenosine triphosphate (ATP) production. As a result, errors in COX1 can trap electrons in the metabolically dangerous Q-cycle of respiratory complex-III (RC-III), where single electrons are frequently transferred between ubiquinones, semiquinone radicals, and ubiquinol molecules. If the ability of RC-IV to accept electrons from cytochrome C is affected, it may result in increased production of reactive oxygen species from RC-III. Additionally, if RC-IV is unable to participate in cytochrome c oxidation, it can lock cytochrome C in its reduced state. With nowhere to go, this increases the potential of cytochrome c leakage, which activates cellular apoptosis pathways when present in the cytoplasm (Garrido et al.)

Cytochrome c oxidase deficiency causes many problems, such as increased ROS production and possible activation of controlled cell death. However, it is possible for organisms to manage these outcomes. Mutations that affect COX1 are one of the most common causes of genetic mitochondrial disorders, presenting as multi-organ, heterogeneous symptoms depending on the level of mitochondrial efficiency (Rak, Malgorzata et al.).

## Purpose - To understand the importance of direction when modeling proteins from translated sequences
**Hypothesis - Changing the direction that the mRNA transcript is read will alter protein structure and function**

## Methodology
Below is the selected reference for the mitochondrial gene encoding subunit I of cytochrome c oxidase. I chose this gene because COX1 contains the catalytic unit of cytochrome c oxidase, housing two copper ions as shown in **FIG-1**.  Therefore, it's correct structure and function are vital for the enzyme's overall activity. I wanted to see how differently the protein would fold depending on the directionality of the translated sequence.
```
Reference genome source of mt-Co1 gene, presented in the 5'-3' direction:
Mus musculus mitochondrion, complete genome
NCBI Reference Sequence: NC_005089.1

GenBank Graphics
>NC_005089.1:5328-6872 Mus musculus mitochondrion, complete genome
ATGTTCATTAATCGTTGATTATTCTCAACCAATCACAAAGATATCGGAACCCTCTATCTACTATTCGGAG
CCTGAGCGGGAATAGTGGGTACTGCACTAAGTATTTTAATTCGAGCAGAATTAGGTCAACCAGGTGCACT
TTTAGGAGATGACCAAATTTACAATGTTATCGTAACTGCCCATGCTTTTGTTATAATTTTCTTCATAGTA
ATACCAATAATAATTGGAGGCTTTGGAAACTGACTTGTCCCACTAATAATCGGAGCCCCAGATATAGCAT
TCCCACGAATAAATAATATAAGTTTTTGACTCCTACCACCATCATTTCTCCTTCTCCTAGCATCATCAAT
AGTAGAAGCAGGAGCAGGAACAGGATGAACAGTCTACCCACCTCTAGCCGGAAATCTAGCCCATGCAGGA
GCATCAGTAGACCTAACAATTTTCTCCCTTCATTTAGCTGGAGTGTCATCTATTTTAGGTGCAATTAATT
TTATTACCACTATTATCAACATGAAACCCCCAGCCATAACACAGTATCAAACTCCACTATTTGTCTGATC
CGTACTTATTACAGCCGTACTGCTCCTATTATCACTACCAGTGCTAGCCGCAGGCATTACTATACTACTA
ACAGACCGCAACCTAAACACAACTTTCTTTGATCCCGCTGGAGGAGGGGACCCAATTCTCTACCAGCATC
TGTTCTGATTCTTTGGGCACCCAGAAGTTTATATTCTTATCCTCCCAGGATTTGGAATTATTTCACATGT
AGTTACTTACTACTCCGGAAAAAAAGAACCTTTCGGCTATATAGGAATAGTATGAGCAATAATGTCTATT
GGCTTTCTAGGCTTTATTGTATGAGCCCACCACATATTCACAGTAGGATTAGATGTAGACACACGAGCTT
ACTTTACATCAGCCACTATAATTATCGCAATTCCTACCGGTGTCAAAGTATTTAGCTGACTTGCAACCCT
ACACGGAGGTAATATTAAATGATCTCCAGCTATACTATGAGCCTTAGGCTTTATTTTCTTATTTACAGTT
GGTGGTCTAACCGGAATTGTTTTATCCAACTCATCCCTTGACATCGTGCTTCACGATACATACTATGTAG
TAGCCCATTTCCACTATGTTCTATCAATGGGAGCAGTGTTTGCTATCATAGCAGGATTTGTTCACTGATT
CCCATTATTTTCAGGCTTCACCCTAGATGACACATGAGCAAAAGCCCACTTCGCCATCATATTCGTAGGA
GTAAACATAACATTCTTCCCTCAACATTTCCTGGGCCTTTCAGGAATACCACGACGCTACTCAGACTACC
CAGATGCTTACACCACATGAAACACTGTCTCTTCTATAGGATCATTTATTTCACTAACAGCTGTTCTCAT
CATGATCTTTATAATTTGAGAGGCCTTTGCTTCAAAACGAGAAGTAATATCAGTATCGTATGCTTCAACA
AATTTAGAATGACTTCATGGCTGCCCTCCACCATATCACACATTCGAGGAACCAACCTATGTAAAAGTAA
AATAA
```

### Predicted model of COX1, presented 5'-3'
Reading the translated sequence 5'-3' is the correct direction, as it generates a defined structure of cytochrome c oxidase I. [ExPASy Translate Tool](https://web.expasy.org/translate/) recognized the entire gene in this direction as a single open reading frame, and [SWISS-MODEL](https://swissmodel.expasy.org/interactive/) generated a 94.16% sequence homology to bovine heart cytochrome c oxidase subunit I. Likewise, it models a successful protein-ligand interaction pipeline between the active site and a peroxide ion, shown as two red dots in the catalytic core.

<p align="center">
<img src="https://user-images.githubusercontent.com/98036665/151746598-19682dd3-16ae-4470-91c8-03b9d1dadf38.png" width="600" height="600">
</p>

```
Open reading frame 1 sequence, translated 5'-3':
MFINRWLFSTNHKDIGTLYLLFGAWAGMVGTALSILIRAELGQPGALLGDDQIYNVIVTAHAFVMIFFMVMPMMIGGFGNWLVPLMIGAPDMAFPRMNNMSFWLLPPSF
LLLLASSMVEAGAGTGWTVYPPLAGNLAHAGASVDLTIFSLHLAGVSSILGAINFITTIINMKPPAMTQYQTPLFVWSVLITAVLLLLSLPVLAAGITMLLTDRNLNTT
FFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGIISHVVTYYSGKKEPFGYMGMVWAMMSIGFLGFIVWAHHMFTVGLDVDTRAYFTSATMIIAIPTGVKVFSWLATL
HGGNIKWSPAMLWALGFIFLFTVGGLTGIVLSNSSLDIVLHDTYYVVAHFHYVLSMGAVFAIMAGFVHWFPLFSGFTLDDTWAKAHFAIMFVGVNMTFFPQHFLGLSGM
PRRYSDYPDAYTTWNTVSSMGSFISLTAVLIMIFMIWEAFASKREVMSVSYASTNLEWLHGCPPPYHTFEEPTYVKVK-
```

### Predicted model of COX1, presented 3'-5'
This is not the correct folding conformation of cytochrome c oxidase I. Reading the translated sequence 3'-5' results in a polypeptide that lacks structural integrity and the binding pockets for heme-a, and the heme-a3-CuB. This protein cannot participate in the transient reduction of the heme-a3/CuB site and stepwise reduction of oxygen to water for successful RC-IV function.

<p align="center">
<img src="https://user-images.githubusercontent.com/98036665/150719390-4b9598a1-a27a-4266-84ca-448126c0dbc9.png" width="600" height="600">
</p>

```
Open reading frame 1 sequence, tranlated 3'-5':
LFYFYMGWFLECVMWW-AAMKSF-ICWSMRYWYYFSFWSKGLSNYKDHDENSC-WNKWSY--DSVSCGVSIWVVWVASWYSWKAQEML-EECYVYSYEYDGEVGFCSCV
I-GEAWK-WESVNKSCYDSKHCSHW-NMVEMGYYMVCIVKHDVKGWVG-NNSG-TTNCK-ENKA-GS-YSW-SFNITSV-GCKSAKYFDTG-NCDNYSGWCKVSSCVYI
-SYCEYVVGSYNKA-KAN-HYCSYYSYMAE-FFFSGVVSNYMWNNSKSWEDKNMNFWVPKESEQMLVENWVPSSSGIKESCV-VAVC--YSNACG-HW-W--EQYGCNK
YGSDK-WSLMLCYGWGFHVDNSGNKINCT-N-WHSS-MKGENC-VYWCSCMG-ISG--WVDCSSCSCSCFYYWWC-EKEKWWW-ESKTYIIYSWECYIWGSDY-WDKSV
SKASNYYWYYYEENYNKSMGSYDNIVNLVIS-KCTWLT-FCSN-NT-CSTHYSRSGSE--MEGSDIFVIGWE-STINEH
```
## Results and Discussion
The data gathered in this experiment supports the idea that directionality of translation is an extremely important factor for protein structure and function. Most mRNA transcripts are meant to be read in the 5'-3' direction, but some genes can encode two or more different functional proteins, one in the 5'-3' direction and one another in the 3'-5' direction. For biologically important genes with evolutionarily conserved functions, it is critical to also conserve directionality. Reading of a 5'-3' transcript in the 3'-5' direction results in a completely different primary sequence with completely different function, if any at all. In the case of COX1, reading in the 3'-5' direction erases function of the gene. Understanding directionality is an important factor when choosing to submit sequences to predictive modeling software.

## Conclusions
Reading the mt-CO1 gene transcript in 5'-3' direction produces a functional protein that can bind its necessary prosthetic groups and participate in the safe acceptance of electrons. Consequently, the 3'-5' direction lacks function given its lack of folding, and would not perform catalysis, as it lacks an active site.

## References
1. Garrido, C, Galluzzi, L, Brunet, M et al. "Mechanisms of cytochrome c release from mitochondria." Cell Death Differ 13, 1423–1433 (2006). doi:10.1038/sj.cdd.4401950
2. Rak, Malgorzata et al. “Mitochondrial cytochrome c oxidase deficiency.” Clinical science (London, England : 1979) vol. 130,6 (2016): 393-407. doi:10.1042/CS20150707

[Return to Table of Contents](https://github.com/bmb0082/lab-project/blob/main/README.md#table-of-contents)

# Modeling insertions, deletions, and INDEL mutations in *Mus musculus* mitochondrial cytochrome c oxidase subunit I
Mutations are changes in DNA sequence that result in a variant form that can be passed down to subsequent generations. Most commonly, alterations in the structure of a gene are caused by single nucleotide transition or transversion mutations called single nucleotide polymorphisms, or SNPs (Hunt, Ryan et al.). However, mutations can also arise from insertion, deletion, or rearrangement of sections of a gene or chromosome.

Mutations do not always have to cause a change in phenotype. These mutations are called neutral mutations, producing *synonymous proteins* that have the same function as the wild type. This can occur when single nucleotide changes encode the same amino acid or a residue with similar enough electrostatic properties to produce a protein with unaltered function. Alternatively, mutations that produce *non-synonymous* proteins, or proteins with altered function, usually, but not always, result from more drastic changes in the genome. For example, frameshift mutations such as insertions or deletions or rearrangement mutations at the chromosomal level can result in a completely different downstream primary structures, causing improper folding and loss of function.

In this project, all mutations were made at the amino acid level by inserting or deleting one or more amino acids (or multiples of 3 nucleotides aligned as codons) that may or may not affect the overall function of the protein.

## Purpose - To investigate the effects of amino acid mutations on predicted protein structure
**Hypothesis - Mutations in the active site of an enzyme will decrease the enzyme's affinity for substrate or cause loss of function**

## Methodology
The wild-type gene encoding COX1 was translated in the 5'-3' direction using [ExPASy Translate Tool](https://web.expasy.org/translate/) and modeled based on homology using [SWISS-MODEL](https://swissmodel.expasy.org). Based on the model of the wild-type control COX1 protein, the residues involved in the active site were identified as H240, V243, H290, and H291. I chose to create mutations in these amino acids in order to observe how they would affect ligand binding, as the active site residues are the most indicative of protein function.

An insertion mutation, deletion mutation, and an insertion/deletion mutation respectively were added directly into the translated sequence and modeled as predicted proteins, producing representational models for all three to observe physical changes in protein structure. Dotplots comparing nucleotide sequences of each mutation against the control were created using [EMBOSS](https://www.bioinformatics.nl/cgi-bin/emboss/dotmatcher) to view the differences between the sequences through a visual representation of what occurs for each mutation and where it occurs in the nucleotide sequence. Alignment scores of the nucleotide and amino acid sequences were producing using [BLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) and [BLASTp](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) respectively to view sequence homology compared to the control.

Below is the selected nucleotide reference sequence for the cytochrome c oxidase subunit I, its 5'-3' translated sequence, and the modeled protein prediction that will serve as the control for the mutated proteins.

<pre>
Reference genome source of mt-Co1 gene, presented in the 5'-3' direction:
Mus musculus mitochondrion, complete genome
NCBI Reference Sequence: NC_005089.1

GenBank Graphics
>NC_005089.1:5328-6872 Mus musculus mitochondrion, complete genome
ATGTTCATTAATCGTTGATTATTCTCAACCAATCACAAAGATATCGGAACCCTCTATCTACTATTCGGAG
CCTGAGCGGGAATAGTGGGTACTGCACTAAGTATTTTAATTCGAGCAGAATTAGGTCAACCAGGTGCACT
TTTAGGAGATGACCAAATTTACAATGTTATCGTAACTGCCCATGCTTTTGTTATAATTTTCTTCATAGTA
ATACCAATAATAATTGGAGGCTTTGGAAACTGACTTGTCCCACTAATAATCGGAGCCCCAGATATAGCAT
TCCCACGAATAAATAATATAAGTTTTTGACTCCTACCACCATCATTTCTCCTTCTCCTAGCATCATCAAT
AGTAGAAGCAGGAGCAGGAACAGGATGAACAGTCTACCCACCTCTAGCCGGAAATCTAGCCCATGCAGGA
GCATCAGTAGACCTAACAATTTTCTCCCTTCATTTAGCTGGAGTGTCATCTATTTTAGGTGCAATTAATT
TTATTACCACTATTATCAACATGAAACCCCCAGCCATAACACAGTATCAAACTCCACTATTTGTCTGATC
CGTACTTATTACAGCCGTACTGCTCCTATTATCACTACCAGTGCTAGCCGCAGGCATTACTATACTACTA
ACAGACCGCAACCTAAACACAACTTTCTTTGATCCCGCTGGAGGAGGGGACCCAATTCTCTACCAGCATC
TGTTCTGATTCTTTGGGCACCCAGAAGTTTATATTCTTATCCTCCCAGGATTTGGAATTATTTCACATGT
AGTTACTTACTACTCCGGAAAAAAAGAACCTTTCGGCTATATAGGAATAGTATGAGCAATAATGTCTATT
GGCTTTCTAGGCTTTATTGTATGAGCC<b>CACCAC</b>ATATTCACAGTAGGATTAGATGTAGACACACGAGCTT
ACTTTACATCAGCCACTATAATTATCGCAATTCCTACCGGTGTCAAAGTATTTAGCTGACTTGCAACCCT
ACACGGAGGTAATATTAAATGATCTCCAGCTATACTATGAGCCTTAGGCTTTATTTTCTTATTTACAGTT
GGTGGTCTAACCGGAATTGTTTTATCCAACTCATCCCTTGACATCGTGCTTCACGATACATACTATGTAG
TAGCCCATTTCCACTATGTTCTATCAATGGGAGCAGTGTTTGCTATCATAGCAGGATTTGTTCACTGATT
CCCATTATTTTCAGGCTTCACCCTAGATGACACATGAGCAAAAGCCCACTTCGCCATCATATTCGTAGGA
GTAAACATAACATTCTTCCCTCAACATTTCCTGGGCCTTTCAGGAATACCACGACGCTACTCAGACTACC
CAGATGCTTACACCACATGAAACACTGTCTCTTCTATAGGATCATTTATTTCACTAACAGCTGTTCTCAT
CATGATCTTTATAATTTGAGAGGCCTTTGCTTCAAAACGAGAAGTAATATCAGTATCGTATGCTTCAACA
AATTTAGAATGACTTCATGGCTGCCCTCCACCATATCACACATTCGAGGAACCAACCTATGTAAAAGTAA
AATAA
</pre>

<p align="center">
<img src="https://user-images.githubusercontent.com/98036665/151746598-19682dd3-16ae-4470-91c8-03b9d1dadf38.png" width="600" height="600">
</p>

<pre>
Control sequence, translated 5'-3'
MFINRWLFSTNHKDIGTLYLLFGAWAGMVGTALSILIRAELGQPGALLGDDQIYNVIVTAHAFVMIFFMVMPMMIGGFGNWLVPLMIGAPDMAFPRMNNMSFWLLPPSF
LLLLASSMVEAGAGTGWTVYPPLAGNLAHAGASVDLTIFSLHLAGVSSILGAINFITTIINMKPPAMTQYQTPLFVWSVLITAVLLLLSLPVLAAGITMLLTDRNLNTT
FFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGIISHVVTYYSGKKEPFGYMGMVWAMMSIGFLGFIVWA<b>HH</b>MFTVGLDVDTRAYFTSATMIIAIPTGVKVFSWLATL
HGGNIKWSPAMLWALGFIFLFTVGGLTGIVLSNSSLDIVLHDTYYVVAHFHYVLSMGAVFAIMAGFVHWFPLFSGFTLDDTWAKAHFAIMFVGVNMTFFPQHFLGLSGM
PRRYSDYPDAYTTWNTVSSMGSFISLTAVLIMIFMIWEAFASKREVMSVSYASTNLEWLHGCPPPYHTFEEPTYVKVK-
</pre>

### Insertion of P291
Single amino acid insertion mutations occur when one codon or three contiguous nucleotides are inserted into a gene, resulting in an additional residue somewhere in the primary sequence of the translated protein. Below, a proline (P) residue was inserted between two wild-type histidines (H290, H291) involved in the the active site. I chose to insert a proline here because it is a notorious "helix breaker" due to its unique structure as the only amino acid where the side chain is connected to the backbone twice, making it rigid and disruptive to the regular α helical backbone conformation (Li, S C et al.). Analysis of the active site using [SWISS-MODEL](https://swissmodel.expasy.org/interactive/) confirms that P291 kinks the chain enough to pull H290 away from interacting in the active site; however, stabilization from other residues in the new active site (H240, V243, H292) still permits ligand binding, represented as two red dots on the model below.

<p align="center">
<img src="https://user-images.githubusercontent.com/98036665/151746314-195cc2b9-c203-438f-a71f-044387dc412d.png" width="600" height="600">
</p>

[BLASTp](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) alignment of the amino acid sequence below yielded the highest alignment score of 1023 with identities of 514/515(99%) compared to the control COX1 translated sequence.

<pre>
Mutated sequence, insertion of P291, translated 5'-3'
MFINRWLFSTNHKDIGTLYLLFGAWAGMVGTALSILIRAELGQPGALLGDDQIYNVIVTAHAFVMIFFMVMPMMIGGFGNWLVPLMIGAPDMAFPRMNNMSFWLLPPSF
LLLLASSMVEAGAGTGWTVYPPLAGNLAHAGASVDLTIFSLHLAGVSSILGAINFITTIINMKPPAMTQYQTPLFVWSVLITAVLLLLSLPVLAAGITMLLTDRNLNTT
FFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGIISHVVTYYSGKKEPFGYMGMVWAMMSIGFLGFIVWA<b>HPH</b>MFTVGLDVDTRAYFTSATMIIAIPTGVKVFSWLAT
LHGGNIKWSPAMLWALGFIFLFTVGGLTGIVLSNSSLDIVLHDTYYVVAHFHYVLSMGAVFAIMAGFVHWFPLFSGFTLDDTWAKAHFAIMFVGVNMTFFPQHFLGLSG
MPRRYSDYPDAYTTWNTVSSMGSFISLTAVLIMIFMIWEAFASKREVMSVSYASTNLEWLHGCPPPYHTFEEPTYVKVK-
</pre>

The figure below compares two dot plots. The control dot plot compared the unmutated sequence plotted against itself to create a baseline. The second dot plot shows the unmutated control COX1 nucleotide sequence on the X-axis and the COX1 sequence mutated with a proline codon insertion on the Y-axis. Flipping back and forth between the two plots reveals differences in their nucleotide sequence around the 850bp mark. The mutated graph expands, or gains nucelotides, compared to the control due to the insertion mutation. This change is most easily observed by watching the domain and range of the mutated graph expand compared to the control.

<p align="center">
<img src="https://user-images.githubusercontent.com/98036665/152700666-b02b93f7-0b49-43ee-a0fb-c5ef9c9a38da.gif" width="60%" height="60%"/>
</p>

Nucleotide alignment created using [BLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) produced the highest alignment score of 2777 with an identity of 1545/1548(99%) compared to the control COX1 gene. As shown below, there is a gap (-) of 3 nucleotides, caused by insertion of of *CCC*, encoding the proline, in the subject sequence after base pair 870. This matches the insertion shown in the dot plot above.

<p align="center">
<img width="653" alt="INS" src="https://user-images.githubusercontent.com/98036665/152700341-fe8f4617-6e03-4c4d-a476-8d3f9df519bb.png">
</p>

### Deletion of H291
Single amino acid deletion mutations occur when one codon or three contiguous nucleotides are removed from a gene, resulting in loss of a residue somewhere in the primary sequence of the translated protein. Below, H291 was removed from the wild-type primary sequence, resulting in a new active site with only three residues (H240, V243, H290). I chose to remove one of the contiguous active site histidines in order to determine if their combined activity is a significant factor in ligand association. Despite the deletion producing a similar active site to the insertion mutation, [SWISS-MODEL](https://swissmodel.expasy.org/interactive/) analysis of the active site shows that the new active site does not bind the ligand. This difference could be because removing the histidine entirely as opposed to displacing it from the active site also removes the peripheral stabilization of the second histidine.

<p align="center">
<img src="https://user-images.githubusercontent.com/98036665/151746630-abef7de7-909d-4c22-a074-f1af242938d7.png" width="600" height="600">
</p>

[BLASTp](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) alignment of the amino acid sequence below yielded the lowest alignment score of 1020 with identities of 513/514(99%) compared to the control COX1 translated sequence.

<pre>
Mutated sequence, deletion of H291, translated 5'-3'
MFINRWLFSTNHKDIGTLYLLFGAWAGMVGTALSILIRAELGQPGALLGDDQIYNVIVTAHAFVMIFFMVMPMMIGGFGNWLVPLMIGAPDMAFPRMNNMSFWLLPPSF
LLLLASSMVEAGAGTGWTVYPPLAGNLAHAGASVDLTIFSLHLAGVSSILGAINFITTIINMKPPAMTQYQTPLFVWSVLITAVLLLLSLPVLAAGITMLLTDRNLNTT
FFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGIISHVVTYYSGKKEPFGYMGMVWAMMSIGFLGFIVWA<b>H</b>MFTVGLDVDTRAYFTSATMIIAIPTGVKVFSWLATLH
GGNIKWSPAMLWALGFIFLFTVGGLTGIVLSNSSLDIVLHDTYYVVAHFHYVLSMGAVFAIMAGFVHWFPLFSGFTLDDTWAKAHFAIMFVGVNMTFFPQHFLGLSGMP
RRYSDYPDAYTTWNTVSSMGSFISLTAVLIMIFMIWEAFASKREVMSVSYASTNLEWLHGCPPPYHTFEEPTYVKVK-
</pre>

As before, a baseline was created by plotting the unmutated control sequence against itself. The second dot plot shows the unmutated control COX1 nucleotide sequence on the X-axis and the COX1 sequence mutated by the H291 codon deletion on the Y-axis. Like the insertion mutation, flipping back and forth between the two plots reveals differences in their nucleotide sequence around the 850bp mark; however, in this instance, the mutated graph contracts, or loses nucleotides, compared to the control due to the deletion mutation. This change is most easily observed by watching the domain and range of the mutated graph shrink compared to the control.

<p align="center">
<img src="https://user-images.githubusercontent.com/98036665/152242715-0fda2d53-6f93-4ab9-a6eb-91a48313cdc8.gif" width="60%" height="60%"/>
</p>


Nucleotide alignment created using [BLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) produced an alignment score of 2772 with an identity of 1545/1548(99%) compared to the control COX1 gene. As shown below, there is a gap (-) of 3 nucleotides, caused by deletion of *CAC*, encoding histidine at amino acid position 291, in the subject sequence starting at base pair 868. This matches the deletion shown in the dot plot above.

<p align="center">
<img width="653" alt="DEL" src="https://user-images.githubusercontent.com/98036665/152701058-ed3fd6b9-17ad-4c9b-a4f0-2337635a928f.png">
</p>

### Insertion of G290, G291 and Deletion of H290, H291
Single amino acid insertion and deletion mutations (INDELS) occur when one codon or three contiguous nucleotides are removed from a gene and different ones are inserted in their place, resulting in a substitution of an amino acid somewhere in the primary sequence of the translated protein. Below, H290 and H291 of the active site are removed and replaced with two glycine residues (G290, G291). I chose to insert glycine, the amino acid with the smallest R group, a single hydrogen, here to observe how completely replacing of the electrostatic effects of the contiguous active site histidines with small, hydrophobic residues would affect binding of a hydrophilic ligand. As shown by analysis of the active site using [SWISS-MODEL](https://swissmodel.expasy.org/interactive/), this substitution results in complete loss of active site integrity and a loss of ligand binding affinity.

<p align="center">
<img src="https://user-images.githubusercontent.com/98036665/151746674-c42968c7-bfd3-4401-9c9b-53f9acd134b4.png" width="600" height="600">
</p>

[BLASTp](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) alignment of the amino acid sequence below produced an alignment score of 1021 with identities of 512/514(99%) compared to the control COX1 translated sequence.

<pre>
Mutated sequence, deletion of H290 H291, insertion of G290 G291, translated 5'-3'
MFINRWLFSTNHKDIGTLYLLFGAWAGMVGTALSILIRAELGQPGALLGDDQIYNVIVTAHAFVMIFFMVMPMMIGGFGNWLVPLMIGAPDMAFPRMNNMSFWLLPPSF
LLLLASSMVEAGAGTGWTVYPPLAGNLAHAGASVDLTIFSLHLAGVSSILGAINFITTIINMKPPAMTQYQTPLFVWSVLITAVLLLLSLPVLAAGITMLLTDRNLNTT
FFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGIISHVVTYYSGKKEPFGYMGMVWAMMSIGFLGFIVWA<b>GG</b>MFTVGLDVDTRAYFTSATMIIAIPTGVKVFSWLATL
HGGNIKWSPAMLWALGFIFLFTVGGLTGIVLSNSSLDIVLHDTYYVVAHFHYVLSMGAVFAIMAGFVHWFPLFSGFTLDDTWAKAHFAIMFVGVNMTFFPQHFLGLSGM
PRRYSDYPDAYTTWNTVSSMGSFISLTAVLIMIFMIWEAFASKREVMSVSYASTNLEWLHGCPPPYHTFEEPTYVKVK-
</pre>

As before, a baseline was created by plotting the unmutated control sequence against itself. The second dot plot shows the unmutated control COX1 nucleotide sequence on the X-axis and the COX1 sequence mutated by a H290, H291 deletion and a G290, G291 insertion on the Y-axis. Similar to both previous mutations, flipping back and forth between the two plots reveals differences in their nucleotide sequence around the 850bp mark. Contrasting the prior mutations, we observe no change in the domain or range of the mutated graph compared to the control because INDEL mutations result in a replacement of nucleotides without affecting sequence length. Compared to the control, this change is most easily observed by watching how the appearance of the dots around the 850bp mark change between the two graphs.

<p align="center">
<img src="https://user-images.githubusercontent.com/98036665/152242496-e104637c-0838-4d8d-8201-2b8d03089255.gif" width="60%" height="60%"/>
</p>

Nucleotide alignment created using [BLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) produced the lowest alignment score of 2760 with an identity of 1539/1545(99%) compared to the control COX1 gene. As shown below, there is a mismatch of of 6 nucleotides, caused by deletion of *CACCAC*, encoding two histidines at amino acid positions 290 and 291, and insertion of *GGGGGG*, encoding two glycines at amino acid positions 290 and 291, in the subject sequence starting at base pair 868. This matches the INDEL mutation shown in the dot plot above.

<p align="center">
<img width="653" alt="INDEL" src="https://user-images.githubusercontent.com/98036665/152701414-e457b9f1-e679-40d5-8ae2-8357b2061c9c.png">
</p>

### Integrative Genomics Viewer Representation
Below is a snapshot generated in [IGV](https://software.broadinstitute.org/software/igv/download) of the three mutated sequences compared to each other and the control using a BLAST-like alignment tool (BLAT). This visual representation confirms the presence of the mutations in the specified regions.

![igv-mutations](https://user-images.githubusercontent.com/98036665/152702859-02b9d887-fa20-476f-ace1-0ab1d44ad265.png)

## Results & Discussion
The data gathered in this experiment supports the hypothesis that if there is mutation in the active site of an enzyme, the enzyme's affinity for substrate will decrease or cause loss of function. Though the insertion of P291 still allows ligand binding, SWISS-MODEL active site analysis shows decreased substrate affinity due to an alpha helix kink that prevents H290's adequate participation in active site stabilization. Likewise, the deletion of H291 eliminates ligand binding due to the removal of an integral active site residue. Alongside, the INDEL mutation, caused by deletion of H290, H291 and insertion of G290, G291 in its place, eliminates recognition of the active site by replacing two integral positively charged residues with hydrophobic glycine residues that lack electrostatic properties. As a result, the insertion mutation produced a relatively synonymous protein, and the deletion and INDEL mutations produced non-synonymous, malfunctioning proteins.

In relation to function, the loss of the activity of the catalytic site of COX1 is detrimental not only to subunit 1, but also the function of respiratory complex IV. In regular respiration and safe acceptance of terminal electrons, cytochrome c delivers electrons to COX one at a time, where they are stored until they can be passed four at a time to molecular oxygen, producing water as a product. As the site of terminal electron acceptance, loss of binding activity in the COX1 active site converts the electron transport chain from an ATP synthesis mechanism to a machine that readily produces reactive oxygen species that contribute to cellular damage and aging. Likewise, inefficient electron transport results in decreased mitochondrial output that affects growth at the organismal level. Often, COX1 mutations that eliminate function result in fetal reabsorption because the embryo lacks the energy it needs to grow (Vondrackova, Alzbeta et al.). On the other hand, COX1 mutations that decrease mitochondrial efficiency rather than eliminating it altogether can result in mitochondrial disorders that will eventually result in lethality, such as [Leigh syndrome](https://pubmed.ncbi.nlm.nih.gov/30743023/) in humans, which has been linked to mutations in mitochondria-encoded COX1 or its nuclear-encoded assembly factor, SURF1.

## Conclusions
The deletion mutation and the INDEL mutation both result in active site remodeling and loss of ligand binding. In a model organism such as *Mus musculus*, these mutations would likely result in embryonic reabsorption, as they achieve complete loss of COX1 function. The insertion mutation results in decreased mitochondrial function and could be an interesting mutation to observe in a *Mus musculus* model to study mitochondrial diseases that have decreased COX activity.

## References
1. Hunt, Ryan et al. “Silent (synonymous) SNPs: should we care about them?.” Methods in molecular biology (Clifton, N.J.) vol. 578 (2009): 23-39. doi:10.1007/978-1-60327-411-1_2
2. Li, S C et al. “Alpha-helical, but not beta-sheet, propensity of proline is determined by peptide environment.” Proceedings of the National Academy of Sciences of the United States of America vol. 93,13 (1996): 6676-81. doi:10.1073/pnas.93.13.6676
3. Vondrackova, Alzbeta et al. “High-resolution melting analysis of 15 genes in 60 patients with cytochrome-c oxidase deficiency.” Journal of human genetics vol. 57,7 (2012): 442-8. doi:10.1038/jhg.2012.49

[Return to Table of Contents](https://github.com/bmb0082/lab-project/blob/main/README.md#table-of-contents)

# Analyzing SNP transition and transversion ratios in three *Mus musculus* mitochondrial genes: cytochrome c oxidase subunit I, cytochrome b subunit, and the D-loop control region
Transitions (Ts) are SNPs that result in substitution between two purines or two pyrimidines, producing in a binding affinity that contains the same number of hydrogen bonds between bases as the wild type. Transversion (Tv) mutations result from substitution between a purine and pyrimidine, creating a variant strand that differs in number of hydrogen bonds from its complimentary strand. Generally, transversions are more likely to alter the amino acid sequence of proteins than transitions due to bucking of the DNA strand that alters interactions with transcription factors.

Due to heteroplasmy of the mitochondira, mtDNA transition/transversion (Ts/Tv) ratios can fall within a wide range. This data can be used to investigate the biological bias that seems to favor transition SNPs over tranversion SNPs (Ts/Tv > 1) due to the increased chance of transversions to detrimentally alter protein structure and function. Therefore, local deviations in the Ts/Tv ratio can be indicative of evolutionary selection of genes.

## Purpose - To determine the the ratio of transition to transversion mutations in two sample mitochondrial genes, COX1 and Cyt-b
**Hypothesis - The Ts/Tv ratios of the genes will be greater than 1, favoring milder transition mutations over more drastic transversion mutations**

## Methodology
A [workflow](https://usegalaxy.org/u/bmb002/w/snp-calling-by-gene--vcf-generation) was created in [Galaxy](https://usegalaxy.org/), an online platform used for data analysis and bioinformatics, outlining the procedure followed for SNP calling by gene and VCF data file generation.

The mitochondrial genome of brain tissue from *Mus musculus* was sequenced and uploaded to Galaxy as mutilple reads in FASTQ formatted files. Quality control using the FastQC program was performed for quality assurance of the sequence reads. Each of the read files were mapped to a the built-in mm10 reference genome using BWA-MEM, then merged to a single BAM file using the MergeSamFiles function. The whole genome BAM was constricted down to a select region containing the COX1 gene using the Slice tool set to restrict to nucleotide coordinates 5328..6872, which were obtained from an [NCBI reference genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_005089.1?report=fasta). Using the bcfTools mpileup program, a Variant Called Format (VCF) file was created from the COX1-restricted BAM in order to perform SNP and INDEL calling. Lastly, statisitcal analysis using the bcfTools Stats program was performed to generate an overview of the specific point mutation SNPs recognized in the dataset. The same procedure was followed to generate a Cyt-b dataset restricted to the mitochondrial coordinates 14145..15288, and a D-loop dataset restricted to the mitochondrial coordinates 1..877.

### Cytochrome c oxidase subunit I (COX1)
COX1 is the mitochondrially-encoded subunit 1 of cytochrome c oxidase, also known as respiratory complex IV. It participates in reduction of peroxide ions to water in terminal electron acceptance. Using the SNP frequencies determined in the VCF, the Ts/Tv ratio for the COX1-restricted dataset was calculated in Excel and graphed using RStudio. The resulting data is shown below.

| **Transitions (Ts)** |  |  |
| --- | --- | --- |
| A-G  | 16.3%  | 96 |
| C-T  | 16.8%  | 96 |
| **Transversions (Tv)** |  |  |
| A-T  | 24.2%  | 143 |
| A-C  | 15.9%  |  94 |
| G-T  | 22.7%  | 134 |
| G-C  | 4.0%  | 24 |
| **Ts/Tv ratio** | | **0.49** |

![COX1-SNPs-by-type](https://user-images.githubusercontent.com/98036665/159791580-14065c71-73d8-428b-8521-8fae109ac02e.png)

### Cytochrome b subunit (Cyt-b)
Cyt-b is the only mitochondrially encoded subunit of the cytochrome bc1 complex, also known as respiratory complex III in the electron transport chain. Its two heme groups participate in electron bifurcation in the lower half of the Q cycle. Using the SNP frequencies determined in the Cyt-b VCF, the Ts/Tv ratio for the Cyt-b-restricted dataset was calculated in Excel and graphed using RStudio. The resulting data is shown below.

| **Transitions (Ts)** |  |  |
| --- | --- | --- |
| A-G  | 13.2%  | 58 |
| C-T  | 17.6%  | 77 |
| **Transversions (Tv)** |  |  |
| A-T  | 28.3%  | 124 |
| A-C  | 18.3%  |  80 |
| G-T  | 20.5%  | 90 |
| G-C  | 2.1%  | 9 |
| **Ts/Tv ratio** | | **0.45** |

![cytb-SNPs-by-type](https://user-images.githubusercontent.com/98036665/159791608-2e523b46-f611-4add-a2d9-cd5422f40fda.png)

### D-loop control region
In the mitochondria, the The D-loop is not a gene, but a non-coding region that acts as a promoter for both the heavy and light strands of the mtDNA. It contains essential transcription and replication elements. Using the SNP frequencies determined in the D-loop VCF, the Ts/Tv ratio for the D-loop-restricted dataset was calculated in Excel and graphed using RStudio. The resulting data is shown below.

| **Transitions (Ts)** |  |  |
| --- | --- | --- |
| A-G | 16.7% | 48 |
| C-T  | 17.0%  | 49 |
| **Transversions (Tv)** |  |  |
| A-T  | 24.7%  | 71 |
| A-C  | 18.4%  |  53 |
| G-T  | 18.8%  | 54 |
| G-C  | 4.5%  | 13 |
| **Ts/Tv ratio** | | **0.51** |

![D-loop-SNPs-by-type](https://user-images.githubusercontent.com/98036665/159810003-af8fd516-b13c-4efc-8eac-a1cbf9e47da1.png)

## Results and Discussion
The data gathered in this experiment does not support the prediction that the Ts/Tv ratio would be greater than 1. The COX1-specific dataset shows an abnormal Ts/Tv ratio of 0.49, which is less than 1, indicating that transversions are favored over transitions at a ratio that is almost exactly equal to 2:1. Looking at the specific nucelotides involved, A-T (24.2%) transversions were the most frequent in the COX1-restricted data, and G-C (4.0%) transversions were the least frequent. Likewise, the Cyt-b-restricted dataset showed an abnormal Ts/Tv of 0.45, indicating an even greater favor of transversions over transitions. Looking at the specific nucelotides involved, A-T (28.3%) transversions were the most frequent in the Cyt-b-restricted data, and G-C (2.1%) transversions were the least frequent. Additionally, the D-loop-specific dataset showed an abnormal Ts/Tv ratio of 0.51. Similar to COX1, this indicates that transversions are favored over transitions at a ratio that is almost exactly equal to 2:1. Looking at the specific nucelotides involved, A-T (24.7%) transversions were the most frequent in the COX1-restricted data, and G-C (4.5%) transversions were the least frequent.

Given the conserved, important function of both cytochrome c oxidase subunit 1 and cytochrome b in their respective respiratory complexes, the unexpected bias for transversions raises many evolutionary and bioinformatics questions. What specific mutuations are present in the datasets and where are they in the consensus sequences? Why are these mutations tolerated and how, if at all, do they affect protein structure and function? Why are A-T transversions so frequent, and by what mechanisms or origins do these mutations arise? Why is the general biologcal bias for transitions not upheld in this dataset, instead favoring transversions that are statistically more likely to cause non-synonymous proteins? How can this dataset be further organized to reveal information about statistical significance?

These questions are subject to further investigation, possibly by exploring these SNPs grouped by their origins. Creation of COX1, Cyt-b, and D-loop consensus sequences may be helpful to elucidate specific SNPs for further investigation of how they may alter protein stucture and function.

## Conclusions
Bioinformatics analysis of a sample from *Mus musculus* brain tissue produced a Ts/Tv ratio of 0.49 for the mitochondrially-encoded gene COX1, 0.45 for mitochondrially-encoded Cyt-b, and 0.51 for the mitochondrial D-loop promoter region.

[Return to Table of Contents](https://github.com/bmb0082/lab-project/blob/main/README.md#table-of-contents)

# Investigating the origins of SNPs by type in three *Mus musculus* mitochondrial genes: cytochrome c oxidase subunit I, cytochrome b subunit, and the D-loop control region
In general, mitochondrial DNA has a 10x greater rate of mutation than nuclear DNA due to less efficient DNA repair, a more mutagenic environment, and higher load of replications per cell division.

The nuclear genome has five main mechanisms of DNA repair mechanisms: base excision repair (BER), nucleotide excision repair (NER), mismatch repair (MMR), homologous recombination (HR) and non-homologous end joining (NHEJ). Mitochondrial BER machinery shares some overlap with that of the nucleus, but there are still not a clear consensus regarding mitochondrial double-strand break repair. Current research shows that a majority of mitochondrial DNA repair is not repair at all, and instead comes from genome sharing that occurs when an injured mitochondria fuses with a healthy one. This heteroplasmy of mitochondrial DNA variants contributes greatly to the increased mutation frequency observed in mtDNA.

Additionally, the more mutagenic local environment of the mitochondria is directly linked to the formation of oxidative radicals (ROS) in the electron transport system. Electron transport in the inner membrane is a complex system that was selected for not because it is the most safe and efficient method, but because it works well enough for the organism to survive. There are many ways that this complex system can go wrong and result in production of ROS in the matrix that can easily interact with nearby mtDNA. Notably, C>A and G>T transversions originate from the ROS-mediated oxidation of guanine (G) to 8-oxoguainine (G`), which pairs with adenine (A) instead of cytosine (C) in the first round of replication. In the second round of replication, adenine will correctly pair with thymine, resulting in a GC>TA transversion. Therefore, any detected C>A and G>T SNPs in the dataset are determined to be ROS-mediated in origin.

![Mutagenic-ROS-med-mtDNA-damage](https://user-images.githubusercontent.com/98036665/160877320-1be47e9c-6d82-495a-8426-09c1a0969d8f.png)

Lastly, the mitochondria has an increased number of required replications per cell division than the nuclear genome. As a result, increased replication speed often increases the instances of DNA polymerase read errors that present genotypically as transition SNPs.

A combination of all these factors contributes to an increased mutation rate in the mitochondrial genome. The main origins of mutation that will be investigated in this experiment are ROS-mediated transversions and polymerase read error-meditated transitions.

## Purpose - To determine the statistical relationships between ROS-mediated SNPs, polymerase read error SNPs, and other SNP origins
**Hypothesis - ROS-mediated SNPs are the most frequent in mitochondrial genes**

## Methodology
Using the data generated in the sliced VCF for each gene from project 3, a purine/pyrimidine origins chart was created in Excel to determine the origins of each SNP type. Per mutagenic ROS-mediated mtDNA damage, C>A and G>T transversions originate from ROS-mediated oxidation. Therefore, any detected C>A and G>T SNPS in the dataset were determined to be ROS-mediated in origin. Transition SNPS, or mutations from purine to purine or pyrimidine to pyrimidine, are due to mitochondrial DNA polymerase gamma read errors. Therefore, any detected A<>G, C<>T SNPs in the dataset were determined to be DNA polymerase read error-mediated in origin.

For each gene, the purine/pyrimidine origin frequences for ROS, read error, and other mutations were plotted using RStudio and a linear regression was run to determine if there is a statistically significant difference in the means of each group. This data is shown below.

### Cytochrome c oxidase subunit I (COX1)

| SNP | Count | TsTv | Type | Origin |
| --- | --- | --- |  --- | --- |
| A>C | 40 | TV | Pur>Pyr | Other Origin |
| A>G | 58 | Ts | Pur>Pur | Read Error |
| A>T | 49 | Tv | Pur>Pyr | Other Origin |
| C>A | 54 | Tv | Pyr>Pur | ROS |
| C>G | 9 | Tv | Pyr>Pur | Other Origin |
| C>T | 46 | Ts | Pyr>Pyr | Read Error |
| G>A | 38 | Ts | Pur>Pur | Read Error |
| G>C | 15 | Tv | Pur>Pyr | Other Origin |
| G>T | 89 | Tv | Pur>Pyr | ROS |
| T>A | 95 | Tv | Pyr>Pur | Other Origin |
| T>C | 53 | Ts | Pyr>Pyr | Read Error |
| T>G | 45 | Tv | Pyr>Pur | Other Origin |

![ROS-COX1](https://user-images.githubusercontent.com/98036665/161463337-75658a1a-82d1-4225-a026-f09c9eed25b2.png)

| Origin | Count | Percent |
| --- | --- | --- |
| Read Error | 195 | 33.0% |
| Other | 253 | 42.8% |
| ROS | 143 | 24.0% |
| Total | 591 |  |

Linear regression of the above COX1 data with the ROS group set as the origin of reference produced a p-value of 0.3842, meaning that despite the visual difference in frequencies represented graphically, there is not a statistically significant difference between the mean frequencies of each origin of mutation.

### Cytochrome b subunit (Cyt-b)

| SNP | Count | TsTv | Type | Origin |
| --- | --- | --- |  --- | --- |
| A>C | 35 | TV | Pur>Pyr | Other Origin |
| A>G | 37 | Ts | Pur>Pur | Read Error |
| A>T | 50 | Tv | Pur>Pyr | Other Origin |
| C>A | 45 | Tv | Pyr>Pur | ROS |
| C>G | 5 | Tv | Pyr>Pur | Other Origin |
| C>T | 36 | Ts | Pyr>Pyr | Read Error |
| G>A | 21 | Ts | Pur>Pur | Read Error |
| G>C | 4 | Tv | Pur>Pyr | Other Origin |
| G>T | 62 | Tv | Pur>Pyr | ROS |
| T>A | 74 | Tv | Pyr>Pur | Other Origin |
| T>C | 41 | Ts | Pyr>Pyr | Read Error |
| T>G | 28 | Tv | Pyr>Pur | Other Origin |

![ROS-CytB](https://user-images.githubusercontent.com/98036665/161463349-fea3ca6e-bc8a-4f06-aa13-1823c3b59040.png)

| Origin | Count | Percent |
| --- | --- | --- |
| Read Error | 135 | 30.8% |
| Other | 196 | 44.8% |
| ROS | 107 | 24.4% |
| Total | 438 |  |

Linear regression of the above Cyt-b data with the ROS group set as the origin of reference produced a p-value of 0.4855, meaning that despite the visual difference in frequencies represented graphically, there is not a statistically significant difference between the mean frequencies of each origin of mutation.

### D-loop control region

| SNP | Count | TsTv | Type | Origin |
| --- | --- | --- |  --- | --- |
| A>C | 26 | TV | Pur>Pyr | Other Origin |
| A>G | 29 | Ts | Pur>Pur | Read Error |
| A>T | 22 | Tv | Pur>Pyr | Other Origin |
| C>A | 27 | Tv | Pyr>Pur | ROS |
| C>G | 7 | Tv | Pyr>Pur | Other Origin |
| C>T | 30 | Ts | Pyr>Pyr | Read Error |
| G>A | 19 | Ts | Pur>Pur | Read Error |
| G>C | 6 | Tv | Pur>Pyr | Other Origin |
| G>T | 45 | Tv | Pur>Pyr | ROS |
| T>A | 49 | Tv | Pyr>Pur | Other Origin |
| T>C | 19 | Ts | Pyr>Pyr | Read Error |
| T>G | 9 | Tv | Pyr>Pur | Other Origin |

![ROS-D-loop](https://user-images.githubusercontent.com/98036665/161464449-f46a1315-d28f-4e4f-b80a-5a73d41ebeb3.png)

| Origin | Count | Percent | 288
| --- | --- | --- |
| Read Error | 97 | 33.7% |
| Other | 119 | 41.3% |
| ROS | 72 | 25.0% |
| Total | 288 |  |

Linear regression of the above D-loop data with the ROS group set as the origin of reference produced a p-value of 0.3808, meaning that despite the visual difference in frequencies represented graphically, there is not a statistically significant difference between the mean frequencies of each origin of mutation.

## Results and Discussion
Though differences in the frequencies of mutation origins were observed, linear regression using RStudio produced statistically insignificant p-values for each gene. However, seperating SNPs by their origins and observing them across seperate genes reveals trends about the general population of mutations. In each gene, the ROS-mediated mutations are the least frequent, which was an unexpected outcome. Through further inspection, this observation has biological relevance. Production of ROS is expected in the mitochondria due to the dangerous nature of electron transport, but ROS-mediated mutations are largely circumvented due to the many mechanisms in place that protect the mtDNA. For example, mitochondria have a multilayer network of anti-oxidant systems, such as superoxide dismutase, catalase, glutathione peroxidase, and glutathione reductase, that protect mtDNA from oxidative damage. Though higher levels of ROS are expected in the mitochondria compared to other subcellular localizations, this does not always translate to increased ROS-mediated mutagenesis due to the subsequent increased level of protective mechanisms found in the mitochondria.

In each gene, the "other" category of transversion origins is the most frequent. This raises further questions that could be of biological relevance. Is this group composed of many small origins that add up to a large number? Or is something big missing in this picture? What are other mechanisms by which transversions can arise, and how can we elucidate our understanding? These questions are subject to further investigation and research into other origins of mutagenesis, both spontaneous and induced.

## Conclusions
There is no statistically significant difference between the frequencies of ROS-mediated mutations, read error-mediated mutations, and other types of mutations. Though many mutations were identified in the aligned read data, this number is drastically reduced during the creation of a consensus sequence. 

[Return to Table of Contents](https://github.com/bmb0082/lab-project/blob/main/README.md#table-of-contents)

# Modeling the effects of SNPs from *Mus musculus* consensus read data in three mitochondrial genes: cytochrome c oxidase subunit I, cytochrome b subunit, and the D-loop control region

## Purpose - To determine through sequence analysis how identified SNPs in collected consensus read data effect protein structure and function
**Hypothesis - The SNPs identified in the colony mice read data will have little to no effect on protein structure and function**

## Methodology
Though there were many mutations identified in the variant called format (VCF) read data, this number is drastically reduced when a consensus sequence is formed. Consensus sequence data shows us the mutations that most accurately characterize the mice in our colony at a given time by condensing reads down to the most frequent nucleotide observed at each position. In this experiment, consensus sequences for each gene (COX1, Cyt-b) and the D-loop were creating using the UseGalaxy IVAR consensus tool from their respective coordinate-restriced BAMs from project 3. Then, BLASTn was run on each sample, comparing the consensus reads to the mm10 reference genome from NCBI in order to identify SNPs in the consensus sequence. These SNPs were then classified by origin. The consensus sequences were translated to amino acid sequences using the [ExPASy Translate Tool](https://web.expasy.org/translate/). Finally, BLASTp was performed on the comparing the translated sequences to a the translated references in order to determine if the SNPs resuled in mutations in the primary amino acid sequence. This protocol was followed for COX1, Cyt-b, and the D-loop control region datasets.

### Cytochrome c oxidase subunit I (COX1)
COX1 is the mitochondrially-encoded subunit 1 of cytochrome c oxidase, also known as respiratory complex IV. COX is a multimeric enzyme with intricate and highly regulated assembly that involves multiple cofactors and associated assmebly proteins. It contains 13 subunits, where the three catalytic subunits (I-III) are mitochondrial coded, and the other 10 are nuclear coded. As a unit, COX is a directly redox-linked proton pump that facilitates the reduction of peroxide ions to water in terminal electron acceptance. Due the essential function of COX1 and its binding of heme-a and heme-a3-CuB, COX1's active sites have conserved structure and function across eukaryotic taxa. However, other regions of the gene are more variable. The MT-CO1 gene is often used as a DNA barcode to identify animal species because its sequence is conserved among conspecifics, but the mutation rate of non-conserved areas is often fast enough to distinguish closely related species.

```
Reference genome source of mt-Co1 gene, presented in the 5'-3' direction:
Mus musculus mitochondrion, complete genome
NCBI Reference Sequence: NC_005089.1

GenBank Graphics
>NC_005089.1:5328-6872 Mus musculus mitochondrion, complete genome
ATGTTCATTAATCGTTGATTATTCTCAACCAATCACAAAGATATCGGAACCCTCTATCTACTATTCGGAG
CCTGAGCGGGAATAGTGGGTACTGCACTAAGTATTTTAATTCGAGCAGAATTAGGTCAACCAGGTGCACT
TTTAGGAGATGACCAAATTTACAATGTTATCGTAACTGCCCATGCTTTTGTTATAATTTTCTTCATAGTA
ATACCAATAATAATTGGAGGCTTTGGAAACTGACTTGTCCCACTAATAATCGGAGCCCCAGATATAGCAT
TCCCACGAATAAATAATATAAGTTTTTGACTCCTACCACCATCATTTCTCCTTCTCCTAGCATCATCAAT
AGTAGAAGCAGGAGCAGGAACAGGATGAACAGTCTACCCACCTCTAGCCGGAAATCTAGCCCATGCAGGA
GCATCAGTAGACCTAACAATTTTCTCCCTTCATTTAGCTGGAGTGTCATCTATTTTAGGTGCAATTAATT
TTATTACCACTATTATCAACATGAAACCCCCAGCCATAACACAGTATCAAACTCCACTATTTGTCTGATC
CGTACTTATTACAGCCGTACTGCTCCTATTATCACTACCAGTGCTAGCCGCAGGCATTACTATACTACTA
ACAGACCGCAACCTAAACACAACTTTCTTTGATCCCGCTGGAGGAGGGGACCCAATTCTCTACCAGCATC
TGTTCTGATTCTTTGGGCACCCAGAAGTTTATATTCTTATCCTCCCAGGATTTGGAATTATTTCACATGT
AGTTACTTACTACTCCGGAAAAAAAGAACCTTTCGGCTATATAGGAATAGTATGAGCAATAATGTCTATT
GGCTTTCTAGGCTTTATTGTATGAGCCCACCACATATTCACAGTAGGATTAGATGTAGACACACGAGCTT
ACTTTACATCAGCCACTATAATTATCGCAATTCCTACCGGTGTCAAAGTATTTAGCTGACTTGCAACCCT
ACACGGAGGTAATATTAAATGATCTCCAGCTATACTATGAGCCTTAGGCTTTATTTTCTTATTTACAGTT
GGTGGTCTAACCGGAATTGTTTTATCCAACTCATCCCTTGACATCGTGCTTCACGATACATACTATGTAG
TAGCCCATTTCCACTATGTTCTATCAATGGGAGCAGTGTTTGCTATCATAGCAGGATTTGTTCACTGATT
CCCATTATTTTCAGGCTTCACCCTAGATGACACATGAGCAAAAGCCCACTTCGCCATCATATTCGTAGGA
GTAAACATAACATTCTTCCCTCAACATTTCCTGGGCCTTTCAGGAATACCACGACGCTACTCAGACTACC
CAGATGCTTACACCACATGAAACACTGTCTCTTCTATAGGATCATTTATTTCACTAACAGCTGTTCTCAT
CATGATCTTTATAATTTGAGAGGCCTTTGCTTCAAAACGAGAAGTAATATCAGTATCGTATGCTTCAACA
AATTTAGAATGACTTCATGGCTGCCCTCCACCATATCACACATTCGAGGAACCAACCTATGTAAAAGTAA
AATAA
```

```
Colony sequence read data of mt-Co1, presented in the 5'-3' direction:
Sliced IVAR consensus read, restricted to nucleotide coordinates 5328..6872

ATGTTCATTAATCGTTGATTATTCTCAACCAATCACAAAGATATCGGAACCCTCTATCTACTATTCGGAG
CCTGAGCGGGAATAGTGGGTACTGCACTAAGTATTTTAATTCGAGCAGAATTAGGTCAACCAGGTGCACT
TTTAGGAGATGACCAAATTTACAATGTTATCGTAACTGCCCATGCTTTTGTTATAATTTTCTTCATAGTA
ATACCAATAATAATCGGAGGCTTTGGAAACTGACTTGTCCCACTAATAATCGGAGCCCCAGATATAGCAT
TCCCACGAATAAATAATATAAGTTTTTGACTCCTACCACCATCATTTCTCCTTCTCCTAGCATCATCAAT
AGTAGAAGCAGGAGCAGGAACAGGATGAACAGTCTACCCACCTCTAGCCGGAAATCTAGCCCATGCAGGA
GCATCAGTAGACCTAACAATTTTCTCCCTTCATTTAGCTGGAGTGTCATCTATTTTAGGTGCAATTAATT
TTATTACCACTATTATCAACATGAAGCCCCCAGCCATAACACAGTATCAAACTCCACTATTTGTCTGATC
CGTACTTATTACAGCCGTACTGCTCCTATTATCACTACCAGTACTAGCCGCAGGCATTACTATACTACTA
ACAGACCGCAACCTAAACACAACTTTCTTTGATCCCGCTGGAGGAGGGGACCCAATTCTCTACCAGCATC
TGTTCTGATTCTTTGGACACCCAGAAGTTTATATTCTTATCCTCCCAGGATTTGGAATTATTTCACATGT
AGTTACTTACTACTCCGGAAAAAAAGAACCTTTCGGCTATATAGGAATAGTATGAGCAATAATGTCTATT
GGCTTTCTAGGCTTTATTGTATGAGCCCACCACATATTCACAGTAGGATTAGATGTAGACACACGAGCTT
ACTTTACATCAGCCACTATAATTATCGCAATTCCTACCGGTGTCAAAGTATTTAGCTGACTTGCAACCCT
ACACGGAGGTAATATTAAATGATCTCCAGCTATACTATGAGCCTTAGGCTTTATTTTCTTATTTACAGTT
GGTGGTCTAACCGGAATTGTTTTATCCAACTCATCCCTTGACATCGTGCTTCACGATACATACTATGTAG
TAGCCCATTTCCACTATGTTCTATCAATGGGAGCAGTGTTTGCTATCATAGCAGGATTTGTTCACTGATT
CCCATTATTTTCAGGCTTCACCCTAGATGACACATGAGCAAAAGCCCACTTCGCCATCATATTCGTAGGA
GTAAACATAACATTCTTCCCTCAACATTTCCTGGGCCTTTCAGGAATACCACGACGCTACTCAGACTACC
CAGATGCTTACACCACATGAAACACTGTCTCTTCTATAGGATCATTTATTTCACTAACAGCTGTTCTCAT
CATGATCTTTATAATTTGAGAGGCCTTTGCTTCAAAACGAGAAGTAATATCAGTATCGTATGCTTCAACA
AATTTAGAATGACTTCATGGCTGCCCTCCACCATATCACACATTCGAGGAACCAACCTATGTAAAAGTAA
AATAA
```

The results of the BLASTp in COX1 identified 4 mismatch SNPs: 225T>C, 516A>G, 603G>A, 717G>A. As shown in the table below, all of the SNPs identified in COX1 are transitions due to polymerase read error.

| Position | Mutation | TsTv | Type | Origin |
| --- | --- | --- |  --- | --- |
| 225 | T>C | Ts | Pyr>Pyr | Read Error |
| 516 | A>G | Ts | Pur>Pur | Read Error |
| 603 | G>A | Ts | Pur>Pur | Read Error |
| 717 | G>A | Ts | Pur>Pur | Read Error |

Next, translated amino acid sequences were generated using [ExPASy Translate Tool](https://web.expasy.org/translate/) for the control COX1 and restricted read data, as shown below

```
COX1 control, open reading frame 1 sequence, translated 5'-3':
MFINRWLFSTNHKDIGTLYLLFGAWAGMVGTALSILIRAELGQPGALLGDDQIYNVIVTAHAFVMIFFMVM
PMMIGGFGNWLVPLMIGAPDMAFPRMNNMSFWLLPPSFLLLLASSMVEAGAGTGWTVYPPLAGNLAHAGAS
VDLTIFSLHLAGVSSILGAINFITTIINMKPPAMTQYQTPLFVWSVLITAVLLLLSLPVLAAGITMLLTDR
NLNTTFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGIISHVVTYYSGKKEPFGYMGMVWAMMSIGFLG
FIVWAHHMFTVGLDVDTRAYFTSATMIIAIPTGVKVFSWLATLHGGNIKWSPAMLWALGFIFLFTVGGLTG
IVLSNSSLDIVLHDTYYVVAHFHYVLSMGAVFAIMAGFVHWFPLFSGFTLDDTWAKAHFAIMFVGVNMTFF
PQHFLGLSGMPRRYSDYPDAYTTWNTVSSMGSFISLTAVLIMIFMIWEAFASKREVMSVSYASTNLEWLHG
CPPPYHTFEEPTYVKVK-
```

```
COX1 read data, open reading frame 1 sequence, translated 5'-3':
MFINRWLFSTNHKDIGTLYLLFGAWAGMVGTALSILIRAELGQPGALLGDDQIYNVIVTAHAFVMIFFMVM
PMMIGGFGNWLVPLMIGAPDMAFPRMNNMSFWLLPPSFLLLLASSMVEAGAGTGWTVYPPLAGNLAHAGAS
VDLTIFSLHLAGVSSILGAINFITTIINMKPPAMTQYQTPLFVWSVLITAVLLLLSLPVLAAGITMLLTDR
NLNTTFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGIISHVVTYYSGKKEPFGYMGMVWAMMSIGFLG
FIVWAHHMFTVGLDVDTRAYFTSATMIIAIPTGVKVFSWLATLHGGNIKWSPAMLWALGFIFLFTVGGLTG
IVLSNSSLDIVLHDTYYVVAHFHYVLSMGAVFAIMAGFVHWFPLFSGFTLDDTWAKAHFAIMFVGVNMTFF
PQHFLGLSGMPRRYSDYPDAYTTWNTVSSMGSFISLTAVLIMIFMIWEAFASKREVMSVSYASTNLEWLHG
CPPPYHTFEEPTYVKVK-
```

Finally, BLASTp was performed in order to determine if the identified SNPs result in downstream mutations in primary amino acid structure, and therefore protein folding and function. Shown below is the basic amino acid local alignment of the control (query) vs the read sequence data (subject) for COX1.

<img width="717" alt="COX1-BLASTp" src="https://user-images.githubusercontent.com/98036665/162488655-681aa5fb-b55d-4a15-a436-a1deb3edafe0.png">
 
BLASTp of colony consensus COX1 vs the control shows 100% amino acid sequence identity despite the presence of 4 mismatch SNPs, thus maintaining the same protein structure and function. A model of COX1 made with [SWISS-MODEL](https://swissmodel.expasy.org/interactive) is shown below.

<img src="https://user-images.githubusercontent.com/98036665/162526503-a4b54859-cb8a-453c-a2bc-eec34396b655.png" width="600" height="600">

### Cytochrome b subunit (Cyt-b)
Cyt-b is the only mitochondrially encoded subunit of the cytochrome c reductase complex, also known as respiratory complex III in the electron transport chain. Its two heme groups, Cyt-bL and Cyt-bH, participate in electron bifurcation in the lower half of the Q cycle. It contains 11 subunits, cytochrome b, cytochrome c1, the Rieske protein, plus two core proteins and 6 low molecular weight proteins. Interestingly, the mitochondrial Cyt-b gene is ideal for phylogenetic studies and species identification, as it shows limited variability within species, and much greater variation between species.

```
Reference genome source of mt-Cyt-b gene, presented in the 5'-3' direction:
Mus musculus mitochondrion, complete genome
NCBI Reference Sequence: NC_005089.1

GenBank Graphics
>NC_005089.1:5328-6872 Mus musculus mitochondrion, complete genome
ATGACAAACATACGAAAAACACACCCATTATTTAAAATTATTAACCACTCATTCATTGACCTACCTGCCC
CATCCAACATTTCATCATGATGAAACTTTGGGTCCCTTCTAGGAGTCTGCCTAATAGTCCAAATCATTAC
AGGTCTTTTCTTAGCCATACACTACACATCAGATACAATAACAGCCTTTTCATCAGTAACACACATTTGT
CGAGACGTAAATTACGGGTGACTAATCCGATATATACACGCAAACGGAGCCTCAATATTTTTTATTTGCT
TATTCCTTCATGTCGGACGAGGCTTATATTATGGATCATATACATTTATAGAAACCTGAAACATTGGAGT
ACTTCTACTGTTCGCAGTCATAGCCACAGCATTTATAGGCTACGTCCTTCCATGAGGACAAATATCATTC
TGAGGTGCCACAGTTATTACAAACCTCCTATCAGCCATCCCATATATTGGAACAACCCTAGTCGAATGAA
TTTGAGGGGGCTTCTCAGTAGACAAAGCCACCTTGACCCGATTCTTCGCTTTCCACTTCATCTTACCATT
TATTATCGCGGCCCTAGCAATCGTTCACCTCCTCTTCCTCCACGAAACAGGATCAAACAACCCAACAGGA
TTAAACTCAGATGCAGATAAAATTCCATTTCACCCCTACTATACAATCAAAGATATCCTAGGTATCCTAA
TCATATTCTTAATTCTCATAACCCTAGTATTATTTTTCCCAGACATACTAGGAGACCCAGACAACTACAT
ACCAGCTAATCCACTAAACACCCCACCCCATATTAAACCCGAATGATATTTCCTATTTGCATACGCCATT
CTACGCTCAATCCCCAATAAACTAGGAGGTGTCCTAGCCTTAATCTTATCTATCCTAATTTTAGCCCTAA
TACCTTTCCTTCATACCTCAAAGCAACGAAGCCTAATATTCCGCCCAATCACACAAATTTTGTACTGAAT
CCTAGTAGCCAACCTACTTATCTTAACCTGAATTGGGGGCCAACCAGTAGAACACCCATTTATTATCATT
GGCCAACTAGCCTCCATCTCATACTTCTCAATCATCTTAATTCTTATACCAATCTCAGGAATTATCGAAG
ACAAAATACTAAAATTATATCCAT
```

```
Colony sequence read data of mt-Cyt-b, presented in the 5'-3' direction:
Sliced IVAR consensus read, restricted to nucleotide coordinates 14145..15288

ATGACAAACATACGAAAAACACACCCATTATTTAAAATTATTAACCACTCATTCATTGACCTACCTGCCC
CATCCAACATTTCATCATGATGAAACTTTGGGTCCCTTCTAGGAGTCTGCCTAATAGTCCAAATCATTAC
AGGTCTTTTCTTAGCCATACACTACACATCAGATACAATAACAGCCTTTTCATCAGTAACACACATTTGT
CGAGACGTAAATTACGGGTGACTAATCCGATATATACACGCAAACGGAGCCTCAATATTTTTTATTTGCT
TATTCCTTCATGTCGGACGAGGCTTATATTATGGATCATATACATTTATAGAAACCTGAAACATTGGAGT
ACTTCTACTGTTCGCAGTCATAGCCACAGCATTTATAGGCTACGTCCTTCCATGAGGACAAATATCATTC
TGAGGTGCCACAGTTATTACAAACCTCCTATCAGCCATCCCATATATTGGAACAACCCTAGTCGAATGAA
TTTGAGGGGGCTTCTCAGTAGACAAAGCCACCTTAACCCGATTCTTCGCTTTCCACTTCATCTTACCATT
TATTATCGCGGCCCTAGCAATCGTTCACCTCCTTTTCCTCCACGAAACAGGATCAAACAACCCAACAGGA
TTAAACTCAGATGCAGATAAAATTCCATTTCACCCCTACTATACAATCAAAGATATCCTAGGTATCCTAA
TCATATTCTTAATTCTCATAACCCTAGTATTATTTTTCCCAGACATACTAGGAGACCCAGACAACTACAT
ACCAGCTAATCCACTAAACACCCCACCCCATATTAAACCCGAATGATATTTCCTATTTGCATACGCCATT
CTACGCTCAATCCCCAATAAACTAGGAGGTGTCCTAGCCTTAATCTTATCTATCCTAATTTTAGCCCTAA
TACCTTTCCTTCATACCTCAAAGCAACGAAGCCTAATATTCCGCCCAATCACACAAATTTTGTACTGAAT
CCTAGTAGCCAACCTACTTATCTTAACCTGAATTGGGGGCCAACCAGTAGAACACCCATTTATTATCATT
GGCCAACTAGCCTCCATCTCATACTTCTCAATCATCTTAATTCTTATACCAATCTCAGGAATTATCGAAG
ACAAAATACTAAAATTATATCCAT
```

Basic nucleotide local alignment of the above control sequence (query) vs the consenus sequence (subject) was performed using BLASTn to locate SNPs in the consensus read data, shown below.

<img width="368" alt="Cyt-b-BLASTn" src="https://user-images.githubusercontent.com/98036665/162497882-caaba962-dfcf-4b51-a01c-4f1fd9ad9184.png">

The results of the BLASTp in Cyt-b identified 2 mismatch SNPs: 525G>A, 594C>T. As shown in the table below, all of the SNPs identified in Cyt-b are transitions due to polymerase read error.

| Position | Mutation | TsTv | Type | Origin |
| --- | --- | --- |  --- | --- |
| 525 | G>A | Ts | Pyr>Pyr | Read Error |
| 594 | C>T | Ts | Pyr>Pyr | Read Error |

Next, translated amino acid sequences were generated using [ExPASy Translate Tool](https://web.expasy.org/translate/) for the control Cyt-b and restricted read data, as shown below

```
Cyt-b control, open reading frame 1 sequence, translated 5'-3':
MTNMRKTHPLFKIINHSFIDLPAPSNISSWWNFGSLLGVCLMVQIITGLFLAMHYTSDTMTAFSSVTHIC
RDVNYGWLIRYMHANGASMFFICLFLHVGRGLYYGSYTFMETWNIGVLLLFAVMATAFMGYVLPWGQMSF
WGATVITNLLSAIPYIGTTLVEWIWGGFSVDKATLTRFFAFHFILPFIIAALAIVHLLFLHETGSNNPTG
LNSDADKIPFHPYYTIKDILGILIMFLILMTLVLFFPDMLGDPDNYMPANPLNTPPHIKPEWYFLFAYAI
LRSIPNKLGGVLALILSILILALMPFLHTSKQRSLMFRPITQILYWILVANLLILTWIGGQPVEHPFIII
GQLASISYFSIILILMPISGIIEDKMLKLYP
```

```
Cyt-b read data, open reading frame 1 sequence, translated 5'-3':
MTNMRKTHPLFKIINHSFIDLPAPSNISSWWNFGSLLGVCLMVQIITGLFLAMHYTSDTMTAFSSVTHIC
RDVNYGWLIRYMHANGASMFFICLFLHVGRGLYYGSYTFMETWNIGVLLLFAVMATAFMGYVLPWGQMSF
WGATVITNLLSAIPYIGTTLVEWIWGGFSVDKATLTRFFAFHFILPFIIAALAIVHLLFLHETGSNNPTG
LNSDADKIPFHPYYTIKDILGILIMFLILMTLVLFFPDMLGDPDNYMPANPLNTPPHIKPEWYFLFAYAI
LRSIPNKLGGVLALILSILILALMPFLHTSKQRSLMFRPITQILYWILVANLLILTWIGGQPVEHPFIII
GQLASISYFSIILILMPISGIIEDKMLKLYP
```

Finally, BLASTp was performed in order to determine if the identified SNPs result in downstream mutations in primary amino acid structure, and therefore protein folding and function. Shown below is the basic amino acid local alignment of the control (query) vs the read sequence data (subject) for Cyt-b.

<img width="717" alt="Cyt-b-BLASTp" src="https://user-images.githubusercontent.com/98036665/162488789-e9dff98d-1916-4ef7-8a88-b1722f948b9e.png">
 
BLASTp of colony consensus Cyt-b vs the control shows 100% amino acid sequence identity despite the presence of 2 mismatch SNPs, thus maintaining the same protein structure and function. A model of Cyt-b made with [SWISS-MODEL](https://swissmodel.expasy.org/interactive) is shown below.

<img src="https://user-images.githubusercontent.com/98036665/162528144-c16b8697-cfa4-4847-9c30-07bfa43dcf78.png" width="600" height="600">

### D-loop control region
The D-loop is technically not a gene, but rather, a non-coding region of mtDNA that acts as a promoter for the heavy strand of the circular mitochondrial genome. As such, it contains essential transcription and replication elements. Mutations in this region may serve as potential indicators of cellular DNA damage, and according to the recent literature, increased rates of D-loop mutation have also been linked to the pathogenesis of many types of cancer.

```
Reference genome source of mt-D-loop, presented in the 5'-3' direction:
Mus musculus mitochondrion, complete genome
NCBI Reference Sequence: NC_005089.1

GenBank Graphics
>NC_005089.1:5328-6872 Mus musculus mitochondrion, complete genome
GTTAATGTAGCTTAATAACAAAGCAAAGCACTGAAAATGCTTAGATGGATAATTGTATCCCATAAACACA
AAGGTTTGGTCCTGGCCTTATAATTAATTAGAGGTAAAATTACACATGCAAACCTCCATAGACCGGTGTA
AAATCCCTTAAACATTTACTTAAAATTTAAGGAGAGGGTATCAAGCACATTAAAATAGCTTAAGACACCT
TGCCTAGCCACACCCCCACGGGACTCAGCAGTGATAAATATTAAGCAATAAACGAAAGTTTGACTAAGTT
ATACCTCTTAGGGTTGGTAAATTTCGTGCCAGCCACCGCGGTCATACGATTAACCCAAACTAATTATCTT
CGGCGTAAAACGTGTCAACTATAAATAAATAAATAGAATTAAAATCCAACTTATATGTGAAAATTCATTG
TTAGGACCTAAACTCAATAACGAAAGTAATTCTAGTCATTTATAATACACGACAGCTAAGACCCAAACTG
GGATTAGATACCCCACTATGCTTAGCCATAAACCTAAATAATTAAATTTAACAAAACTATTTGCCAGAGA
ACTACTAGCCATAGCTTAAAACTCAAAGGACTTGGCGGTACTTTATATCCATCTAGAGGAGCCTGTTCTA
TAATCGATAAACCCCGCTCTACCTCACCATCTCTTGCTAATTCAGCCTATATACCGCCATCTTCAGCAAA
CCCTAAAAAGGTATTAAAGTAAGCAAAAGAATCAAACATAAAAACGTTAGGTCAAGGTGTAGCCAATGAA
ATGGGAAGAAATGGGCTACATTTTCTTATAAAAGAACATTACTATACCCTTTATGAAACTAAAGGACTAA
GGAGGATTTAGTAGTAAATTAAGAATAGAGAGCTTAA
```

```
Colony sequence read data of mt-D-loop, presented in the 5'-3' direction:
Sliced IVAR consensus read, restricted to nucleotide coordinates 1..877

GTTAATGTAGCTTAATAACAAAGCAAAGCACTGAAAATGCTTAGATGGATAGTTATATCCCATAAACACA
AAGGTTTGGTCCTGGCCTTATAATTAATTAGAGGTAAAATTACACATGCAAACCTCCATAGACCGGTGTA
AAATCCCTTAAACATTTACTTAAAATTTAAGGAGAGGGTATCAAGCACATTAAAATAGCTTAAGACACCT
TGCCTAGCCACACCCCCACGGGACTCAGCAGTGATAAATATTAAGCAATAAACGAAAGTTTGACTAAGTT
ATACCTCTTAGGGTTGGTAAATTTCGTGCCAGCCACCGCGGTCATACGATTAACCCAAACTAATTATCTT
CGGCGTAAAACGTGTCAACTATAAATAAATTAATAGAATTAAAATCCAACTTATATGTGAAAATTCATTG
TTAGGACCTAAACTCAATAACGAAAGTAATTCTAATCATTTATAATACACGACAGCTAAGACCCAAACTG
GGATTAGATACCCCACTATGCTTAGCCATAAACCTAAATAATTAAATTTAACAAAACTATTTGCCAGAGA
ACTACTAGCCATAGCTTAAAACTCAAAGGACTTGGCGGTACTTTATATCCATCTAGAGGAGCCTGTTCTA
TAATCGATAAACCCCGCTCTACCTCACCATCTCTTGCTAATTCAGCCTATATACCGCCATCTTCAGCAAA
CCCTAAAAAGGTATTAAAGTAAGCAAAAGAATCAAACATAAAAACGTTAGGTCAAGGTGTAGCCAATGAA
ATGGGAAGAAATGGGCTACATTTTCTTATAAAAGAACATTACTATACCCTTTATGAAACTAAAGGACTAA
GGAGGATTTAGTAGTAAATTAAGAATAGAGAGCTTAATTGAATTGAGCAATGAAGTACGCACACACCGCC
CGTCACCCTCCTCAAATTAAATTAAACTTAACATAATTAATTTCTAGACATCCGTTTATGAGAGGAGATA
AGTCGTAACAAGGTAAGCATACTGGAAAGTGTGCTTGGAATAATCMT
```

Basic nucleotide local alignment of the above control sequence (query) vs the consenus sequence (subject) was performed using BLASTn to locate SNPs in the consensus read data, shown below.

The results of the BLASTp in the D-loop identified 4 mismatch SNPs: 52A>G, 55G>A, 381A>T, 465G>A. As shown in the table below, 3 of the 4 of the SNPs identified in the D-loop are transitions due to polymerase read error. The remaining SNP is a transversion with an unknown origin.

| Position | Mutation | TsTv | Type | Origin |
| --- | --- | --- |  --- | --- |
| 52 | A>G | Ts | Pur>Pur | Read Error |
| 55 | G>A | Ts | Pur>Pur | Read Error |
| 381 | A>T | Tv | Pur>Pyr | Other Origin |
| 465 | G>A | Ts | Pur>Pur | Read Error |

Since the D-loop is a regulatory region that is not transcribed or translated, BLASTp and SWISS-MODEL, which are protein bioinformatics tools, cannot be used. Therefore, only BLASTn was performed to identify SNPs. Predictions regarding change in the function the D-loop and the H strand promoter it contains are based on whether or not the number of mutations present affect polymerase binding for transcription or mitochondrial genome replication.

## Results and Discussion
The 4 SNPs located in COX1 (225T>C, 516A>G, 603G>A, 717G>A) were all transition mutations that originated from polymerase read error that did not produce a change in amino acid sequence compared to the control reference sequence. Similarly, the 2 SNPs located in Cyt-b (525G>A, 594C>T) were also transitions due to read error, and they did not result in a change the downstream amino acid sequence. In both cases, the mutations were located in the 3rd position of their respective codons; therefore, they were able to maintain the same primary sequence because, in many cases, only the first two positions of a codon are read due to the degeneracy of the genetic code. This is especially true in the mitochondia, where tRNAs are limited. In addition, the 4 SNPs identified in the D-loop (52A>G, 55G>A, 381A>T, 465G>A) are unlikely to affect DNA and RNA polymerase binding. In general, point mutations, both transitions and transversions, are not as detrimental to D-loop function as insertion or deletion mutations.

## Conclusions
The SNPs found in COX1 and Cyt-b do not have a downstream effect in protein structure to the wobble base. Additionally, the low number of mutations in the D-loop are unlikely to cause dysfunctional polymerase binding, but do not exclude change in function if more mutations were to accumulate.
