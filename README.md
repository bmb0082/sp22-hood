# Table of Contents

<details>
<summary>Projects</summary>

 1. [Modeling transcript directionality in *Mus musculus* mitochondrial cytochrome c oxidase subunit I](https://github.com/bmb0082/lab-project/blob/main/README.md#modeling-transcript-directionality-in-mus-musculus-mitochondrial-cytochrome-c-oxidase-subunit-i)
 2. [Modeling insertions, deletions, and INDEL mutations in *Mus musculus* mitochondrial cytochrome c oxidase subunit I](https://github.com/bmb0082/lab-project/blob/main/README.md#modeling-insertions-deletions-and-indel-mutations-in-mus-musculus-mitochondrial-cytochrome-c-oxidase-subunit-i)
 3. [Analyzing SNP transition and transversion ratios in *Mus musculus* cytochrome c oxidase subunit I](https://github.com/bmb0082/lab-project/blob/main/README.md#analyzing-snp-transition-and-transversion-ratios-in-mus-musculus-cytochrome-c-oxidase-subunit-i)
 
</details>


# Modeling transcript directionality in *Mus musculus* mitochondrial cytochrome c oxidase subunit I
Mitochrondrial cytochrome c oxidase I (COX1), encoded by the mt-Co1 gene, is a subunit of cytochrome c oxidase (RC-IV), the last enzyme in the mitochondiral electron transport chain. Correct native protein structure is essential for safe acceptance of terminal electons and overall mitochondiral efficiency. As a result, errors in COX1 can trap electrons in respiratory complex-III (RC-III), resulting in increased production of reactive oxygen species. Cytochrome c oxidase deficiency is the most common cause of genetic mitochondiral disorders, presenting as multi-organ, heterogeneous symptoms depending on the level of mitochondiral efficiency.

## Purpose - To understand the importance of direction when modeling proteins from translated sequences

Below is the selected reference for the mitochondiral gene encoding subunit I of cytochrome c oxidase. I chose this gene because COX1 contains the catalytic unit of cytochrome c oxidase; therefore, it's correct stucture and function are vital for the enzyme's overall activity. I wanted to see how differently the protein would fold depending on the directionality of the translated sequence.
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
Reading the translated sequence 5'-3' is the correct direction, as it generates a defined structure of cytochrome c oxidase I. [ExPASy Translate Tool](https://web.expasy.org/translate/) recognized the entire gene in this direction as a single open reading frame, and [SWISS-MODEL](https://swissmodel.expasy.org/interactive/) generated a 94.16% sequence homology to bovine heart cytochrome c oxidase subunit I. Likewise, it models a protein-ligand interaction pipeline between the active site and a peroxide ion, shown as two red dots in the catalytic core.

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
This is not the correct folding conformation of cytochrome c oxidase I. Reading the translated sequence 3'-5' results in a polypeptide that lacks structural integrity and the binding pockets for heme-a, heme-a3, and CuB copper (II) cofactors that are necessary for successful electron transport.

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

# Modeling insertions, deletions, and INDEL mutations in *Mus musculus* mitochondrial cytochrome c oxidase subunit I
Mutations are changes in DNA sequence that result in a variant form that can be passed down to subsequent generations. Most commonly, alterations in the structure of a gene are caused by single nucleotide transition or transversion mutations called single nucleotide polymorphisms (SNPs). However, mutations can also arise from insersion, deletion, or rearrangement of sections of a gene or chromosome.

Mutations do not always have to cause a change in phenotype. These mutations, called neutral mutations, produce *synonymous proteins* that have the same function as the wild type. This can occur when single nucleotide changes encode the same amino acid or a residue with similar enough electrostatic properties to produce a protein with unaltered function. Alternatively, mutations that produce *non-synonymous* proteins, or proteins with altered function, usually, but not always, result from more drastic changes in the genome. For example, frameshift mutations such as insertions or deletions or rearrangement mutations at the chromosomal level can result in a completely different downstream primary structures, causing improper folding and loss of function.

In this project, all mutations will be performed at the amino acid level, such as changes, insertions, or deletions of one or more amino acids (or multiples of 3 nucleotides aligned as codons) that may or may not affect the overall function of the protein.

## Purpose - To investigate the effects of amino acid mutations on predicted protein structure
Hypothesis - Mutations in the active site of an enzyme will decrease the enzyme's affinity for substrate or cause loss of function.

## Methodology
The wild-type gene encoding COX1 was translated in the 5'-3' direction using [ExPASy Translate Tool](https://web.expasy.org/translate/) and modeled based on homology using [SWISS-MODEL](https://swissmodel.expasy.org). Based on the model of the wild-type control COX1 protein, the residues involved in the active site were identified as H240, V243, H290, and H291. I chose to create mutations in these amino acids in order to observe how they would affect ligand binding, as the active site residues are the most indicative of protein function.

An insertion mutation, deletion mutation, and an insertion/deletion mutation respectively were added directly into the translated sequence and modeled as predicted proteins, producing representational models for all three in order to observe physical changes in protein structure. Dotplots comparing nucelotide sequences of each mutation against the control were created using [EMBOSS](https://www.bioinformatics.nl/cgi-bin/emboss/dotmatcher) in order to view the differences between the sequences through a visual representation of what occurs for each mutation and where it occurs in the nucleotide sequence. Alignment scores of the nucleotide and amino acid sequences were producing using [BLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) and [BLASTp](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) respectively in order to view sequence homology compared to the control.

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
Single amino acid insertion mutations occur when one codon or three contiguous nucleotides are inserted into a gene, resulting in an additional residue somewhere in the primary sequence of the translated protein. Below, a proline (P) residue was inserted between the two wild-type histidines (H290, H291) involved in the the active site. I chose to insert a proline here because it is a notorious "helix breaker" due to its unique structure as the only amino acid where the side chain is connected to the backbone twice, making it rigid and disruptive to the regular α helical backbone conformation. Analysis of the active site using [SWISS-MODEL](https://swissmodel.expasy.org/interactive/) confirms that P291 kinks the chain enough to pull H290 away from interacting in the active site; however, stabilization from other residues in the new active site (H240, V243, H292) still permits ligand binding, represented as two red dots on the model below.

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

The figure below compares two dot plots. The control dot plot compared the unmutated sequence plotted agaisnt itself to create a baseline. The second dot plot shows the unmutated control COX1 nucleotide sequence on the X-axis and the COX1 sequence mutated with a proline codon insertion on the Y-axis. Flipping back and forth between the two plots reveals differences in their nucleotide sequence around the 850bp mark. The mutated graph expands, or gains nucelotides, compared to the control due to the insertion mutation. This change is most easily observed by watching the domain and range of the mutated graph expand compared to the control.

<p align="center">
<img src="https://user-images.githubusercontent.com/98036665/152700666-b02b93f7-0b49-43ee-a0fb-c5ef9c9a38da.gif" width="60%" height="60%"/>
</p>

Nucleotide alignment created using [BLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) produced the highest alignment score of 2777 with an identity of 1545/1548(99%) compared to the control COX1 gene. As shown below, there is a gap (-) of 3 nucleotides, caused by insertion of of *CCC*, encoding the proline, in the subject sequence after base pair 870. This matches the insertion shown in the dot plot above.

<p align="center">
<img width="653" alt="INS" src="https://user-images.githubusercontent.com/98036665/152700341-fe8f4617-6e03-4c4d-a476-8d3f9df519bb.png">
</p>

### Deletion of H291
Single amino acid deletion mutations occur when one codon or three contiguous nucleotides are removed from a gene, resulting in loss of a residue somewhere in the primary sequence of the translated protein. Below, H291 was removed from the wild-type primary sequence, resuling in a new active site with only three residues (H240, V243, H290). I chose to remove one of the contiguous active site histidines in order to see if their combined activity is a significant factor in ligand association. Despite the deletion producing a similar active site to the insertion mutation, [SWISS-MODEL](https://swissmodel.expasy.org/interactive/) analysis of the active site shows that the new active site does not bind the ligand. This difference could be because removing the histidine entirely as opposed to displacing it from the active site also removes the peripheral stabilization of the second histidine.

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

As before, a baseline was created by plotting the unmutated control sequence agaisnt itself. The second dot plot shows the unmutated control COX1 nucleotide sequence on the X-axis and the COX1 sequence mutated by the H291 codon deletion on the Y-axis. Similar to the insertion mutation, flipping back and forth between the two plots reveals differences in their nucleotide sequence around the 850bp mark; however, in this instance, the mutated graph contracts, or loses nucelotides, compared to the control due to the deletion mutation. This change is most easily observed by watching the domain and range of the mutated graph shrink compared to the control.

<p align="center">
<img src="https://user-images.githubusercontent.com/98036665/152242715-0fda2d53-6f93-4ab9-a6eb-91a48313cdc8.gif" width="60%" height="60%"/>
</p>


Nucleotide alignment created using [BLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) produced an alignment score of 2772 with an identity of 1545/1548(99%) compared to the control COX1 gene. As shown below, there is a gap (-) of 3 nucleotides, caused by deletion of *CAC*, encoding histidine at amino acid position 291, in the subject sequence starting at base pair 868. This matches the deletion shown in the dot plot above.

<p align="center">
<img width="653" alt="DEL" src="https://user-images.githubusercontent.com/98036665/152701058-ed3fd6b9-17ad-4c9b-a4f0-2337635a928f.png">
</p>

### Insertion of G290, G291 and Deletion of H290, H291
Single amino acid insertion and deletion mutations (INDELS) occur when one codon or three contiguous nucleotides are removed from a gene and different ones are inserted in their place, resulting in a substitution of an amino acid somewhere in the primary sequence of the translated protein. Below, H290 and H291 of the active site are removed and replaced with two glycine residues (G290, G291). I chose to insert glycine, the amino acid with the smallest R group, a single hydrogen, here in order to observe how completely replacing of the electrostatic effects of the contiguous active site histidines with small, hydrophobic residues would affect binding of a hydrophilic ligand. As shown by analysis of the active site using [SWISS-MODEL](https://swissmodel.expasy.org/interactive/), this substitution results in complete loss of active site integrity and a loss of ligand binding affinity.

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

As before, a baseline was created by plotting the unmutated control sequence agaisnt itself. The second dot plot shows the unmutated control COX1 nucleotide sequence on the X-axis and the COX1 sequence mutated by a H290, H291 deletion and a G290, G291 insertion on the Y-axis. Similar to the both previous mutations, flipping back and forth between the two plots reveals differences in their nucleotide sequence around the 850bp mark. Contrasting the two previous mutations, there is no change in the domain or range of the mutated graph compared to the control because INDEL mutations result in a replacement of nucleotides rather than addition or subtraction. Compared to the control, this change is most easily observed by watching how the appearance of the dots around the 850bp mark change between the two graphs.

<p align="center">
<img src="https://user-images.githubusercontent.com/98036665/152242496-e104637c-0838-4d8d-8201-2b8d03089255.gif" width="60%" height="60%"/>
</p>

Nucleotide alignment created using [BLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) produced the lowest alignment score of 2760 with an identity of 1539/1545(99%) compared to the control COX1 gene. As shown below, there is a mismatch of of 6 nucleotides, caused by deletion of *CACCAC*, encoding two histidines at amino acid positions 290 and 291, and insertion of *GGGGGG*, encoding two glycines at amino acid positions 290 and 291, in the subject sequence starting at base pair 868. This matches the INDEL mutation shown in the dot plot above.

<p align="center">
<img width="653" alt="INDEL" src="https://user-images.githubusercontent.com/98036665/152701414-e457b9f1-e679-40d5-8ae2-8357b2061c9c.png">
</p>

### Integrative Genomics Viewer Representation
Below is a snapshot generated in [IGV](https://software.broadinstitute.org/software/igv/download) of a visual represenation of the three mutated sequences compared to each other and the control using a BLAST-like alignment tool (BLAT).

![igv-mutations](https://user-images.githubusercontent.com/98036665/152702859-02b9d887-fa20-476f-ace1-0ab1d44ad265.png)

## Results & Discussion
The data gathered in this experiment supports the hypothesis that if there is mutation in the active site of an enzyme, the enzyme's affinity for substrate will decrease or cause loss of function. Though the insertion of P291 still allows ligand binding, SWISS-MODEL active site analysis shows decreased substrate affinity due to an alpha helix kink that prevents H290's adequate participation in active site stabilization. Likewise, the deletion of H291 eliminates ligand binding due to the removal of an integral active site residue. Alongside, the INDEL mutation, caused by deletion of H290, H291 and insertion of G290, G291 in its place, eliminates recognition of the active site by replacing two integral positively charged residues with hydrophobic glycine residues that lack electrostatic properties. As a result, the insertion mutuation produced a relatively synonymous protein, and the deletion and INDEL mutations produced non-synonymous, malfunctioning proteins.

In relation to function, the loss of the activity of the catalytic site of COX1 is detrimental not only to subunit 1, but also the function of respiratory complex IV as a whole. In regular respiration and safe acceptance of terminal electrons, cytochrome c delivers electrons  to COX one at a time, where they are stored until they can be passed four at a time to molecular oxygen, producing water as a product. As the site of terminal electron acceptance, loss of binding activity in the COX1 active site converts the electron transport chain from an ATP synthesis mechanism to a machine that readily produces reactive oxygen species that contribute to cellular damage and aging. Likewise, inefficient electron transport results in decreased mitochondrial output that affects growth at the organismal level. Often times, COX1 mutations that eliminate function result in fetal reabsorption because the embryo lacks the energy it needs to grow. On the other hand, COX1 mutations that decrease mitochondrial efficiency rather than eliminating it altogether will often result in mitochondrial disorders that will eventually result in lethality, such as [Leigh syndrome](https://pubmed.ncbi.nlm.nih.gov/30743023/) in humans. 

## Conclusions
The deletion mutation and the INDEL mutation both result in active site remodeling and loss of ligand binding. In a model organism such as *Mus musculus*, these mutations would likely result in embryonic reabsorption, as they achieve complete loss of COX1 function. The insertion mutation results in decreased mitochondrial function and could be an interesting mutation to observe in studying mitochondiral disease in a *Mus musculus* model.

# Analyzing SNP transition and transversion ratios in *Mus musculus* cytochrome c oxidase subunit I
Introduction text

## Purpose - To investigate the the ratio of transition to transversion mutations in a COX1 nucleotide sequence and map specific SNPs to a COX1 consensus sequence

## Methodology
The mitochondrial genome of Mus musculus brain tissue was sequenced and uploaded to [Galaxy](https://usegalaxy.org/), an online workflow platform used for data analysis and bioinformatics, as mutilple FASTQ formatted files. Quality control using the FastQC program was performed for quality assurance of the sequence reads. Each of the read files were mapped to a the built-in mm10 reference genome using BWA-MEM, then merged to a single BAM file using the MergeSamFiles function. The whole genome BAM was constricted down to a select region containing the COX1 gene using the Slice tool set to restrict to nucleotide coordinates 5328..6872, which were obtained from an [NCBI reference genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_005089.1?report=fasta). A Variant Called Format (VCF) file was created from the COX1-restricted BAM in order to perform SNP and INDEL using the bcfTools mpileup program. Lastly, statisitcal analysis using bcfTools stats was performed to generate overview of specific point mutation SNPs.

A [workflow](https://usegalaxy.org/u/bmb002/w/snp-calling-by-gene--vcf-generation) was created in Galaxy outlining the procedure followed for SNP calling by gene and VCF data file generation.

## Conclusions

