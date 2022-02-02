# Project 1 - Modeling *Mus musculus* mitochondrial cytochrome c oxidase subunit I
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

### Predicted model of the mt-Co1 protein, presented 5'-3'
Reading the translated sequence 5'-3' is the correct direction, as it generates a defined structure of cytochrome c oxidase I. [ExPASy Translate Tool](https://web.expasy.org/translate/) recognized the entire gene in this direction as a single open reading frame, and [SWISS-MODEL](https://swissmodel.expasy.org/interactive/SAhTnr/models/) generated a 94.16% sequence homology to bovine heart cytochrome c oxidase subunit I. Likewise, it models a protein-ligand interaction pipeline between the active site and a peroxide ion, shown as two red dots in the catalytic core.

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

### Predicted model of the mt-Co1 protein, presented 3'-5'
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

# Project 2 - Modeling mutations of *Mus musculus* mitochondrial cytochrome c oxidase subunit I
Mutations are changes in DNA sequence that result in a variant form that can be passed down to subsequent generations. Most commonly, alterations in the structure of a gene are caused by single nucleotide transversion or transversion mutations called single nucleotide polymorphisms (SNPs). However, mutations can also arise from insersion, deletion, or rearrangement of sections of a gene or chromosome.

Mutations do not always have to cause a change in phenotype. These mutations, called neutral mutations, produce *synonymous proteins* that have the same function as the wild type. This can occur when single nucleotide changes encode the same amino acid or a residue with similar enough electrostatic properties to produce a protein with unaltered function. Alternatively, mutations that produce *non-synonymous* proteins, or proteins with altered function, usually, but not always, result from more drastic changes in the genome. For example, frameshift mutations such as insertions or deletions or rearrangement mutations at the chromosomal level can result in a completely different downstream primary structure, causing improper folding and complete loss of function.

In this project, all mutations will be performed at the amino acid level, such as changes, insertions, or deletions of one or more amino acids that may or may not affect the overall function of the protein.

## Purpose - To investigate the effects of amino acid mutations on predicted protein structure
Hypothesis - If there is mutation in the active site of an enzyme, the enzyme's affinity for substrate will decrease or cause loss of function.

## Methodology
The wild-type gene encoding COX1 was translated in the 5'-3' direction using [ExPASy Translate Tool](https://web.expasy.org/translate/) and modeled based on homology using [SWISS-MODEL](https://swissmodel.expasy.org). Based on the model of the wild-type control COX1 protein, the residues involved in the active site were identified as H240, V243, H290, and H291. I chose to create mutations in these amino acids in order to observe how they would affect ligand binding, as the active site residues are the most indicative of protein function.

An insertion mutation, deletion mutation, and an insertion/deletion mutation respectively were added directly into the translated sequence at the amino acid level and modeled, producing models for all three in order to observe physical changes in protein structure. Dotplots comparing nucelotide seqeucens of each mutation agaisnt the control were created using [JDotter](https://4virology.net/virology-ca-tools/jdotter/) in order to view the differences between the sequences and have a visual representation of what occurs for each mutation.

Below is the selected nucelotide reference sequence for the cytochrome c oxidase subunit I, its 5'-3' translated sequence, and the modeled protein prediction that will serve as the control for the mutated proteins.

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
GGCTTTCTAGGCTTTATTGTATGAGCC/CACCAC/ATATTCACAGTAGGATTAGATGTAGACACACGAGCTT
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

<p align="center">
<img src="https://user-images.githubusercontent.com/98036665/151746598-19682dd3-16ae-4470-91c8-03b9d1dadf38.png" width="600" height="600">
</p>

```
Control sequence, translated 5'-3'
MFINRWLFSTNHKDIGTLYLLFGAWAGMVGTALSILIRAELGQPGALLGDDQIYNVIVTAHAFVMIFFMVMPMMIGGFGNWLVPLMIGAPDMAFPRMNNMSFWLLPPSF
LLLLASSMVEAGAGTGWTVYPPLAGNLAHAGASVDLTIFSLHLAGVSSILGAINFITTIINMKPPAMTQYQTPLFVWSVLITAVLLLLSLPVLAAGITMLLTDRNLNTT
FFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGIISHVVTYYSGKKEPFGYMGMVWAMMSIGFLGFIVWA/HH/MFTVGLDVDTRAYFTSATMIIAIPTGVKVFSWLATL
HGGNIKWSPAMLWALGFIFLFTVGGLTGIVLSNSSLDIVLHDTYYVVAHFHYVLSMGAVFAIMAGFVHWFPLFSGFTLDDTWAKAHFAIMFVGVNMTFFPQHFLGLSGM
PRRYSDYPDAYTTWNTVSSMGSFISLTAVLIMIFMIWEAFASKREVMSVSYASTNLEWLHGCPPPYHTFEEPTYVKVK-
```

### Insertion Mutation
Single amino acid insertion mutations occur when one codon or three contiguous nucleotides are inserted into a gene, resulting in an additional residue somewhere in the primary sequence of the translated protein. Below, a proline (P) residue was inserted between the two wild-type histidines (H290, H291) involved in the the active site. I chose to insert a proline here because it is a notorious "helix breaker" due to its unique structure as the only amino acid where the side chain is connected to the backbone twice, making it rigid and disruptive to the regular α helical backbone conformation. Analysis of the active site using [SWISS-MODEL](https://swissmodel.expasy.org/interactive/6DHyNG/models/) confirms that P291 kinks the chain enough to pull H290 away from interacting in the active site; however, stabilization from other residues in the new active site (H240, V243, H292) still permits ligand binding, represented as two red dots on the model below.

<p align="center">
<img src="https://user-images.githubusercontent.com/98036665/151746314-195cc2b9-c203-438f-a71f-044387dc412d.png" width="600" height="600">
</p>

```
Mutated sequence, insertion of P291, translated 5'-3'
MFINRWLFSTNHKDIGTLYLLFGAWAGMVGTALSILIRAELGQPGALLGDDQIYNVIVTAHAFVMIFFMVMPMMIGGFGNWLVPLMIGAPDMAFPRMNNMSFWLLPPSF
LLLLASSMVEAGAGTGWTVYPPLAGNLAHAGASVDLTIFSLHLAGVSSILGAINFITTIINMKPPAMTQYQTPLFVWSVLITAVLLLLSLPVLAAGITMLLTDRNLNTT
FFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGIISHVVTYYSGKKEPFGYMGMVWAMMSIGFLGFIVWA/HPH/MFTVGLDVDTRAYFTSATMIIAIPTGVKVFSWLAT
LHGGNIKWSPAMLWALGFIFLFTVGGLTGIVLSNSSLDIVLHDTYYVVAHFHYVLSMGAVFAIMAGFVHWFPLFSGFTLDDTWAKAHFAIMFVGVNMTFFPQHFLGLSG
MPRRYSDYPDAYTTWNTVSSMGSFISLTAVLIMIFMIWEAFASKREVMSVSYASTNLEWLHGCPPPYHTFEEPTYVKVK-
```

Explanation of INS dot plot changes observed and why.

Add centered & resized INS nucelotide dot plot GIF here.

### Deletion Mutation
Single amino acid deletion mutations occur when one codon or three contiguous nucleotides are removed from a gene, resulting in loss of a residue somewhere in the primary sequence of the translated protein. Below, H291 was removed from the wild-type primary sequence, resuling in a new active site with only three residues (H240, V243, H290). I chose to remove one of the contiguous active site histidines in order to see if their combined activity is a significant factor in ligand association. Despite the deletion producing a similar active site to the insertion mutation, [SWISS-MODEL](https://swissmodel.expasy.org/interactive/MqvkuH/models/) analysis of the active site shows that the new active site does not bind the ligand. This difference could be because removing the histidine entirely as opposed to displacing it from the active site also removes the peripheral stabilization of the second histidine.

<p align="center">
<img src="https://user-images.githubusercontent.com/98036665/151746630-abef7de7-909d-4c22-a074-f1af242938d7.png" width="600" height="600">
</p>

```
Mutated sequence, deletion of H291, translated 5'-3'
MFINRWLFSTNHKDIGTLYLLFGAWAGMVGTALSILIRAELGQPGALLGDDQIYNVIVTAHAFVMIFFMVMPMMIGGFGNWLVPLMIGAPDMAFPRMNNMSFWLLPPSF
LLLLASSMVEAGAGTGWTVYPPLAGNLAHAGASVDLTIFSLHLAGVSSILGAINFITTIINMKPPAMTQYQTPLFVWSVLITAVLLLLSLPVLAAGITMLLTDRNLNTT
FFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGIISHVVTYYSGKKEPFGYMGMVWAMMSIGFLGFIVWA/H/MFTVGLDVDTRAYFTSATMIIAIPTGVKVFSWLATLH
GGNIKWSPAMLWALGFIFLFTVGGLTGIVLSNSSLDIVLHDTYYVVAHFHYVLSMGAVFAIMAGFVHWFPLFSGFTLDDTWAKAHFAIMFVGVNMTFFPQHFLGLSGMP
RRYSDYPDAYTTWNTVSSMGSFISLTAVLIMIFMIWEAFASKREVMSVSYASTNLEWLHGCPPPYHTFEEPTYVKVK-
```

Explanation of DEL dot plot changes observed and why.

Add centered & resized DEL nucelotide dot plot GIF here.

### Insertion and Deletion Mutation
Single amino acid insertion and deletion mutations (INDELS) occur when one codon or three contiguous nucleotides are removed from a gene and different ones are inserted in their place, resulting in a substitution of an amino acid somewhere in the primary sequence of the translated protein. Below, H290 and H291 of the active site are removed and replaced with two glycine residues (G290, G291). I chose to insert glycine, the amino acid with the smallest R group, a single hydrogen, here in order to observe how completely replacing of the electrostatic effects of the contiguous active site histidines with small, hydrophobic residues would affect binding of a hydrophilic ligand. As shown by analysis of the active site using [SWISS-MODEL](https://swissmodel.expasy.org/interactive/4XSm2G/models/), this substitution results in complete loss of active site integrity and a loss of ligand binding affinity.

<p align="center">
<img src="https://user-images.githubusercontent.com/98036665/151746674-c42968c7-bfd3-4401-9c9b-53f9acd134b4.png" width="600" height="600">
</p>

```
Mutated sequence, deletion of H290 and H291, insertion of G290 G291, translated 5'-3'
MFINRWLFSTNHKDIGTLYLLFGAWAGMVGTALSILIRAELGQPGALLGDDQIYNVIVTAHAFVMIFFMVMPMMIGGFGNWLVPLMIGAPDMAFPRMNNMSFWLLPPSF
LLLLASSMVEAGAGTGWTVYPPLAGNLAHAGASVDLTIFSLHLAGVSSILGAINFITTIINMKPPAMTQYQTPLFVWSVLITAVLLLLSLPVLAAGITMLLTDRNLNTT
FFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGIISHVVTYYSGKKEPFGYMGMVWAMMSIGFLGFIVWA/GG/MFTVGLDVDTRAYFTSATMIIAIPTGVKVFSWLATL
HGGNIKWSPAMLWALGFIFLFTVGGLTGIVLSNSSLDIVLHDTYYVVAHFHYVLSMGAVFAIMAGFVHWFPLFSGFTLDDTWAKAHFAIMFVGVNMTFFPQHFLGLSGM
PRRYSDYPDAYTTWNTVSSMGSFISLTAVLIMIFMIWEAFASKREVMSVSYASTNLEWLHGCPPPYHTFEEPTYVKVK-
```

Explanation of INDEL dot plot changes observed and why.

Add centered & resized INDEL nucelotide dot plot GIF here.

## Results
Discuss your findings. Do they support or reject your hypothesis? Why or why not?

## Conclusions
End the section by stating which mutation you think would achieve a change in function and why (remember what the protein does).
