# Modeling *Mus musculus* mitochondrial cytochrome c oxidase I
Mitochrondrial cytochrome c oxidase I (mt-Co1) is a subunit of cytochrome c oxidase (RC-IV), the last enzyme in the mitochondiral electron transport chain. Correct native protein structure is essential for safe acceptance of terminal electons and overall mitochondiral efficiency. As a result, errors in mt-Co1 can trap electrons in respiratory complex-III (RC-III), resulting in increased production of reactive oxygen species. Cytochrome c oxidase deficiency is the most common cause of genetic mitochondiral disorders, presenting as multi-organ, heterogeneous symptoms depending on the level of mitochondiral efficiency.

## **Purpose: To understand the importance of direction when modeling proteins from translated sequences**

Below is the selected mitochondiral gene encoding subunit I of for cytochrome c oxidase. I chose this gene because mt-Co1 contains the catalytic unit of cytochrome c oxidase; therefore, it's correct stucture and function are vital for the enzyme's overall activity. I wanted to see how different the protein would fold depending on the directionality of the translated sequence.
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
Reading the translated sequence 5'-3' is the correct direction, as it generates a defined structure of cytochrome c oxidase I. [ExPASy Translate Tool](https://web.expasy.org/translate/) recognized the entire gene in this direction as a single open reading frame, and [SWISS-MODEL](https://swissmodel.expasy.org/interactive/SAhTnr/models/) generated a 94.16% sequence homology to bovine heart cytochrome c oxidase subunit I. Likewise, it models a protein-ligand interaction pipeline between the active site and a peroxide ion.

<p align="center">
<img src="https://user-images.githubusercontent.com/98036665/151239026-dd14ff6f-7f56-4f55-8523-3b206c1305ee.png" width="600" height="600">
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

# Modeling mutations of *Mus musculus* mitochondrial cytochrome c oxidase I
Background information goes here

## **Purpose: To investigate the effect of mutations on predicted protein structure**
