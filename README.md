# Modeling *Mus musculus* mitochondrial cytochrome c oxidase I
Mitochrondrial cytochrome c oxidase I (mt-Co1) is a subunit of cytochrome c oxidase (RC-IV), the last enzyme in the mitochondiral electron transport chain. Correct native protein structure is essential for safe acceptance of terminal electons and overall mitochondiral efficiency. As a result, errors in mt-Co1 can trap electrons in respiratory complex-III (RC-III), resulting in increased production of reactive oxygen species.

## **Purpose: To understand the importance of direction when modelling proteins from translated sequences**

Below is selected mitochondiral gene encoding subunit I of for cytochrome c oxidase. I chose this gene because mutations that affect function of subunit I are some of the most common causes of genetic mitochondiral disorders. Subunit I of cytochome c oxidase is the home of the catalytic unit; therefore, it's correct stucture and function are vital for cytochrome c oxidase activity.
```
Source mt-Co1 gene, presented in the 5'-3' direction
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
Reading the translated sequence 5'-3' is the correct direction, as it generates a defined structure of cytochrome c oxidase I. ExPASy translate tool recognized the entire gene in this direction as a single open reading frame, and SWISS-MODEL generated a sequence identity value of 94.16% compared a known sequence for bovine heart cytochrome c oxidase subunit I. Likewise, it models a protein-ligand interaction pipeline between the active site and a peroxide ion.
![COX1 5-3](https://user-images.githubusercontent.com/98036665/151239026-dd14ff6f-7f56-4f55-8523-3b206c1305ee.png)
```
Open reading frame 1 translated sequence, 5'-3':
MFINRWLFSTNHKDIGTLYLLFGAWAGMVGTALSILIRAELGQPGALLGDDQIYNVIVTAHAFVMIFFMVMPMMIGGFGNWLVPLMIGAPDMAFPRMNNMSFWLLPPSF
LLLLASSMVEAGAGTGWTVYPPLAGNLAHAGASVDLTIFSLHLAGVSSILGAINFITTIINMKPPAMTQYQTPLFVWSVLITAVLLLLSLPVLAAGITMLLTDRNLNTT
FFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGIISHVVTYYSGKKEPFGYMGMVWAMMSIGFLGFIVWAHHMFTVGLDVDTRAYFTSATMIIAIPTGVKVFSWLATL
HGGNIKWSPAMLWALGFIFLFTVGGLTGIVLSNSSLDIVLHDTYYVVAHFHYVLSMGAVFAIMAGFVHWFPLFSGFTLDDTWAKAHFAIMFVGVNMTFFPQHFLGLSGM
PRRYSDYPDAYTTWNTVSSMGSFISLTAVLIMIFMIWEAFASKREVMSVSYASTNLEWLHGCPPPYHTFEEPTYVKVK-
```

### Predicted mt-Co1 protein, ORF 3'-5'
This is not the correct folding conformation of cytochrome c oxidase I. Reading the translated sequence 3'-5' results in a polypeptide that lacks structural integrity and the binding pockets for heme-a, heme-a3, and CuB copper (II) cofactors that are necessary for successful electron transport.
![COX1 3-5](https://user-images.githubusercontent.com/98036665/150719390-4b9598a1-a27a-4266-84ca-448126c0dbc9.png)

```
Reference genome:
Mus musculus mitochondrion, complete genome
NCBI Reference Sequence: NC_005089.1

GenBank Graphics
>NC_005089.1 Mus musculus mitochondrion, complete genome
```
