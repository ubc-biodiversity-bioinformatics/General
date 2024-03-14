# Common file formats


## Fasta files
- Assemblies (chromosome level/scaffold level), Protein files, gene annotations 
- extensions: .fasta, .fa, .cds.fa 

#### Format:
>Scaffold_name
Sequence
>Scaffold_name
Sequence

#### Notes:
- Some software require for there to be a character limit for each line of the sequence 
- Some sofrware require for there to be no special characters (eg *) to be run 
    - Eg. Interproscan: to replace * with X : 
    ```awk '{ gsub(/\*/, "X"); print }' Protein_file.fasta > protein_file_interproscan_input.fasta```
- To get number of scaffolds

    ```grep ">" file.fasta | wc -l ```


## GFF and GFF3 files
- Annotation output, a specific type of tsv file
- Extensions .gff3, .gff 

| GFF | GFF3 |
|-----| ----- |
| 0 or 1 based indices | 1 based indices|
| start and end order can depend | start < end | 
| ... | ...| 

#### Gff3 format is by column:
1. Sequence ID (has to match fasta )
2. Source  (Maker/Augustus/Genebank etc)
3. CDS/mRNA/gene/exon/intron etc - (http://www.sequenceontology.org/)
4. Start
5. End
6. Score (e values or p values)
7. Strand (+/-)
8. Phase (0/1/2)
9. Id, name, parent tags

#### Helpful notes on formatting:

#### Resources for gff and gff3 files
- https://plastid.readthedocs.io/en/latest/concepts/gff3.html
- https://learn.gencore.bio.nyu.edu/ngs-file-formats/gff3-format/

## BED files
- extensions: .bed, bedGraph

#### Format:
Chr Start End 

#### Resources 
https://bedtools.readthedocs.io/en/latest/