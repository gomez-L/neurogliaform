# Code for replicating reported findings on: "Developmental Emergence of Neurogliaform cell Diversity"

![cover_5](https://github.com/gomez-L/neurogliaform/assets/8470916/8b42a0af-3e9c-4fcb-86e4-8207312a8901)
#### Due to GitHub file storage limitations, source data for replicating this work is only accessible at https://osf.io/se7gu, NCBI GEO GSE225639 or as supplementary material for the research paper:

Lucia Gomez, Christelle Cadilhac, Julien Prados, Nand Mule, Denis Jabaudon, Alexandre Dayer; Developmental emergence of cortical neurogliaform cell diversity. Development 2023; dev.201830. doi: https://doi.org/10.1242/dev.201830

ORCID ID Lucia Gomez: https://orcid.org/0000-0002-3914-8868

GABAergic interneurons are key regulators of cortical circuit function. Among the dozens of reported transcriptionally distinct subtypes of cortical interneurons, neurogliaform cells (NGCs) are unique: they are recruited by long-range excitatory inputs, are a source of slow cortical inhibition and are able to modulate the activity of large neuronal populations. Despite their functional relevance, the developmental emergence and diversity of NGCs remains unclear. Here, by combining single-cell transcriptomics, genetic fate mapping, electrophysiological and morphological characterization, we reveal that discrete molecular subtypes of NGCs, with distinctive anatomical and molecular profiles, populate the mouse neocortex. Furthermore, we show that NGC subtypes emerge gradually through development, as incipient discriminant molecular signatures are apparent in preoptic area (POA)-born NGC precursors. By identifying NGC developmentally conserved transcriptional programs, we report that the transcription factor Tox2 constitutes an identity hallmark across NGC subtypes. Using CRISPR-Cas9-mediated genetic loss of function, we show that Tox2 is essential for NGC development: POA-born cells lacking Tox2 fail to differentiate into NGCs. Together, these results reveal that NGCs are born from a spatially restricted pool of Tox2+ POA precursors, after which intra-type diverging molecular programs are gradually acquired post-mitotically and result in
functionally and molecularly discrete NGC cortical subtypes.

Programming Languages: R, Python, Bash

Raw data for this project is available in NCBI Gene Expression Omnibus under accession number GSE225639 

Contents here provided are also available at Open Science Framework https://osf.io/se7gu 

## Guidelines for replication
1. Download all contents preferably from OSF repository, here only files under 25MB (https://osf.io/se7gu) </br>
2. Use R to execute .Rmd or .R files (scripts): </br>
    * Order of scripts: </br>
          * CCA_KNN_SVM_T18_P15.Rmd </br>
          * CCA_KNN_SVM_T18_P30.Rmd </br>
          * Postnatal_Type_Architecture.Rmd </br>
          * Postnatal_Subtype_Architecture.Rmd </br>
          * lib/go_hgnc.R </br>
          * Embryonic_pseudotime.Rmd </br>
          * NGC_timeconserved_genes.Rmd </br>

Contact the main author and developer at: lucy._.y@hotmail.com
      
