# Variant STARR-seq paper
This repository contains the code used to process all the data and reproduce the results and figures from the manuscript "Enhancers display constrained sequence flexibility and context-specific modulation of motif function".

For more information, see the manuscript:  
[*<ins>Enhancers display constrained sequence flexibility and context-specific modulation of motif function</ins>*](https://www.biorxiv.org/content/10.1101/2022.08.31.506061v1)  
Franziska Reiter\*, Bernardo P. de Almeida\*, Alexander Stark. bioRxiv  
\* Equal contribution.

The raw sequencing data are available from GEO under accession number [GSE211659](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211659).  
The data necessary to reproduce all results and figures are available on zenodo at https://doi.org/10.5281/zenodo.7010528. Please download the respective datasets and follow the scripts below.    

## Random variant STARRseq
Random 8nt variants tested at seven enhancer positions in Drosophila S2 cells  
		- [Pipeline for mapping sequencing reads](Random_variant_STARRseq/Read_mapping_pipeline.sh)  
		- [R markdown to reproduce results](Random_variant_STARRseq/Random_variant_STARRseq_analysis.Rmd)  
		- [Results in html](https://rawcdn.githack.com/bernardo-de-almeida/Variant_STARRseq/77a09172310b16b021d5e0309619a5e934a26d21/Random_variant_STARRseq/Random_variant_STARRseq_analysis.html)  

<p align="center">
	<img src="img/Random_Variants_STARRseq_screen.png" width="700" style="margin-bottom:0;margin-top:0;"/>
</p>

## Drosophila motif pasting STARRseq
Systematic pasting of Drosophila TF motifs in hundreds of enhancer positions  
		- Sequencing reads were processed as in [here](https://github.com/bernardo-de-almeida/DeepSTARR/tree/main/Oligo_UMISTARRseq)  
		- [R markdown to reproduce results](Drosophila_motif_pasting_STARRseq/Drosophila_motif_pasting_STARRseq_analysis.Rmd)  
		- [Results in html](https://rawcdn.githack.com/bernardo-de-almeida/Variant_STARRseq/77a09172310b16b021d5e0309619a5e934a26d21/Drosophila_motif_pasting_STARRseq/Drosophila_motif_pasting_STARRseq_analysis.html)  

<p align="center">
	<img src="img/Drosophila_motif_pasting_experiment.png" width="700" style="margin-bottom:0;margin-top:0;"/>
</p>

## Human motif pasting STARRseq
Systematic pasting of human TF motifs in hundreds of enhancer positions  
		- Sequencing reads were processed as in [here](https://github.com/bernardo-de-almeida/DeepSTARR/tree/main/Oligo_UMISTARRseq)  
		- [R markdown to reproduce results](Human_motif_pasting_STARRseq/Human_motif_pasting_STARRseq_analysis.Rmd)  
		- [Results in html](https://rawcdn.githack.com/bernardo-de-almeida/Variant_STARRseq/77a09172310b16b021d5e0309619a5e934a26d21/Human_motif_pasting_STARRseq/Human_motif_pasting_STARRseq_analysis.html)  

<p align="center">
	<img src="img/Human_motif_pasting_result.png" width="400" style="margin-bottom:0;margin-top:0;"/>
</p>

## Questions
If you have any questions/requests/comments please contact me at [bernardo.almeida94@gmail.com](mailto:bernardo.almeida94@gmail.com).
