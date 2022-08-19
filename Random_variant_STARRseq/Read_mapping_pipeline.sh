#######################
## prepare files
######################

folder=Random_variant_STARRseq
cd $folder

# folder to write results
dataFolder=$folder/data
mkdir -p $dataFolder


# experiment file
experimentFile=$dataFolder/experiment.txt
head $experimentFile

# to submit to the cluster
bsub=bsub_gridengine
mkdir -p log_mapping

#######################
## Extract reads using mapping script + mapp reads
#######################

ALIGNER=src/bowtie_pe.sh

grep -v -E "^#" $experimentFile | \
    while read line; do
    INFILE=$( echo $line | awk '{print $2}' )
    outdir=$dataFolder
    BARCODES14=$( echo $line | awk '{print $9}' )
    BARCODES12=$( echo $line | awk '{print $11}' )
    OUTFILE=$( echo $line | awk '{print $10}')
    BARCODE14_LEN=$( echo $BARCODES14 | awk '{print length($1)}' )
    BARCODE12_LEN=$( echo $BARCODES12 | awk '{print length($1)}' )
    GENOME=$( echo $line | awk '{print $5}' )

    BARCODE14_LEN=10 # length of UMI

    # extract reads
    $bsub -o log_mapping -C 1 -n "${OUTFILE}_mapping" "$ALIGNER -i $INFILE -o ${outdir}/${OUTFILE}.bb -B $BARCODES12 -L $BARCODE12_LEN -l $BARCODE14_LEN --umi -f A -g $GENOME -F 1" > log_mapping/msg.${OUTFILE}_reads.tmp
    
done

# look at fasta files
zcat data/S21717_REEFmix2_Rep1_1.fasta.gz | grep ">" |wc -l
zcat data/S21717_REEFmix2_Rep1_2.fasta.gz | grep ">" |wc -l
# 9719001


######################################
## Adapted mapping - fw read 36bp 3mismatch / rv read 150bp 0mismatch
######################################

# 1 - cut fw reads to 36bp
# 2 - first change headers of fasta sequences to match, still keep UMI
grep -v -E "^#" $experimentFile | \
    while read line; do
    OUTFILE=$( echo $line | awk '{print $10}')

    bsub -o log_cutadapt "module load cutadapt/1.18-foss-2018b-python-3.6.6; cutadapt -f fasta -l 36 -o data/${OUTFILE}_1_36bp.fasta.gz data/${OUTFILE}_1.fasta.gz; \
    zcat data/${OUTFILE}_1_36bp.fasta.gz | awk -F '[/^>|_]' 'NF>1{print \">seq\" ++i \"_\" \$3} {print \$1}' | awk NF | gzip > data/${OUTFILE}_1_36bp_header.fasta.gz; \
    zcat data/${OUTFILE}_2.fasta.gz | awk -F '[/^>|_]' 'NF>1{print \">seq\" ++i \"_\" \$3} {print \$1}' | awk NF | gzip > data/${OUTFILE}_2_header.fasta.gz"
done

### map forward read allowing 3MM and multimapping (no -m --strara), just to see if reads map correctly to first sequence position
grep -v -E "^#" $experimentFile | \
    while read line; do
    OUTFILE=$( echo $line | awk '{print $10}')
    GENOME=$( echo $line | awk '{print $5}' )

    BOWTIE_MM=3
    bsub -o log_bowtie -C 1 -n "${OUTFILE}_mapping" "bowtie -p 1 --un data/${OUTFILE}_unmapped_${BOWTIE_MM}MM_fw.bowtie -f --best -v $BOWTIE_MM /groups/stark/indices/bowtie/${GENOME}/${GENOME} data/${OUTFILE}_1_36bp_header.fasta.gz > data/${OUTFILE}_${BOWTIE_MM}MM_fw.bowtie 2> data/${OUTFILE}_${BOWTIE_MM}MM_fw.log"

done

# compress bowtie files
BOWTIE_MM=3
grep -v -E "^#" $experimentFile | \
    while read line; do
    OUTFILE=$( echo $line | awk '{print $10}')
    bsub "gzip data/${OUTFILE}_${BOWTIE_MM}MM_fw.bowtie"
done

# check mapping rate
grep -v -E "^#" $experimentFile | \
    while read line; do
    OUTFILE=$( echo $line | awk '{print $10}')
    
    echo 
    echo $OUTFILE
    cat data/${OUTFILE}_3MM_fw.log

done



### select mates of forward reads that mapped to correct position (position 0 or 214)
grep -v -E "^#" $experimentFile | \
    while read line; do
    OUTFILE=$( echo $line | awk '{print $10}')
    bsub -o log_select_mates "zcat data/${OUTFILE}_3MM_fw.bowtie.gz | awk '{if(\$4==0 || \$4==213)print \$0}' | cut -f1 > data/${OUTFILE}_3MM_fw_reads; \
    zcat data/${OUTFILE}_2_header.fasta.gz | seqtk subseq - data/${OUTFILE}_3MM_fw_reads | gzip > data/${OUTFILE}_2_header_fw_mapped.fasta.gz; \
    rm data/${OUTFILE}_3MM_fw_reads"
    echo $OUTFILE
done

### map reverse read 150bp with no mismatches
grep -v -E "^#" $experimentFile | \
    while read line; do
    OUTFILE=$( echo $line | awk '{print $10}')
    GENOME=$( echo $line | awk '{print $5}' )
    BOWTIE_MM=0

    bsub -o log_bowtie -C 1 -n "${OUTFILE}_mapping_rv" "bowtie -p 1 --un data/${OUTFILE}_unmapped_${BOWTIE_MM}MM_rv.bowtie -f -m 1 --best --strata -v $BOWTIE_MM /groups/stark/indices/bowtie/${GENOME}/${GENOME} data/${OUTFILE}_2_header_fw_mapped.fasta.gz > data/${OUTFILE}_${BOWTIE_MM}MM_rv.bowtie 2> data/${OUTFILE}_${BOWTIE_MM}MM_rv.log"
done

# check mapping rate
grep -v -E "^#" $experimentFile | \
    while read line; do
    OUTFILE=$( echo $line | awk '{print $10}')
    
    echo 
    echo $OUTFILE
    cat data/${OUTFILE}_0MM_rv.log

done


######################################
## UMI collapsing
######################################

BOWTIE_MM=0

grep -v -E "^#" $experimentFile | \
    while read line; do
    OUTFILE=$( echo $line | awk '{print $10}')
    # remove sequence name from first column
    # collapse by UMI and position
    cat data/${OUTFILE}_${BOWTIE_MM}MM_rv.bowtie | awk '{sub(/\_/,"\t",$1)};1' OFS="\t" | cut -f1 --complement | awk -vOFS='\t' '{print $3,$4,$4+length($5),$2,$1}' | \
    sort -k1,1 -k2,2n -k3,3nr -k4,4 | \
    uniq -c | awk -vOFS='\t' '{print $2,$3,$4,$6,$1,$5}' | \
    sort -k1,1 -k2,2n -k3,3nr -k6,6 -k5,5nr > data/tmp_${OUTFILE}.reads.bowtie

    # run
    Rexec="singularity run --app Rscript singularity.R.with_pckg.simg"
    DIR=src
    bsub -o log_UMI_R -C 5 "$Rexec ${DIR}/STARRseq_UMI_collapsing.R -i data/tmp_${OUTFILE}.reads.bowtie -m 1 -c 5 -o data/${OUTFILE}_final_mapped.UMI.bed; \
    rm data/tmp_${OUTFILE}.reads.bowtie"

    echo $OUTFILE
done


### make all reads bed file
grep -v -E "^#" $experimentFile | \
    while read line; do
    OUTFILE=$( echo $line | awk '{print $10}')
    # remove sequence name from first column
    cat data/${OUTFILE}_${BOWTIE_MM}MM_rv.bowtie | awk '{sub(/\_/,"\t",$1)};1' OFS="\t" | cut -f1 --complement | awk -vOFS='\t' '{print $3,$4,$4+length($5),$1,0,$2}' > data/${OUTFILE}_final_mapped.all.bed

    echo $OUTFILE
done

