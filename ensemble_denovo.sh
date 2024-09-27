#conda create --name env bioconda::cutadapt bioconda::seqkit bioconda::porechop bioconda::spades bioconda::minimap2 bioconda::samtools bioconda::centrifuge -y  
inputFQ=$1
output_name=$2

primer_=GTTTCCCACTGGAGGATA
primer_rev_com=TATCCTCCAGTGGGAAAC


echo "###### STEP 1: pull SISPA PCR product ######"
mkdir $output_name
zcat $inputFQ | seqkit grep -R 1:50 -m 0 -s -i -P -p $primer_ | seqkit fx2tab | cut -d' ' -f 1 | sort | uniq > ${output_name}/Freads.txt
seqkit grep -f ${output_name}/Freads.txt $inputFQ | seqkit grep -R -50:-1 -m 0 -s -i -P -p $primer_rev_com -o ${output_name}/sispa_read.fq


echo "###### STEP 2: trim primer + nine random mers  ######"
cutadapt -g GTTTCCCACTGGAGGATA...TATCCTCCAGTGGGAAAC -o ${output_name}/trimmed_read.fq ${output_name}/sispa_read.fq
cutadapt -u 9 -o ${output_name}/trimmed_read_1.fq ${output_name}/trimmed_read.fq
cutadapt -u -9 -o ${output_name}/trimmed_read_2.fq ${output_name}/trimmed_read_1.fq

echo "###### STEP 3: discard adapter  ######"
porechop -i ${output_name}/trimmed_read_2.fq -o ${output_name}/clean_read.fq --discard_middle 
    
    
echo "###### STEP 4: select dengue 2 reads from centrifuge  ######"	
centrifuge -x db/p_compressed+h+v -U ${output_name}/clean_read.fq --report-file ${output_name}/report.txt -S ${output_name}/results.txt 

#modify the header name
path=$(echo ${output_name}/report.txt)
id=$(cat $path | grep "Dengue" | sort -t$'\t' -k5 -n | tail -n 1 | cut -f2)
den_type=$(cat $path | grep "Dengue" | sort -t$'\t' -k5 -n | tail -n 1 | cut -f1 | tr ' ' '_')
fl_=$(echo $path | sed 's\report\results\' | sed 's\Dengue virus 2\Dengue_virus_2\' )
grep $id $fl_ | cut -f1 > ${output_name}/${den_type}.txt
seqkit grep -f ${output_name}/${den_type}.txt ${output_name}/clean_read.fq | gzip -9  >${output_name}/viral_read.fq.gz

echo "###### STEP 5: genome assembly +  grid search ######"
for quality in 20 18 16 14 12 10;
do   
    for kmers in 127 111 99 77 55 33 11;
    do         
	    #5.1 assemble
	    output_SPAdes=$(echo ${output_name}/assembly/${quality}_${kmers})
	    
	    gunzip -c ${output_name}/viral_read.fq.gz | NanoFilt -q $quality | gzip -9 > ${output_name}/viral_read_Q_${quality}.fq.gz
        spades.py -s ${output_name}/viral_read_Q_${quality}.fq.gz -k $kmers -o  ${output_SPAdes}/
        
        #5.2 filter the contig with length over 1k bps
        seqkit seq -m 10000 ${output_SPAdes}/scaffolds.fasta > ${output_SPAdes}/QC_scaffold.fasta
        sed -e "s/length/K${kmers}/g" -e "s/cov/${quality}/g" ${output_SPAdes}/QC_scaffold.fasta > ${output_SPAdes}/QC_scaffold_final.fasta
        
        #5.3  AA translation from six possible ORF
        for ORF in 1 2 3 -1 -2 -3;
        do
            seqkit translate ${output_SPAdes}/QC_scaffold_final.fasta -f $ORF -s -M >> ${output_SPAdes}/AA.fasta 
        done 
        
        seqkit sort -r -l ${output_SPAdes}/AA.fasta | seqkit range -r 1:1 > ${output_SPAdes}/AA_final.fasta 
        #5.4 ensure that AA start with M 
        Num=$(($(cat ${output_SPAdes}/AA_final.fasta  | tail -n +2 | grep M -bo | head -n +1 | cut -d':' -f1)+1))
        seqkit subseq -r $Num:-1 ${output_SPAdes}/AA_final.fasta > ${output_SPAdes}/QC_AA_final.fasta     
        
        #5.5 estimate depth of assembled genome
        minimap2 -a ${output_SPAdes}/QC_scaffold_final.fasta ${output_name}/viral_read_Q_${quality}.fq.gz > ${output_SPAdes}/map_read_virus.sam
        samtools sort ${output_SPAdes}/map_read_virus.sam -o ${output_SPAdes}/map_read_virus.bam
        samtools index ${output_SPAdes}/map_read_virus.bam
        samtools depth ${output_SPAdes}/map_read_virus.bam > ${output_SPAdes}/depth_virus.txt        
    done
done   


echo "###### STEP 6: ensemble the assembled genome ######"      
Rscript ensemble.R ${output_name}















