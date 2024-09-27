We recommend to create the environment before running the script.
``` 
conda create --name env bioconda::cutadapt bioconda::seqkit bioconda::porechop bioconda::spades bioconda::minimap2 bioconda::samtools bioconda::centrifuge -y
```

Please install the following R pacakges before running the script
1. msa (https://github.com/UBod/msa)
2. Biostrings (https://github.com/Bioconductor/Biostrings)
3. gtools (https://github.com/r-gregmisc/gtools)
4. data.table (https://github.com/Rdatatable/data.table)

Also, the database file must be downloaded from the centrifuge website and keep in `db directory`.

```
wget https://genome-idx.s3.amazonaws.com/centrifuge/p_compressed%2Bh%2Bv.tar.gz
mkdir db
tar -xvzf p_compressed+h+v.tar.gz
mv p_compressed* db/
```


The script can be run with toy example dataset by the following command where the first and second argument is compress fastq file and output directory, respectively. 

```

./ensembled_denovo.sh sample/sample.fq.gz output

```

The output consists of three files as shown in the `example_output` directory

1. AA_best.fa
2. DNA_best.fa
3. Depth_final.txt

In this example, 4 out of 42 sets of hyperpameters were selected as the optimal hyperparameters as it can assemble the genome and translate into full AA length in this sample.
Any of them can be used as the representative of this sample for the downstream analysis.
