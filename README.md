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

