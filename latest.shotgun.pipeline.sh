# shotgun pipeline in bash on wonderwm

# create a rename.list from meta sheet from Anitha

####### set up output folder ########
# renamed folder to hold renamed reads #
mkdir renamed
mkdir deHost
mkdir contigs
mkdir contigs/spades
mkdir contigs/megahit
mkdir contigs/quast
mkdir contigs/prokka
mkdir contigs/emapper
mkdir contigs/prodigal
mkdir humann
mkdir MAGs
mkdir MAGs/spades
mkdir MAGs/megahit
mkdir taxonomy
mkdir taxonomy/readLevel
mkdir taxonomy/metaphlan
mkdir taxonomy/contigs
mkdir taxonomy/readLevel/bracken
mkdir targeted
mkdir targeted/card_bai
mkdir targeted/cazyme
mkdir targeted/VFDB
mkdir targeted/Blauticin_LCA

#######  rename fastqs  ####### 
for l in `cat rename.list`
do
	n=$(echo $l | cut -d ":" -f 1)
	rn=$(echo $l | cut -d ":" -f 2)
	f=$(ls rawReads/${n}_*R1*fastq*)
	echo "mv $f renamed/${rn}_R1.fastq.gz"
	r=$(ls rawReads/${n}_*R2*fastq*)
	echo "mv $r renamed/${rn}_R2.fastq.gz"
done > rename.sh
chmod 755 rename.sh
./rename.sh

#######   de-Host using kneaddata ########
conda activate base
for f in renamed/*R1.fastq.gz
do
    r=$(echo $f | sed 's/R1.fastq/R2.fastq/')
    n=$(basename $f | cut -d "_" -f 1)
echo "kneaddata --input ${f} --input ${r} --reference-db ~/Downloads/kneaddataDB/hg37dec_v0.1 --reference-db ~/Downloads/kneaddataDB/mouse_C57BL_6NJ --output ./deHost/${n} --bypass-trf --trimmomatic ~/Downloads/Trimmomatic-0.39/ -t 32 --trimmomatic-options=\"ILLUMINACLIP:~/Downloads/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:1:TRUE\" --trimmomatic-options=\"LEADING:3\" --trimmomatic-options=\"TRAILING:3\" --trimmomatic-options=\"SLIDINGWINDOW:4:15\" --trimmomatic-options=\"MINLEN:50\"  "
echo "rm deHost/${n}/*contam*fastq"
echo "rm deHost/${n}/*unmatched_*fastq"
echo "rm deHost/${n}/*trimmed*fastq"
echo "gzip deHost/${n}/*kneaddata_paired_*fastq"
done > kneaddata.trim.sh
chmod 755 kneaddata.trim.sh
./kneaddata.trim.sh # this has to be run on command line, no background

# gather kneaddate stats
grep "READ COUNT:" deHost/*/*log > trim.stats.txt

# after copy trim.stats.txt to local shotgun folder: run
# $ Rscript ~/OneDrive\ -\ The\ University\ of\ Chicago/DFIutility/produceShotgunQual.R 
# it will produce a trimStats.rds file and a quality plot for the entire run

#######   readLevel taxonomy profiling ########
# kraken2 and bracken
conda activate kraken2
for f in deHost/*/*paired_1.fastq.gz
do
n=$(basename $f | cut -d "_" -f 1-2 )
r=$(echo $f | sed 's/paired_1/paired_2/')
kraken2 --use-names --confidence 0.1 --report taxonomy/readLevel/${n}_report.txt --paired -threads 20 --db ~/Downloads/kraken2db --output taxonomy/readLevel/${n}.standard.txt $f $r
# since bracken is not installed on wonderwm, the report files need to copied to superman /mnt/Data/sdb/shotgun for bracken
bracken -d /mnt/Data/sda/Databases/brackenDB -i taxonomy/readLevel/${n}_report.txt -r 100 -t 10 -l S -o taxonomy/readLevel/bracken/${n}_report_species.txt
bracken -d /mnt/Data/sda/Databases/brackenDB -i taxonomy/readLevel/${n}_report.txt -r 100 -t 10 -l G -o taxonomy/readLevel/bracken/${n}_report_genus.txt
bracken -d /mnt/Data/sda/Databases/brackenDB -i taxonomy/readLevel/${n}_report.txt -r 100 -t 10 -l F -o taxonomy/readLevel/bracken/${n}_report_family.txt
bracken -d /mnt/Data/sda/Databases/brackenDB -i taxonomy/readLevel/${n}_report.txt -r 100 -t 10 -l O -o taxonomy/readLevel/bracken/${n}_report_order.txt
bracken -d /mnt/Data/sda/Databases/brackenDB -i taxonomy/readLevel/${n}_report.txt -r 100 -t 10 -l P -o taxonomy/readLevel/bracken/${n}_report_phylum.txt
bracken -d /mnt/Data/sda/Databases/brackenDB -i taxonomy/readLevel/${n}_report.txt -r 100 -t 10 -l C -o taxonomy/readLevel/bracken/${n}_report_class.txt
done
# parse bracken and kraken2 results
# $ Rscript ~/OneDrive\ -\ The\ University\ of\ Chicago/DFIutility/getKraken2std.R ./taxonomy/readLevel/

# if run only on report for bracken on Superman
for r in *_report.txt
do
n=$(echo $r | cut -d "_" -f 1-2)
bracken -d /mnt/Data/sda/Databases/brackenDB -i ${r} -r 100 -t 10 -l S -o ./bracken/${n}_report_species.txt
bracken -d /mnt/Data/sda/Databases/brackenDB -i ${r} -r 100 -t 10 -l G -o ./bracken/${n}_report_genus.txt
bracken -d /mnt/Data/sda/Databases/brackenDB -i ${r} -r 100 -t 10 -l F -o ./bracken/${n}_report_family.txt
bracken -d /mnt/Data/sda/Databases/brackenDB -i ${r} -r 100 -t 10 -l O -o ./bracken/${n}_report_order.txt
bracken -d /mnt/Data/sda/Databases/brackenDB -i ${r} -r 100 -t 10 -l P -o ./bracken/${n}_report_phylum.txt
bracken -d /mnt/Data/sda/Databases/brackenDB -i ${r} -r 100 -t 10 -l C -o ./bracken/${n}_report_class.txt
done

# parse bracken and kraken2 results
# $ Rscript ~/Documents/Eddi/DFIutility/getBracken.R ./taxonomy/readLevel/bracken/

#metaphlan4
conda activate mpa
for f in deHost/*/*paired_1.fastq.gz
do
n=$(basename $f | cut -d "_" -f 1-2 )
r=$(echo $f | sed 's/paired_1/paired_2/')
metaphlan ${f},${r} -t rel_ab_w_read_stats --bowtie2out taxonomy/metaphlan/${n}.bowtie2.bz2 --nproc 24 --input_type fastq --unclassified_estimation -o taxonomy/metaphlan/${n}_metagenome.txt --bowtie2db ~/Downloads/mpaDB
done

# parse metaphlan4 results
# $ Rscript ~/OneDrive\ -\ The\ University\ of\ Chicago/DFIutility/getMetaphlan.R ./taxonomy/metaphlan

#######  still in development:  humann3 functional annotation ####### 
conda activate humann
for f in deHost/*/*paired_1.fastq.gz
do
n=$(basename $f | cut -d "_" -f 1 )
r=$(echo $f | sed 's/paired_1/paired_2/')
humann --input $f --output humann/${n}_R1 --threads 32 --metaphlan-options "--bowtie2db /Users/PamerLabUser/Downloads/mpaDB"              
humann --input $r --output humann/${n}_R2 --threads 32 --metaphlan-options "--bowtie2db /Users/PamerLabUser/Downloads/mpaDB"   
done

#######   targeted panel #######   
for f in deHost/*/*paired_1.fastq.gz
do
n=$(basename $f | cut -d "_" -f 1-2 )
r=$(echo $f | sed 's/paired_1/paired_2/')
diamond blastx --threads 20 --db ~/Downloads/diamondDB/card_bai --sensitive -q $f -o targeted/card_bai/${n}.R1.card_bai.dmnd --log --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
diamond blastx --threads 20 --db ~/Downloads/diamondDB/card_bai --sensitive -q $r -o targeted/card_bai/${n}.R2.card_bai.dmnd --log --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
diamond blastx --threads 20 --db ~/Downloads/cazymeDB/CAZyDB.07312020 --sensitive -q $f -o targeted/cazyme/${n}.R1.cazyme.dmnd --log --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
diamond blastx --threads 20 --db ~/Downloads/cazymeDB/CAZyDB.07312020 --sensitive -q $r -o targeted/cazyme/${n}.R2.cazyme.dmnd --log --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
diamond blastx --threads 28 --db ~/Downloads/VFDB/VFDB_setB_pro --sensitive -q $f -o targeted/VFDB/${n}.R1.VFDB.dmnd --log --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
diamond blastx --threads 28 --db ~/Downloads/VFDB/VFDB_setB_pro  --sensitive -q $r -o targeted/VFDB/${n}.R2.VFDB.dmnd --log --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
done

# parse targeted results
# $ Rscript ~/OneDrive\ -\ The\ University\ of\ Chicago/DFIutility/getTargeted.R ./targeted

########### optional: blauticin + LCA ##############
for f in deHost/*/*paired_1.fastq.gz
do
n=$(basename $f | cut -d "_" -f 1-2 )
r=$(echo $f | sed 's/paired_1/paired_2/')
diamond blastx --threads 20 --db /Users/PamerLabUser/Downloads/exDB/Blauticin_LCA --sensitive -q $f -o targeted/Blauticin_LCA/${n}.R1.card_bai.dmnd --log --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
diamond blastx --threads 20 --db /Users/PamerLabUser/Downloads/exDB/Blauticin_LCA --sensitive -q $r -o targeted/Blauticin_LCA/${n}.R2.card_bai.dmnd --log --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
done

#######   MAG route ########
# megahit assembly
for f in deHost/*/*paired_1.fastq.gz
do
n=$(basename $f | cut -d "_" -f 1-2 )
r=$(echo $f | sed 's/paired_1/paired_2/')
megahit -t 30 -1 $f -2 $r -o contigs/megahit/${n}
mv contigs/megahit/${n}/final.contigs.fa contigs/megahit/${n}.megahit.contigs.fasta
prodigal -i contigs/megahit/${n}.megahit.contigs.fasta -c -m -g 11 -p meta -o contigs/prodigal/${n}.pr -q -d contigs/prodigal/${n}.fna -a contigs/prodigal/${n}.faa
rm -rf contigs/megahit/${n}
done

# emapper from prodigal
for c in contigs/prodigal/*.faa
do
n=$(echo $c | cut -d "/" -f 3 | sed 's/.faa//')
echo "python ~/Downloads/eggnog-mapper-2.0.1b/emapper.py -i $c --cpu 16 --output contigs/emapper/${n} -m diamond -d ~/Downloads/eggnog-mapper-2.0.1b/data/eggnog_proteins.dmnd"
# echo "python ~/Downloads/eggnog-mapper/emapper.py -i $c --cpu 16 --output contigs/emapper/${n} -m diamond -d ~/Downloads/eggnog_db/eggnog_proteins.dmnd"
done > emapper.prodigal.sh
chmod 775 emapper.prodigal.sh
./emapper.prodigal.sh # this has to be run without conda, python v2.7

# parse targeted results
# $ Rscript ~/OneDrive\ -\ The\ University\ of\ Chicago/DFIutility/getEmapper.R ./contigs/emapper

# after removing standard kraken2 output, emapper, and targeted diamond output, transfer the entire folder to /gpfs/data/dfi-cores/Bioinfo/shotgun
