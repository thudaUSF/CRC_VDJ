#!/bin/bash
#SBATCH --job-name=SRAbamslice
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=8GB
#SBATCH --output=output.%j.SRAbamslice
#SBATCH --partition=rra
#SBATCH --qos=rra
#change this number to the total number of files in your manifest, ie the number 372 do not change the %100
#SBATCH --array=0-3%1

module purge
module add apps/sratoolkit/2.10.7
module add apps/samtools/1.3.1
module add apps/bwa/0.7.17
module add apps/sambamba/0.6.6

NCORE=8 #Total nodes divided by number of things going on at the same time (6 nodes / 3 tasks)
MEM=64GB #

#THE BELOW MUST BE CHANGED PER PROJECT; nothing else must be changed
Path="/work/pi_gblanck/Taha/"
PathToManifest="/work/pi_gblanck/Taha/nofqnobam.txt" #First, setup listunfinished.sh (After editting script to reference appropriate directories, edit sraruntable by filtering everything except WXS files and spotlength>140 (for xHLA min. 70bp), also just keep two columns, SRRID and subject ID) then run 'sh listunfinished.sh' to get nofqnobam.txt; should be in Unix format, dos2unix 
Token="/work/pi_gblanck/Taha/prj_20312_D25671.ngc" #make sure token is appropriate for the project and placed in folder (no need to set up with 'vdb-c	onfig -i')
srafolder="/home/t/thuda/sratoolkit/sra/"
OutputFolder="/work/pi_gblanck/Taha/phs001169fq/" #Where fastq files, full aligned bams go
FinalOutput="/work/pi_gblanck/Taha/phs001169bams/" #Where final sliced bam files go

#Set array of all SRRIDs to downloads
mapfile -t myArray < $PathToManifest
NumberOfBams=${#myArray[@]}
InputString=${myArray[$SLURM_ARRAY_TASK_ID]}
SRRID=$InputString #$(cut -d' ' -f1 <<< $InputString) if there's a comma afterwards

normalsize=500000 #at least 0.5MB
#Switch to SRA folder, because prefetch AND FASTER-Q only works from that folder (if not, fasterq will try to download and dump to fastq at the same time)
cd $srafolder

#Check if file is there, if "sizeable", don't do anything if the final file is already there
if [ -f "${FinalOutput}/${SRRID}.bam" ]; then
	filesize=$(stat -c%s "${FinalOutput}/${SRRID}.bam")
fi

if (( filesize > normalsize)); then
	echo "${SRRID}.bam already located with filesize:${filesize}, exitting"
	exit 0
fi
	
echo $SRRID
#prefetch fails the first time, so just run it first once and then re-run it
prefetch --ngc $Token $SRRID

#Point specifically to the full sra file
srafile="${srafolder}${SRRID}/${SRRID}_dbGaP-20312.sra" #change to end with project#
cd $OutputFolder
time fasterq-dump --ngc $Token $SRRID -O $OutputFolder -t $srafolder -e 8 # -e options is threads to work, 6 is default, you can remove time here, but it outputs time which is nice

#Check filesize if aligned bam exists; this stuff is a precaution against re-aligning if something substantial is already there
if [ -f "${OutputFolder}/${SRRID}.bam" ]; then 
	filesize=$(stat -c%s "${OutputFolder}/${SRRID}.bam")
fi

cd $OutputFolder 
#Only run if the bam is not there or not past minimum size (usually failed because prefetch didn't work so file doesn't exist)
if [ ! -f "${OutputFolder}/${SRRID}.bam" ] || (( filesize < normalsize)); then
	if [ -f "${OutputFolder}/${SRRID}.fastq" ]; then
	bwa mem -t $NCORE "/work/pi_gblanck/Taha/supportfiles/hg38.fna" "${SRRID}.fastq" | samtools view -b - | sambamba sort -m $MEM -t $NCORE -o "${SRRID}.bam" /dev/stdin
	fi
	if [ -f "${OutputFolder}/${SRRID}_1.fastq" ]; then
	bwa mem -t $NCORE "/work/pi_gblanck/Taha/supportfiles/hg38.fna" "${SRRID}_1.fastq" "${SRRID}_2.fastq" | samtools view -b - | sambamba sort -m $MEM -t $NCORE -o "${SRRID}.bam" /dev/stdin #Some SRRs are paired with _1,_2, some are just one - if I just let one or the other fail, it'll write over the initial thing, so I have to separate it
	fi
	
	sambamba index "${SRRID}.bam"

	samtools view -o "${SRRID}HLA.bam" -b "${SRRID}.bam" chr6:29844528-33100696  
	samtools view -o "${SRRID}TRA.bam" -b "${SRRID}.bam" chr14:21621904-22552132 
	samtools view -o "${SRRID}TRB.bam" -b "${SRRID}.bam" chr7:142299011-142813287 
	samtools view -o "${SRRID}TRG.bam" -b "${SRRID}.bam" chr7:38240024-38368055 
	samtools view -o "${SRRID}IGH.bam" -b "${SRRID}.bam" chr14:105586437-106879844 
	samtools view -o "${SRRID}IGK.bam" -b "${SRRID}.bam" chr2:88857361-90235368 
	samtools view -o "${SRRID}IGL.bam" -b "${SRRID}.bam" chr22:22026076-22922913 
	#WRITE IN UNNMAPPED

	samtools merge -f $FinalOutput$SRRID.bam "${SRRID}TRA".bam "${SRRID}TRB".bam "${SRRID}TRG".bam "${SRRID}IGH".bam "${SRRID}IGK".bam "${SRRID}IGL".bam "${SRRID}HLA".bam

	cd $FinalOutput
	sambamba index "${SRRID}".bam
fi

rm $srafolder$SRRID/* || true
rmdir $srafolder$SRRID || true
rm "${OutputFolder}/${SRRID}.fastq" || true
rm "${OutputFolder}/${SRRID}_1.fastq" || true
rm "${OutputFolder}/${SRRID}_2.fastq" || true
rm "${OutputFolder}/${SRRID}"* || true
cd $OutputFolder
find . -type f -name "${SRRID}"\* 
find . -type f -name "${SRRID}"\* -exec rm {} \

#Old Receptor Chromosome Locations
#chr14:21621904-22552132 TRA
#chr7:142299011-142813287 TRB
#chr14:22422546-22466577 TRD (within TRA, so no need to be specific)
#chr7:38240024-38368055 TRG
#chr14:105586437-106879844 IGH
#chr2:88857361-90235368 IGK
#chr22:22026076-22922913 IGL
#chr6:29941260-29945884 HLA-A 
#chr6:31269491-31357188 HLA-B
#chr6:31268749-31272130 HLA-C
#chr6:33075990-33089696 HLA-DPB1
#chr6:32659467-32668383 HLA-DQB1
#chr6:32578769-32589848 HLA-DRB1
