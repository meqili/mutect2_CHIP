#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -V
#$ -m bea
#$ -pe threaded 2
#$ -t 1-32:1
#$ -tc 12
#$ -N M2_snpEff
#$ -o /nfs/projects/dbGap/mutect2/mutect_CHIP/logs
#$ -e /nfs/projects/dbGap/mutect2/mutect_CHIP/logs

gatk_path="/home/ql2387/software/gatk-4.2.5.0/gatk"
ref_fasta="/nfs/seqscratch09/AZ-IPF/reference/hs37d5.fa"
BAM_LIST="/nfs/projects/dbGap/mutect2/mutect_CHIP/all_IPF_recalibrated.bam.list"

seq $SGE_TASK_ID $(($SGE_TASK_ID + $SGE_TASK_STEPSIZE - 1)) | while read -r line; do input_bam=$(sed -n "$line p" $BAM_LIST); base_file_name=$(basename ${input_bam} '.realn.recal.bam'); 
mkdir ${base_file_name};

${gatk_path} Mutect2 \
   -R $ref_fasta \
   -pon /nfs/projects/dbGap/mutect2/resource/Mutect2-WGS-panel-b37.vcf \
   -I  $input_bam \
   --germline-resource /nfs/projects/dbGap/mutect2/resource/af-only-gnomad.raw.sites.vcf \
   --af-of-alleles-not-in-resource 0.00000005 \
   -L 19:17935489-17958980 -L 4:55523985-55606981 -L 17:29421845-29709234 -L 10:89622770-89731787 -L 20:57414673-57486347 -L 1:36931544-36948979 -L 21:36159998-37377065 -L 10:8095467-8117261 -L X:133507183-133562920 -L 12:49961901-50038549 -L X:39908968-40036682 -L 17:7564997-7590956 -L 4:106066932-106201073 -L 1:150898639-150937313 -L 3:128198170-128212128 -L 7:148504375-148581513 -L 22:41487690-41576181 -L 12:11802688-12048436 -L X:48644862-48652816 -L 1:1716629-1822595 -L 15:90626177-90645836 -L 1:65298812-65432287 -L 5:170814020-170838241 -L 1:115246990-115259615 -L 12:112856055-112947817 -L 17:74730097-74733556 -L 12:49412658-49453657 -L 17:58677444-58741949 -L 11:85955486-85989955 -L 13:28577311-28674829 -L 7:140419027-140624664 -L 19:33790740-33793570 -L 16:3774955-3930827 -L 11:64531978-64546358 -L 9:4984933-5128283 -L X:15808495-15841483 -L 11:119076652-119178959 -L 5:149432754-149493035 -L 1:43803378-43818543 -L 10:112327349-112364494 -L 3:136054977-136471320 -L 22:30727877-30753036 -L 2:209100851-209130898 -L 11:118307105-118397639 -L 7:50343620-50472899 -L 8:117858074-117887205 -L X:53400970-53449777 -L 17:1553823-1588276 -L 16:67596210-67673186 -L 17:63006733-63053057 -L X:123093962-123556614 -L 21:44512966-44527797 -L 19:56165412-56186181 -L 2:198254408-198299915 -L 17:37921098-38020541 -L 3:47057819-47205557 -L 13:33160464-33352257 -L 20:30946055-31027222 -L 17:30263937-30328164 -L 18:42260038-42648575 -L 7:139025005-139108298 -L 2:213864329-214017251 -L 6:79645484-79788053 -L 2:25956522-26101485 -L 12:22777909-22843699 -L 6:107473661-107780868 -L X:44732657-44971947 -L X:129114983-129192158 -L 2:25455745-25565559 -L X:154299595-154351449 -L 11:32409221-32457276 -L 12:25357623-25403970 -L 3:105374205-105588496 -L 7:101458859-101927349 \
   -O ${base_file_name}/${base_file_name}"-unfiltered.vcf.gz"

${gatk_path} FilterMutectCalls \
   -R $ref_fasta \
   -V ${base_file_name}/${base_file_name}"-unfiltered.vcf.gz" \
   -O ${base_file_name}/${base_file_name}"-filtered.vcf.gz"

SNPEFF_FILE=${base_file_name}-filtered.ann.vcf.gz;  
if [ ! -f $SNPEFF_FILE ]; then 
   java -Xmx8g -jar /nfs/goldstein/software/snpEff/4.3/snpEff.jar ann -v GRCh37.75 \
     ${base_file_name}/${base_file_name}"-filtered.vcf.gz"  >  $SNPEFF_FILE
fi

done
