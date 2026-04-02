#!/bin/bash
#SBATCH --mem=50g
#SBATCH --cpus-per-task=24

INPUT=$1
WORK_DIR=$2
PRESET=$3
REF=$4

HitSV_TOOL="HitSV call "

mkdir ${WORK_DIR}
cd ${WORK_DIR}

BAM_ORI=${INPUT}

VCF=${WORK_DIR}/HitSV.vcf
LOG=${WORK_DIR}/HitSV.log

#call_single_chr
task_call_single_chr() {
    echo "Starting task task_call_single_chr-$1"
    ${HitSV_TOOL} -S $1 -E $1 -p ${PRESET} -l ${BAM_ORI} -r ${REF} -o ${WORK_DIR}/PART_$1_D.vcf >${WORK_DIR}/PART_$1_D_1.log 2> ${WORK_DIR}/PART_$1_D_2.log
    echo "Finished task task_call_single_chr-$1"
}

#combine
task_combine() {
    echo "Starting task combine"
    cat ${WORK_DIR}/PART_0_D.vcf > ${VCF}
    for i in {1..23}; do
        echo ${WORK_DIR}/PART_${i}_D.vcf
        cat ${WORK_DIR}/PART_${i}_D.vcf | grep -v "#" >> ${VCF}
    done
    echo "Finished task combine"
}

for i in {0..23}; do
    task_call_single_chr $i &
done
wait
task_combine
