#!/bin/bash

MODE=""
HitSV_TOOL=""
WORK_DIR=""
REF=""
FA_IDX=""
PRESET=""
PRESET_KSW=""
INPUT_SRS=""
INPUT_LRS=""
MAX_PARALLEL=20
MAX_CHR=2000

usage() {
    cat <<EOF
Usage: $0 <MODE> [options]

Modes:
  LRS     Long-read sequencing only
  SRS     Short-read sequencing only
  Hybrid  Combined SRS + LRS analysis

Options:
  -t TOOL       Path to HitSV tool
  -w WORK_DIR   Working directory
  -r REF        Reference genome
  -i FA_IDX     FA index (required for SRS/Hybrid)
  -p PRESET     Preset (required for LRS/Hybrid)
  -P PRESET_KSW KSW preset (required for all modes)
  -n INPUT_SRS  SRS input BAM (required for SRS/Hybrid)
  -l INPUT_LRS  LRS input BAM (required for LRS/Hybrid)
  -m MAX_PARALLEL  Max parallel jobs (default: 20)
  -c MAX_CHR       Max chromosome index (default: 2000)

Examples:
  LRS mode:
    $0 LRS -t /path/to/HitSV -p preset -P ksw_preset -l input.bam -w work_dir -r ref.fa

  SRS mode:
    $0 SRS -t /path/to/HitSV -P ksw_preset -n input.bam -w work_dir -r ref.fa -i ref.fa.fai

  Hybrid mode:
    $0 Hybrid -t /path/to/HitSV -p preset -P ksw_preset -n srs.bam -l lrs.bam -w work_dir -r ref.fa -i ref.fa.fai
EOF
    exit 1
}

if [ $# -lt 1 ]; then
    usage
fi

MODE=$1
shift

case "${MODE}" in
    LRS|SRS|Hybrid) ;;
    *)
        echo "Error: Invalid mode '${MODE}'. Must be LRS, SRS, or Hybrid."
        usage
        ;;
esac

while getopts "t:w:r:i:p:P:n:l:m:c:" opt; do
    case ${opt} in
        t) HitSV_TOOL=${OPTARG} ;;
        w) WORK_DIR=${OPTARG} ;;
        r) REF=${OPTARG} ;;
        i) FA_IDX=${OPTARG} ;;
        p) PRESET=${OPTARG} ;;
        P) PRESET_KSW=${OPTARG} ;;
        n) INPUT_SRS=${OPTARG} ;;
        l) INPUT_LRS=${OPTARG} ;;
        m) MAX_PARALLEL=${OPTARG} ;;
        c) MAX_CHR=${OPTARG} ;;
        *) usage ;;
    esac
done

if [ -z "${HitSV_TOOL}" ] || [ -z "${WORK_DIR}" ] || [ -z "${REF}" ] || [ -z "${PRESET_KSW}" ]; then
    echo "Error: Missing required parameters (-t, -w, -r, -P are required for all modes)."
    usage
fi

case "${MODE}" in
    LRS)
        if [ -z "${INPUT_LRS}" ] || [ -z "${PRESET}" ]; then
            echo "Error: LRS mode requires -l (INPUT_LRS) and -p (PRESET)."
            usage
        fi
        ;;
    SRS)
        if [ -z "${INPUT_SRS}" ] || [ -z "${FA_IDX}" ]; then
            echo "Error: SRS mode requires -n (INPUT_SRS) and -i (FA_IDX)."
            usage
        fi
        ;;
    Hybrid)
        if [ -z "${INPUT_SRS}" ] || [ -z "${INPUT_LRS}" ] || [ -z "${FA_IDX}" ] || [ -z "${PRESET}" ]; then
            echo "Error: Hybrid mode requires -n (INPUT_SRS), -l (INPUT_LRS), -i (FA_IDX), and -p (PRESET)."
            usage
        fi
        ;;
esac

chmod 777 ${HitSV_TOOL}

mkdir -p ${WORK_DIR}
cd ${WORK_DIR}

VCF=${WORK_DIR}/HitSV.vcf

if [ "${MODE}" = "SRS" ] || [ "${MODE}" = "Hybrid" ]; then
    STAT=${WORK_DIR}/stat.json
    ${HitSV_TOOL} srs_bam_stat ${REF} ${INPUT_SRS} > ${STAT}

    TL_ORI=${WORK_DIR}/TL.bam
    TL=${WORK_DIR}/TL.sort.bam

    ${HitSV_TOOL} srs_trans_reads ${REF} ${INPUT_SRS} ${TL_ORI}
    samtools sort -@ 2 --output-fmt=BAM -m5G -o ${TL} ${TL_ORI}
    samtools index ${TL}

    echo ${STAT}
fi

echo "Begin ALL [${MODE}]" > ${WORK_DIR}/D_ALL.log

task_call_single_chr() {
    local i="$1"
    echo "Starting task $i" >> ${WORK_DIR}/D_ALL.log

    case "${MODE}" in
        LRS)
            ${HitSV_TOOL} call -S $i -E $i -p ${PRESET} -P ${PRESET_KSW} -l ${INPUT_LRS} -r ${REF} \
                -o ${WORK_DIR}/PART_${i}_D.vcf 2> /dev/null
            ;;
        SRS)
            ${HitSV_TOOL} call -S $i -E $i -n ${INPUT_SRS} -r ${REF} -I ${FA_IDX} -T ${STAT} -L ${TL} -P ${PRESET_KSW} \
                -o ${WORK_DIR}/PART_${i}_D.vcf 2> /dev/null
            ;;
        Hybrid)
            ${HitSV_TOOL} call -S $i -E $i -n ${INPUT_SRS} -l ${INPUT_LRS} -r ${REF} -I ${FA_IDX} -T ${STAT} -L ${TL} -p ${PRESET} -P ${PRESET_KSW} \
                -o ${WORK_DIR}/PART_${i}_D.vcf 2> /dev/null
            ;;
    esac

    echo "Finished task $i" >> ${WORK_DIR}/D_ALL.log
}
export -f task_call_single_chr
export HitSV_TOOL WORK_DIR INPUT_SRS INPUT_LRS FA_IDX REF STAT TL PRESET PRESET_KSW MODE

seq 0 ${MAX_CHR} | xargs -P ${MAX_PARALLEL} -I {} bash -c 'task_call_single_chr {}'

task_combine() {
    cat ${WORK_DIR}/PART_0_D.vcf > ${VCF}
    for i in $(seq 1 ${MAX_CHR}); do
        grep -v "#" ${WORK_DIR}/PART_${i}_D.vcf >> ${VCF}
    done
}
task_combine
echo "Finish ALL [${MODE}]" >> ${WORK_DIR}/D_ALL.log
