.rc_iupac <- function(s) {
  map <- c(A="T", C="G", G="C", T="A", M="K", R="Y", Y="R", W="W", S="S", K="M",
           V="B", H="D", D="H", B="V", N="N",
           a="t", c="g", g="c", t="a", m="k", r="y", y="r", w="w", s="s", k="m",
           v="b", h="d", d="h", b="v", n="n")
  ch <- strsplit(s, "", fixed = TRUE)[[1]]
  rc <- vapply(rev(ch), function(x) if (!is.na(map[x])) map[x] else "N", character(1))
  paste0(rc, collapse = "")
}

PRIMER_F_RC <- .rc_iupac(PRIMERSEQ_F)
PRIMER_R_RC <- .rc_iupac(PRIMERSEQ_R)

scripttext <- paste0(
  "
#!/bin/bash
set -e

PRIMERSEQ_F=", PRIMERSEQ_F, "
PRIMERSEQ_R=", PRIMERSEQ_R, "
PRIMERNAME=", PRIMERNAME, "

CUTADAPT=", shQuote(CUTADAPT), "
FASTQ_LOCATION=", shQuote(FASTQ_LOCATION), "

# reverse-complements for anchored second pass
PRIMER_F_RC=", PRIMER_F_RC, "
PRIMER_R_RC=", PRIMER_R_RC, "

cd \"${FASTQ_LOCATION}\"
FILES=$(ls ./${PRIMERNAME}-*R1*.fastq.gz)
mkdir -p '../for_dada2/'
mkdir -p '../logs/'

for i in $FILES
do
  FILE_NAME=$(echo ${i} | cut -d _ -f 1,2,3)
  echo \"[cutadapt pass 1] ${FILE_NAME}\"
  R1=${i}
  R2=$(echo ${i} | sed 's/R1/R2/g')

  # temp outputs for pass 1 inside for_dada2
  R1_BASE=$(basename \"${R1}\" .fastq.gz)
  TMP_R1=\"../for_dada2/${R1_BASE}_p1.fastq.gz\"

  if [ -f \"${R2}\" ]; then
    # -------------------- PAIRED-END --------------------
    R2_BASE=$(basename \"${R2}\" .fastq.gz)
    TMP_R2=\"../for_dada2/${R2_BASE}_p1.fastq.gz\"

    # ---- PASS 1: original, unanchored ----
    ${CUTADAPT} -g ${PRIMERSEQ_F} \\
                -G ${PRIMERSEQ_R} \\
                --rc \\
                -o \"${TMP_R1}\" \\
                -p \"${TMP_R2}\" \\
                --discard-untrimmed \\
                -j 0 \\
                \"${R1}\" \"${R2}\" \\
                >> '../logs/cutadapt_trim_report.txt' 2>&1 || true

    echo \"[cutadapt pass 2 (anchored)] ${FILE_NAME}\"
    # ---- PASS 2: anchored cleanup only; no discard ----
    R1_OUT_BASE=$(basename \"${R1}\")
    R2_OUT_BASE=$(basename \"${R2}\")
    ${CUTADAPT} -g \"^${PRIMERSEQ_F}\" \\
                -G \"^${PRIMERSEQ_R}\" \\
                -a \"${PRIMER_R_RC}\\$\" \\
                -A \"${PRIMER_F_RC}\\$\" \\
                -o \"../for_dada2/${R1_OUT_BASE}\" \\
                -p \"../for_dada2/${R2_OUT_BASE}\" \\
                -j 0 \\
                \"${TMP_R1}\" \"${TMP_R2}\" \\
                >> '../logs/cutadapt_trim_report.txt' 2>&1

    rm -f \"${TMP_R1}\" \"${TMP_R2}\"

  else
    # -------------------- SINGLE-END --------------------

    # ---- PASS 1: original, unanchored ----
    ${CUTADAPT} -g ${PRIMERSEQ_F} \\
                -a \"${PRIMER_R_RC}\\$\" \\
                --rc \\
                -o \"${TMP_R1}\" \\
                --discard-untrimmed \\
                -j 0 \\
                \"${R1}\" \\
                >> '../logs/cutadapt_trim_report.txt' 2>&1 || true

    echo \"[cutadapt pass 2 (anchored)] ${FILE_NAME} (single-end)\"
    # ---- PASS 2: anchored cleanup only; no discard ----
    R1_OUT_BASE=$(basename \"${R1}\")
    ${CUTADAPT} -g \"^${PRIMERSEQ_F}\" \\
                -a \"${PRIMER_R_RC}\\$\" \\
                -o \"../for_dada2/${R1_OUT_BASE}\" \\
                -j 0 \\
                \"${TMP_R1}\" \\
                >> '../logs/cutadapt_trim_report.txt' 2>&1

    rm -f \"${TMP_R1}\"
  fi
done
"
)

cat(scripttext, file = paste0(CODE_LOCATION, "/trim_primers.sh"))
