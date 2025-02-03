
scripttext <- paste0(
"
#!/bin/bash

PRIMERSEQ_F=", PRIMERSEQ_F,"
PRIMERSEQ_R=", PRIMERSEQ_R,"
PRIMERNAME=", PRIMERNAME,"

CUTADAPT=", shQuote(CUTADAPT),"
FASTQ_LOCATION=",shQuote(FASTQ_LOCATION),"


cd \"${FASTQ_LOCATION}\"
FILES=$(ls ./${PRIMERNAME}'-'*R1*.fastq.gz)
mkdir '../for_dada2/'
mkdir '../logs/'

for i in $FILES
  do
FILE_NAME=$(echo ${i} | cut -d _ -f 1,2,3) #grab everything before R1 (aka sample name)
echo $FILE_NAME
R1=${i}
R2=$(echo ${i} | sed 's/R1/R2/g')
${CUTADAPT} -g ${PRIMERSEQ_F} \\\
-G ${PRIMERSEQ_R} \\\
-o ../for_dada2/${R1} \\\
-p ../for_dada2/${R2} \\\
--discard-untrimmed \\\
-j 0 \\\
${R1} ${R2} 1>> '../logs/cutadapt_trim_report.txt'
done
")

cat(scripttext, 
    file = paste0(CODE_LOCATION, "/trim_primers.sh")) #note, R code, not bash code, so no shQuote()
