#!/bin/bash

vcf=$1
vcf2=$2

grep "^##" ${vcf}
grep "^##contig" ${vcf2}
grep "^##FORMAT" ${vcf2}
grep "^#C" ${vcf}
grep -v "^#" ${vcf} | grep -v "SVLEN=NULL"