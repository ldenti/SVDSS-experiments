# Truthsets
We got these callsets by combining truvari `tp-base.vcf` and `fn.vcf` computed on full genome (i.e., no `.bed` included):
```
cat tp-base.vcf <(grep -v "^#" fn.vcf) | bcftools sort -T. -Oz -o callset.vcf.gz
tabix -p vcf callset.vcf.gz
```
