#!/bin/bash
#SBATCH -p batch
#SBATCH --mem=10G
#SBATCH -t 48:00:00
#SBATCH -n 1
#SBATCH -N 1

##cluster specific module load
##this script only uses bcftools
module load bcftools

VCF="data/filtered_snps.pass.vcf"
VAR_TABLE="data/filtered_var.nohet.table"

FRQ_FILE="data/joined.sync.MAF01.frq"
OUTFILE="data/var_frq.nohet.tsv"

# =============================================================================
# STEP 1: VCF to genotype table
# =============================================================================
# Extract CHROM, POS, REF, ALT and per-sample GT from VCF.
# Homozygous REF (0/0) -> 1, homozygous ALT -> 0, heterozygous or missing -> NA.

printf "CHROM\tPOS\tREF\tALT" > "$VAR_TABLE"
bcftools query -l "$VCF" | while read S; do
    printf "\t%s" "$S" >> "$VAR_TABLE"
done
printf "\n" >> "$VAR_TABLE"

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$VCF" \
| awk -F'\t' 'BEGIN{OFS="\t"}
NR==1 {
    next
}
{
    n = split($0, a, FS)

    for(i=5; i<=n; i++){
        gt = a[i]

        if(gt=="./." || gt==".|." || gt=="." || gt==""){
            a[i] = "NA"
            continue
        }

        gsub(/\|/, "/", gt)
        split(gt, al, "/")

        if(al[1]=="." || al[2]=="."){
            a[i] = "NA"
            continue
        }

        if(al[1]=="0" && al[2]=="0"){
            a[i] = 1
        } else if(al[1] != al[2]){
            a[i] = "NA"
        } else {
            a[i] = 0
        }
    }

    for(i=1; i<=n; i++){
        printf "%s", a[i]
        if(i<n) printf OFS; else printf ORS
    }
}' >> "$VAR_TABLE"

# =============================================================================
# STEP 2: Join genotype table with pooled allele frequencies
# =============================================================================
# Match rows on CHROM+POS and append allele frequency columns to the right.

awk 'BEGIN{FS=OFS="\t"}

    FNR==1 {
        if (NR==1) {
            header1=$0
        } else {
            n=split(header1, h_fields, FS)
            frq_header=""
            for(i=3; i<=n; i++) frq_header = frq_header (i>3 ? FS : "") h_fields[i]
            print $0, frq_header > outfile
        }
        next
    }

    NR==FNR {
        key=$1 FS $2
        frq[key]=$0
        next
    }

    {
        key=$1 FS $2
        if (key in frq) {
            n=split(frq[key], frq_fields, FS)
            frq_tail=""
            for(i=3; i<=n; i++) frq_tail = frq_tail (i>3 ? FS : "") frq_fields[i]
            print $0, frq_tail > outfile
        }
    }

' outfile="$OUTFILE" "$FRQ_FILE" "$VAR_TABLE"
