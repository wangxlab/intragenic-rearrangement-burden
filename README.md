# intragenic-rearrangement-burden
Calculation of intragenic rearrangement burden based on structural variant calling files

# vcf file specifications
Required: vcf file must be in GRCh37 genome build.  
Required: vcf file must contain the following columns #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOUR.  
Required: vcf file ALT column must contain rearrangement annotation like "A]1:193674914]".  
Required: vcf file must contain "PASS" in the FILTER column.  

## usage for CountIGRBurden.exe
Users need to type in four arguments for the CountIGRBurden.exe to run the calculation:  
1. input vcf file match the specified format above;  
2. sample name   
3. output IGR count file  
4. exon database provided with this package (UCSCAndGencodeCompV27lift37.annotation.gtf.exon_intron.merge.an)   

```
./CountIGRBurden.exe input.vcf sample_name output exondb
```

## Example
```
./CountIGRBurden.exe example.sv.vcf example example.IGRcount.tsv UCSCAndGencodeCompV27lift37.annotation.gtf.exon_intron.merge.an
```

## Example Output
SampleID | Intergenic.AGR(n) | Intergenic.Inter-chr(n) | Intergenic.Intro-chr(n) | Intragenic.DEL(n) | Intragenic.DUP(n) | IGR.count(DUP+DEL) | IGR.burden(sqrt)
--- | --- | --- | --- | --- | --- | --- | --- 
example | 1 | 74 | 6 | 8 | 113 | 121 | 11

Notes to headers:  
Intergenic.AGR(n): adjacent intergenic rearrangements.  
Intergenic.Inter-chr(n): interchrosomal intergenic rearrangements.  
Intergenic.Intro-chr(n): intrachrosomal intergenic rearrangements.  
Intragenic.DEL(n): intragenic deletion counts.  
Intragenic.DUP(n): intragenic duplication counts.  
IGR.count(DUP+DEL): IGR counts including intragenic duplications and deletions  
IGR.burden(sqrt): square root transformed IGR counts.  
