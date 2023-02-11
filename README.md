# intragenic-rearrangement-burden
calculation of intragenic rearrangement burden based on structural variant calling files

## usage
Users need to type in four arguments (1. input vcf file; 2. sample name; 3. output IGR count file; 4. exon database) for the CountIGRBurden.exe to run the calculation. 

./CountIGRBurden.exe input.vcf sample_name output exondb
  
## Example
./CountIGRBurden.exe SP78921.e9d607e8-41c3-4210-870b-a9ab9b1d1c8c.broad-dRanger_snowman.svfix.20151023.somatic.sv.vcf SP78921 SP78921.IGRcount.tsv UCSCAndGencodeCompV27lift37.annotation.gtf.exon_intron.merge.an
