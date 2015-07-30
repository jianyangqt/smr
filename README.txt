
####SMR#######

# dense binary ESD file
./SMR --bfile aricEur_bg_hwe_c1 --gwas-summary height_step1_refale.raw --beqtl-summary blood_eqtl --out test

# dense text ESD file
./SMR --bfile aricEur_bg_hwe_c1 --gwas-summary height_step1_refale.raw --eqtl-summary blood_eqtl --out test

# sparse binary ESD file
./SMR --bfile aricEur_bg_hwe_c1 --gwas-summary height_step1_refale.raw --beqtl-summary blood_eqtl.s --out test

# with --keep   and --extract-snp parameters (—keep can only affect Plink data, —extract-snp can affect both Plink data and eQTL data)
./SMR --bfile aricEur_bg_hwe_c1 --gwas-summary height_step1_refale.raw --beqtl-summary blood_eqtl --out test --keep indi.list --extract-snp snp.list


# with --extract-probe parameter
./SMR --bfile aricEur_bg_hwe_c1 --gwas-summary height_step1_refale.raw --beqtl-summary blood_eqtl --out test --extract-probe probe.list


# with --maf parameter
./SMR --bfile aricEur_bg_hwe_c1 --gwas-summary height_step1_refale.raw --beqtl-summary blood_eqtl --out test --maf 0.01

# with all the computing parameters
./SMR --bfile aricEur_bg_hwe_c1 --gwas-summary height_step1_refale.raw --beqtl-summary blood_eqtl --out test --keep indi.list --extract-snp snp.list --maf 0.01 --extract-probe probe.list



####Make binary or text eQTL data files############

# convert from besd to esd
./SMR --beqtl-summary blood_eqtl --make-esd --out test


# convert from esd to besd
./SMR --eqtl-summary blood_eqtl --make-besd --out test


# extract from eQTL data according to a specified SNP list
#./SMR --beqtl-summary blood_eqtl --make-esd --out test --extract-snp snp.list
./SMR --beqtl-summary blood_eqtl --make-besd --out test --extract-snp snp.list


# extract from eQTL data according to a specified probe list
#./SMR --beqtl-summary blood_eqtl --make-esd --out test --extract-probe probe.list
./SMR --beqtl-summary blood_eqtl --make-besd --out test --extract-probe probe.list


# extract and specify both SNP list and probe list
#./SMR --beqtl-summary blood_eqtl --make-esd --out test --extract-probe probe.list --extract-snp snp.list
./SMR --beqtl-summary blood_eqtl --make-besd --out test --extract-probe probe.list --extract-snp snp.list


# extract common SNPs between GWAS data and eQTL data then save in binary (or text) format
#./SMR --gwas-summary height_step1_refale.raw --beqtl-summary blood_eqtl --out test --make-esd
./SMR --gwas-summary height_step1_refale.raw --beqtl-summary blood_eqtl --out test --make-besd 


# extract common SNPs between Plink data and eQTL data then save in binary (or text) format
#./SMR --bfile aricEur_bg_hwe_c1 --beqtl-summary blood_eqtl --out test --make-esd
./SMR --bfile aricEur_bg_hwe_c1 --beqtl-summary blood_eqtl --out test --make-besd


# extract common SNPs among Plink data, GWAS data and eQTL data then save in binary (or text) format
#./SMR --bfile aricEur_bg_hwe_c1 --gwas-summary height_step1_refale.raw --beqtl-summary blood_eqtl --out test --make-esd
./SMR --bfile aricEur_bg_hwe_c1 --gwas-summary height_step1_refale.raw --beqtl-summary blood_eqtl --out test --make-besd

# extract common SNPs with specified parameters
./SMR --bfile aricEur_bg_hwe_c1 --gwas-summary height_step1_refale.raw --beqtl-summary blood_eqtl --out test --keep indi.list --extract-snp snp.list --maf 0.01 --extract-probe probe.list --make-besd

# --make-besd and --make-esd can be combined together
./SMR --bfile aricEur_bg_hwe_c1 --gwas-summary height_step1_refale.raw --beqtl-summary blood_eqtl --out test --keep indi.list --extract-snp snp.list --maf 0.01 --extract-probe probe.list --make-besd --make-esd

