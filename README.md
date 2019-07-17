# gconv

Genotype file format conversion between [Variant Call Format (VCF)](https://samtools.github.io/hts-specs/), [PED/MAP](http://zzz.bwh.harvard.edu/plink/data.shtml#ped), [HapMap](https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Load/Load#markdown-header-hapmap) and the general genotype file format (see below).

## Download

https://github.com/njau-sri/gconv/releases

### 中国下载镜像

腾讯微云分享链接 https://share.weiyun.com/5NKUT0T

## Command line options

```
usage: gconv [options]
  --geno  <>    Input general genotype file
  --hmp   <>    Input HapMap genotype file
  --out   <>    Output file with format suffix (.vcf/.ped/.hmp/.geno)
  --ped   <>    Input PLINK ped file (map file has same basename)
  --vcf   <>    Input VCF genotype file
  --sort        sorting loci in ascending chromosome position order
```

## General genotype file format (.geno)

Each row is a marker, each column is an individual. The first row contains column names and individual names. The first three columns are marker names, chromosome labels and genome positions, respectively.

- arbitrary allele code: `number` `character` `string`

- supported allele separator: `space` `tab` `/` `:`

- supported missing genotype: `N` `-` `.` `?`

### Example: two individuals typed at five SSR markers

| Locus | Chromosome | Position | Ind1    | Ind2    |
|-------|------------|----------|---------|---------|
| Mk1   | 1          | 100      | 630/630 | 909/909 |
| Mk2   | 1          | 200      | 557/711 | 711/711 |
| Mk3   | 1          | 300      | 445/445 | 445/668 |
| Mk4   | 2          | 1001     | 307/307 | 340/340 |
| Mk5   | 2          | 1002     | 264/273 | 264/264 |
