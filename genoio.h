#pragma once

#include "vcfio.h"

// General Genotype File Format
//    - arbitrary allele code         number   character   string
//    - supported allele separator    WHITESPACE   /   :
//    - supported missing genotype    N   -   .   ?
//
// Example  2 Individuals  5 SSR Markers
//    Locus  Chromosome  Position   Ind1      Ind2
//    Mk1    1           100        630/630   909/909
//    Mk2    1           200        557/711   711/711
//    Mk3    1           300        445/445   445/668
//    Mk4    2           1001       307/307   340/340
//    Mk5    2           1002       264/273   264/264

int read_geno(const std::string &filename, Genotype &gt);

int write_geno(const Genotype &gt, const std::string &filename);
