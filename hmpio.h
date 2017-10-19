#pragma once

#include <string>
#include <vector>
#include "vcfio.h"

struct HmpEntry
{
    std::string chr;
    std::string id;
    std::vector<std::string> as;
    std::vector<allele_t> gt;
    int pos = -1;
};

int parse_hmp_header(const std::string &s, std::vector<std::string> &v);

int parse_hmp_entry(const std::string &s, HmpEntry &e);

int read_hmp(const std::string &filename, Genotype &gt);

int write_hmp(const Genotype &gt, const std::string &filename);
