#pragma once

#include <string>
#include <vector>
#include "vcfio.h"

struct PedEntry
{
    std::string fid;
    std::string iid;
    std::string pid;
    std::string mid;
    std::vector<allele_t> gt;
    double pheno = 0;
    int sex = 0;
};

struct MapEntry
{
    std::string chr;
    std::string id;
    double dist;
    int pos;
};

int parse_ped_entry(const std::string &s, PedEntry &e);

int parse_map_entry(const std::string &s, MapEntry &e);

int read_ped(const std::string &filename, Genotype &gt);

int write_ped(const Genotype &gt, const std::string &filename);
