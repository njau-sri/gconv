#include <fstream>
#include <iostream>
#include <algorithm>
#include "pedio.h"
#include "split.h"

using std::size_t;

namespace
{
    int parse_ped_gt(const char *s, size_t n, char &a)
    {
        if (n != 1) {
            std::cerr << "ERROR: PED genotype must be represented by 1 character: " << std::string(s, n) << "\n";
            return 1;
        }

        a = s[0];

             if (a == '1') a = 'A';
        else if (a == '2') a = 'C';
        else if (a == '3') a = 'G';
        else if (a == '4') a = 'T';
        else if (a == '0') a = 'N';
        else if (a != 'A' && a != 'C' && a != 'G' && a != 'T') {
            std::cerr << "ERROR: invalid PED genotype code: " << std::string(s, n) << "\n";
            return 2;
        }

        return 0;
    }

    int check_compatibility_ped(const Genotype &gt)
    {
        for (auto &v : gt.allele) {
            if (v.size() > 2)
                return 1;

            for (auto &e : v) {
                if (e != "A" || e != "C" || e != "G" || e != "T")
                    return 2;
            }
        }

        return 0;
    }

} // namespace

int parse_ped_entry(const std::string &s, PedEntry &e)
{
    std::vector<Token> v;
    split(s, " \t", v);

    if (v.size() < 6) {
        std::cerr << "ERROR: incorrect number of columns at PED entry line: " << v.size() << "\n";
        return 1;
    }

    e.fid = v[0].to_string();
    e.iid = v[1].to_string();
    e.pid = v[2].to_string();
    e.mid = v[3].to_string();
    e.sex = std::stoi(v[4].to_string());
    e.pheno = std::stod(v[5].to_string());

    e.gt.clear();
    for (auto itr = v.begin() + 6; itr != v.end(); ++itr) {
        char a;
        if (parse_ped_gt(itr->data(), itr->size(), a) != 0)
            return 1;
        e.gt.push_back(a);
    }

    return 0;
}

int parse_map_entry(const std::string &s, MapEntry &e)
{
    std::vector<std::string> v;
    split(s, " \t", v);

    if (v.size() != 4) {
        std::cerr << "ERROR: expected 4 columns at MAP entry line: " << v.size() << "\n";
        return 1;
    }

    e.chr = v[0];
    e.id = v[1];
    e.dist = std::stod(v[2]);
    e.pos = std::stoi(v[3]);

    return 0;
}

int read_ped(const std::string &filename, Genotype &gt)
{
    std::ifstream ifs1(filename + ".map");
    if (!ifs1) {
        std::cerr << "ERROR: can't open file for reading: " << filename << ".map\n";
        return 1;
    }

    MapEntry me;

    for (std::string line; std::getline(ifs1, line); ) {
        line.erase(line.find_last_not_of("\r\n") + 1);

        if (parse_map_entry(line, me) != 0)
            return 1;

        gt.loc.push_back(me.id);
        gt.chr.push_back(me.chr);
        gt.pos.push_back(me.pos);
    }

    std::ifstream ifs2(filename + ".ped");
    if (!ifs2) {
        std::cerr << "ERROR: can't open file for reading: " << filename << ".ped\n";
        return 1;
    }

    PedEntry pe;
    std::vector< std::vector<allele_t> > dat;

    for (std::string line; std::getline(ifs2, line); ) {
        line.erase(line.find_last_not_of("\r\n") + 1);

        if (parse_ped_entry(line, pe) != 0)
            return 1;

        if (pe.gt.size() != 2 * gt.loc.size()) {
            std::cerr << "ERROR: column count doesn't match at " << pe.fid << ", " << pe.iid << "\n";
            return 1;
        }

        dat.push_back(pe.gt);
    }

    auto m = gt.loc.size();
    auto n = dat.size();
    std::vector<std::string> as;

    for (size_t j = 0; j < m; ++j) {
        std::vector<allele_t> v;
        for (size_t i = 0; i < n; ++i) {
            v.push_back(dat[i][j * 2]);
            v.push_back(dat[i][j * 2 + 1]);
        }

        auto z = v;
        std::sort(z.begin(), z.end());
        z.erase(std::unique(z.begin(), std::remove(z.begin(), z.end(), 'N')), z.end());

        if (z.size() > 2) {
            std::cerr << "ERROR: PED variant must be bi-allelic: " << z[0];
            for (size_t i = 1; i < z.size(); ++i)
                std::cerr << "/" << z[i];
            std::cerr << "\n";
            return 1;
        }

        as.clear();
        allele_t mis = 'N', ref = 'N', alt = 'N';
        if (!z.empty()) {
            ref = z[0];
            as.emplace_back(1, z[0]);
        }
        if (z.size() == 2) {
            alt = z[1];
            as.emplace_back(1, z[1]);
        }

        for (auto &a : v) {
            if (a == mis)
                a = 0;
            else if (a == ref)
                a = 1;
            else if (a == alt)
                a = 2;
        }

        gt.dat.push_back(v);
        gt.allele.push_back(as);
    }

    gt.ploidy = 2;

    return 0;
}

int write_ped(const Genotype &gt, const std::string &filename)
{
    int info = check_compatibility_ped(gt);
    if (info != 0) {
        std::cerr << "ERROR: genotype data is not compatible with PED format: " << info << "\n";
        return 1;
    }

    std::ofstream ofs1(filename + ".ped");
    if (!ofs1) {
        std::cerr << "ERROR: can't open file for writing: " << filename << ".ped\n";
        return 1;
    }

    std::ofstream ofs2(filename + ".map");
    if (!ofs2) {
        std::cerr << "ERROR: can't open file for writing: " << filename << ".map\n";
        return 1;
    }

    auto m = gt.loc.size();
    auto n = gt.ind.size();
    bool haploid = gt.ploidy == 1;

    std::string line;

    for (size_t i = 0; i < n; ++i) {
        line.clear();
        line.append(gt.ind[i]).append(" 0 0 1 1");
        auto k1 = haploid ? i : i * 2;
        auto k2 = haploid ? i : i * 2 + 1;
        for (size_t j = 0; j < m; ++j) {
            auto a = gt.dat[j][k1];
            auto b = gt.dat[j][k2];
            if (a && b) {
                line.push_back(' ');
                line.append(gt.allele[j][a - 1]);
                line.push_back(' ');
                line.append(gt.allele[j][b - 1]);
            }
            else
                line.append(" 0 0");
        }
        ofs1 << i + 1 << " " << line << "\n";
    }

    for (size_t j = 0; j < m; ++j)
        ofs2 << gt.chr[j] << " " << gt.loc[j] << " 0 " << gt.pos[j] << "\n";

    return 0;
}
