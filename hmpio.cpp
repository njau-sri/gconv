#include <fstream>
#include <iostream>
#include <algorithm>
#include "hmpio.h"
#include "split.h"

using std::size_t;

namespace
{
    int parse_hmp_gt(const char *s, size_t n, char &a, char &b)
    {
        if (n != 2) {
            std::cerr << "ERROR: HapMap genotype must be represented by 2 characters: " << std::string(s, n) << "\n";
            return 1;
        }

        a = s[0];
        b = s[1];

        if ((a != 'A' && a != 'C' && a != 'G' && a != 'T' && a != 'I' && a != 'D' && a != 'N') ||
            (b != 'A' && b != 'C' && b != 'G' && b != 'T' && b != 'I' && b != 'D' && b != 'N')) {
            std::cerr << "ERROR: invalid HapMap genotype code: " << std::string(s, n) << "\n";
            return 2;
        }

        return 0;
    }

    int check_compatibility_hmp(const Genotype &gt)
    {
        for (auto &v : gt.allele) {
            if (v.size() > 2)
                return 1;

            if (v.size() == 2 && v[0].size() > 1 && v[1].size() > 1)
                return 2;

            for (auto &e : v) {
                if (e.size() == 1 && e != "A" && e != "C" && e != "G" && e != "T" && e != "-")
                    return 3;

                if (e.size() != 1 && e.find_first_not_of("ACGT") != std::string::npos)
                    return 4;
            }
        }

        return 0;
    }

} // namespace

int parse_hmp_header(const std::string &s, std::vector<std::string> &v)
{
    v.clear();

    split(s, " \t", v);

    if (v.size() < 11) {
        std::cerr << "ERROR: incorrect number of columns in HapMap header line: " << v.size() << "\n";
        return 1;
    }

    v.erase(v.begin(), v.begin() + 11);

    return 0;
}

int parse_hmp_entry(const std::string &s, HmpEntry &e)
{
    std::vector<Token> v;
    split(s, " \t", v);

    if (v.size() < 11) {
        std::cerr << "ERROR: incorrect number of columns at HapMap entry line: " << v.size() << "\n";
        return 1;
    }

    e.id = v[0].to_string();
    e.chr = v[2].to_string();
    e.pos = std::stoi(v[3].to_string());

    e.gt.clear();
    for (auto itr = v.begin() + 11; itr != v.end(); ++itr) {
        char a, b;
        if (parse_hmp_gt(itr->data(), itr->size(), a, b) != 0)
            return 1;
        e.gt.push_back(a);
        e.gt.push_back(b);
    }

    e.as.clear();
    split(v[1].to_string(), "/", e.as);

    if (e.as.size() != 2 || e.as[0] == "N" || e.as[1] == "N") {
        auto z = e.gt;
        std::sort(z.begin(), z.end());
        z.erase(std::unique(z.begin(), std::remove(z.begin(), z.end(), 'N')), z.end());

        if (z.size() > 2) {
            std::cerr << "ERROR: HapMap variant must be bi-allelic: " << z[0];
            for (size_t i = 1; i < z.size(); ++i)
                std::cerr << "/" << z[i];
            std::cerr << "\n";
            return 1;
        }

        e.as.clear();
        if (!z.empty())
            e.as.emplace_back(1, z[0]);
        if (z.size() == 2)
            e.as.emplace_back(1, z[1]);
    }

    if (e.as.empty()) {
        std::fill(e.gt.begin(), e.gt.end(), 0);
        return 0;
    }

    if (e.as.size() == 1) {
        std::fill(e.gt.begin(), e.gt.end(), 1);
        return 0;
    }

    allele_t ref = 'N', alt = 'N';
    if (e.as[0] == "A" || e.as[0] == "C" || e.as[0] == "G" || e.as[0] == "T")
        ref = e.as[0][0];
    if (e.as[1] == "A" || e.as[1] == "C" || e.as[1] == "G" || e.as[1] == "T")
        alt = e.as[1][0];

    static const allele_t mis = 'N', ins = 'I', del = 'D';

    for (auto &a : e.gt) {
        if (a == mis)
            a = 0;
        else if (a == ref)
            a = 1;
        else if (a == alt)
            a = 2;
        else if (a == ins)
            a = e.as[0] == "-" ? 2 : 1;
        else if (a == del)
            a = e.as[0] == "-" ? 1 : 2;
        else {
            std::cerr << "ERROR: inconsistent allele code: " << e.id << ", " << v[1].to_string() << ", " << a << "\n";
            return 1;
        }
    }

    return 0;
}

int read_hmp(const std::string &filename, Genotype &gt)
{
    std::ifstream ifs(filename);
    if (!ifs) {
        std::cerr << "ERROR: can't open file for reading: " << filename << "\n";
        return 1;
    }

    for (std::string line; std::getline(ifs, line);) {
        line.erase(line.find_last_not_of("\r\n") + 1);
        if (parse_hmp_header(line, gt.ind) != 0)
            return 1;
        break;
    }

    HmpEntry e;

    for (std::string line; std::getline(ifs, line);) {
        line.erase(line.find_last_not_of("\r\n") + 1);

        if (parse_hmp_entry(line, e) != 0)
            return 1;

        if (e.gt.size() != 2 * gt.ind.size()) {
            std::cerr << "ERROR: column count doesn't match at " << e.id << "\n";
            return 1;
        }

        gt.loc.push_back(e.id);
        gt.chr.push_back(e.chr);
        gt.pos.push_back(e.pos);
        gt.allele.push_back(e.as);
        gt.dat.push_back(e.gt);
    }

    gt.ploidy = 2;

    return 0;
}

int write_hmp(const Genotype &gt, const std::string &filename)
{
    int info = check_compatibility_hmp(gt);
    if (info != 0) {
        std::cerr << "ERROR: genotype data is not compatible with HapMap format: " << info << "\n";
        return 1;
    }

    std::ofstream ofs(filename);
    if (!ofs) {
        std::cerr << "ERROR: can't open file for writing: " << filename << "\n";
        return 1;
    }

    auto m = gt.loc.size();
    auto n = gt.ind.size();
    bool haploid = gt.ploidy == 1;

    std::vector<bool> indel(m, false);
    for (size_t j = 0; j < m; ++j) {
        auto& as = gt.allele[j];
        if (as.size() == 1 && (as[0] == "-" || as[0].size() > 1))
            indel[j] = true;
        if (as.size() == 2 && (as[0] == "-" || as[1] == "-"))
            indel[j] = true;
    }

    ofs << "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode";
    for (size_t i = 0; i < n; ++i)
        ofs << "\t" << gt.ind[i];
    ofs << "\n";

    std::string line;

    for (size_t j = 0; j < m; ++j) {
        line.clear();
        line.append("NA\tNA\tNA\tNA\tNA\tNA\tNA");

        if (indel[j]) {
            auto ref = "D", alt = "I";
            if (gt.allele[j][0] != "-")
                std::swap(ref, alt);
            for (size_t i = 0; i < n; ++i) {
                auto a = haploid ? gt.dat[j][i] : gt.dat[j][i * 2];
                auto b = haploid ? gt.dat[j][i] : gt.dat[j][i * 2 + 1];
                if (a && b) {
                    line.push_back('\t');
                    line.append(a == 1 ? ref : alt).append(b == 1 ? ref : alt);
                }
                else
                    line.append("\tNN");
            }
        }
        else {
            for (size_t i = 0; i < n; ++i) {
                auto a = haploid ? gt.dat[j][i] : gt.dat[j][i * 2];
                auto b = haploid ? gt.dat[j][i] : gt.dat[j][i * 2 + 1];
                if (a && b) {
                    line.push_back('\t');
                    line.append(gt.allele[j][a - 1]).append(gt.allele[j][b - 1]);
                }
                else
                    line.append("\tNN");
            }
        }

        ofs << gt.loc[j] << "\t";

        if (gt.allele[j].empty())
            ofs << "N/N";
        else {
            ofs << gt.allele[j][0];
            if (gt.allele[j].size() == 1)
                ofs << "/N";
            else
                ofs << "/" << gt.allele[j][1];
        }

        ofs << "\t" << gt.chr[j] << "\t" << gt.pos[j] << "\t" << line << "\n";
    }

    return 0;
}
