#include <string>
#include <numeric>
#include <iostream>
#include <algorithm>
#include "cmdline.h"
#include "vcfio.h"
#include "pedio.h"
#include "hmpio.h"

namespace
{
    struct Par
    {
        std::string vcf;
        std::string ped;
        std::string hmp;
        std::string out;
        bool sort = false;
    };

    Par par;

    bool ends_with(const std::string &s1, const std::string &s2)
    {
        return s1.size() >= s2.size() && s1.compare(s1.size() - s2.size(), s2.size(), s2) == 0;
    }

    template<typename T>
    std::vector<T> unique(std::vector<T> v)
    {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
        return v;
    }

    template<typename T>
    std::vector<size_t> order(const std::vector<T> &v)
    {
        std::vector<size_t> z(v.size());
        std::iota(z.begin(), z.end(), size_t(0));
        std::sort(z.begin(), z.end(), [&v](size_t i, size_t j) { return v[i] < v[j]; });
        return z;
    }

    template<typename T1, typename T2>
    std::vector<T1> subset(const std::vector<T1> &v, const std::vector<T2> &idx)
    {
        std::vector<T1> z;
        z.reserve(idx.size());
        for (auto i : idx)
            z.push_back(v[i]);
        return z;
    }

    void sort_chrpos(Genotype &gt)
    {
        auto chrid = unique(gt.chr);

        auto m = gt.loc.size();
        std::vector<size_t> z;
        z.reserve(m);

        for (auto &e : chrid) {
            std::vector<size_t> idx;
            for (size_t i = 0; i < m; ++i) {
                if (gt.chr[i] == e)
                    idx.push_back(i);
            }

            auto pos = subset(gt.pos, idx);
            if (!std::is_sorted(pos.begin(), pos.end())) {
                auto ord = order(pos);
                subset(idx, ord).swap(idx);
            }

            z.insert(z.end(), idx.begin(), idx.end());
        }

        subset(gt.loc, z).swap(gt.loc);
        subset(gt.chr, z).swap(gt.chr);
        subset(gt.pos, z).swap(gt.pos);
        subset(gt.dat, z).swap(gt.dat);
        subset(gt.allele, z).swap(gt.allele);
    }

} // namespace

int gconv(int argc, char *argv[])
{
    std::cerr << "GCONV (Built on " __DATE__ " " __TIME__ ")\n";

    CmdLine cmd("gconv [options]");

    cmd.add("--vcf", "VCF genotype file", "");
    cmd.add("--ped", "PLINK ped file (map file has same basename)", "");
    cmd.add("--hmp", "HapMap genotype file", "");
    cmd.add("--out", "output file with format suffix (.vcf/.ped/.hmp)", "");
    cmd.add("--sort", "sorting loci in ascending chromosome position order");

    if (argc < 2) {
        cmd.help();
        return 1;
    }

    cmd.parse(argc, argv);

    par.vcf = cmd.get("--vcf");
    par.ped = cmd.get("--ped");
    par.hmp = cmd.get("--hmp");
    par.out = cmd.get("--out");
    par.sort = cmd.has("--sort");

    Genotype gt;

    std::cerr << "INFO: reading genotype file...\n";

    if (!par.vcf.empty()) {
        if (read_vcf(par.vcf, gt) != 0)
            return 1;
    }
    else if (!par.ped.empty()) {
        auto prefix = par.ped;
        if (ends_with(prefix, ".ped"))
            prefix = prefix.substr(0, prefix.size() - 4);
        if (read_ped(prefix, gt) != 0)
            return 1;
    }
    else if (!par.hmp.empty()) {
        if (read_hmp(par.hmp, gt) != 0)
            return 1;
    }

    std::cerr << "INFO: " << gt.ind.size() << " individuals, " << gt.loc.size() << " loci\n";

    if (par.sort)
        sort_chrpos(gt);

    if (ends_with(par.out, ".vcf")) {
        if (write_vcf(gt, par.out) != 0)
            return 1;
    }
    else if (ends_with(par.out, ".ped")) {
        auto prefix = par.out.substr(0, par.out.size() - 4);
        if (write_ped(gt, prefix) != 0)
            return 1;
    }
    else if (ends_with(par.out, ".hmp")) {
        if (write_hmp(gt, par.out) != 0)
            return 1;
    }
    else {
        std::cerr << "ERROR: unrecognized output format: " << par.out << "\n";
        return 1;
    }

    return 0;
}
