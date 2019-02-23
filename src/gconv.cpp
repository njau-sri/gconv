#include <string>
#include <numeric>
#include <iostream>
#include <algorithm>
#include "cmdline.h"
#include "vcf.h"
#include "ped.h"
#include "hmp.h"
#include "geno.h"
#include "util.h"


using std::size_t;


namespace {


struct Parameter
{
    std::string vcf;
    std::string ped;
    std::string hmp;
    std::string geno;
    std::string out;
    bool sort = false;
} par;


void sort_chrpos(Genotype &gt)
{
    auto chr = unique(gt.chr);

    auto m = gt.loc.size();
    std::vector<size_t> z;
    z.reserve(m);

    for (auto &e : chr) {
        std::vector<size_t> idx;
        for (size_t i = 0; i < m; ++i) {
            if (gt.chr[i] == e)
                idx.push_back(i);
        }

        auto pos = subset(gt.pos, idx);
        if ( ! std::is_sorted(pos.begin(), pos.end()) ) {
            auto ord = order(pos);
            subset(idx,ord).swap(idx);
        }

        z.insert(z.end(), idx.begin(), idx.end());
    }

    subset(gt.loc,z).swap(gt.loc);
    subset(gt.chr,z).swap(gt.chr);
    subset(gt.pos,z).swap(gt.pos);
    subset(gt.dat,z).swap(gt.dat);
    subset(gt.allele,z).swap(gt.allele);
}


} // namespace


int gconv(int argc, char *argv[])
{
    std::cerr << "GCONV 2019.0.dev (Built on " __DATE__ " " __TIME__ ")\n";

    CmdLine cmd;

    cmd.add("--vcf", "VCF genotype file", "");
    cmd.add("--ped", "PLINK ped file (map file has same basename)", "");
    cmd.add("--hmp", "HapMap genotype file", "");
    cmd.add("--geno", "General genotype file", "");
    cmd.add("--out", "output file with format suffix (.vcf/.ped/.hmp/.geno)", "");
    cmd.add("--sort", "sorting loci in ascending chromosome position order");

    cmd.parse(argc, argv);

    if (argc < 2) {
        cmd.show();
        return 1;
    }

    par.vcf = cmd.get("--vcf");
    par.ped = cmd.get("--ped");
    par.hmp = cmd.get("--hmp");
    par.geno = cmd.get("--geno");
    par.out = cmd.get("--out");
    par.sort = cmd.has("--sort");

    Genotype gt;

    std::cerr << "INFO: reading genotype file...\n";

    if ( ! par.vcf.empty() ) {
        if (read_vcf(par.vcf, gt) != 0)
            return 1;
    }
    else if ( ! par.ped.empty() ) {
        auto prefix = par.ped;
        if (ends_with(prefix, ".ped"))
            prefix = prefix.substr(0, prefix.size() - 4);
        if (read_ped(prefix, gt) != 0)
            return 1;
    }
    else if ( ! par.hmp.empty() ) {
        if (read_hmp(par.hmp, gt) != 0)
            return 1;
    }
    else if ( ! par.geno.empty() ) {
        if (read_geno(par.geno, gt) != 0)
            return 1;
    }

    std::cerr << "INFO: " << gt.ind.size() << " individuals, " << gt.loc.size() << " loci\n";

    if (gt.ind.empty() && gt.loc.empty())
        return 1;

    if (par.sort)
        sort_chrpos(gt);

    if ( ends_with(par.out, ".vcf") ) {
        if (write_vcf(gt, par.out) != 0)
            return 1;
    }
    else if ( ends_with(par.out, ".ped") ) {
        auto prefix = par.out.substr(0, par.out.size() - 4);
        if (write_ped(gt, prefix) != 0)
            return 1;
    }
    else if ( ends_with(par.out, ".hmp") ) {
        if (write_hmp(gt, par.out) != 0)
            return 1;
    }
    else if ( ends_with(par.out, ".geno") ) {
        if (write_geno(gt, par.out) != 0)
            return 1;
    }
    else {
        std::cerr << "ERROR: unrecognized output format: " << par.out << "\n";
        return 1;
    }

    return 0;
}
