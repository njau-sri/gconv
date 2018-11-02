#include <fstream>
#include <iostream>
#include <algorithm>
#include "geno.h"
#include "split.h"


using std::size_t;


namespace {

template<typename T1, typename T2>
size_t index(const std::vector<T1> &v, const T2 &a)
{
    size_t i = 0, n = v.size();

    for (i = 0; i < n; ++i) {
        if (v[i] == a)
            break;
    }

    return i;
}

// http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html
//
//   W --- A/T
//   S --- C/G
//   M --- A/C
//   K --- G/T
//   R --- A/G
//   Y --- C/T
//

char encode_iupac(char a, char b)
{
    if (a == b && (a == 'A' || a == 'C' || a == 'G' || a == 'T'))
        return a;

    if (a > b)
        std::swap(a, b);

    if (a == 'A') {
        if (b == 'C') return 'M';
        if (b == 'G') return 'R';
        if (b == 'T') return 'W';
    }

    if (a == 'C') {
        if (b == 'G') return 'S';
        if (b == 'T') return 'Y';
    }

    if (a == 'G' && b == 'T')
        return 'K';

    return 'N';
}

std::pair<char, char> decode_iupac(char c)
{
    char a = 0, b = 0;

    if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
        a = b = c;
    else if (c == 'W') {
        a = 'A';
        b = 'T';
    }
    else if (c == 'S') {
        a = 'C';
        b = 'G';
    }
    else if (c == 'M') {
        a = 'A';
        b = 'C';
    }
    else if (c == 'K') {
        a = 'G';
        b = 'T';
    }
    else if (c == 'R') {
        a = 'A';
        b = 'G';
    }
    else if (c == 'Y') {
        a = 'C';
        b = 'T';
    }

    return { a, b };
}

bool is_iupac(const std::vector< std::vector<allele_t> > &dat)
{
    bool ret = false;

    for (auto &v : dat) {
        for (auto &a : v) {
            switch (a) {
            case 'A': case 'C': case 'G': case 'T': case 0:
                break;
            case 'W': case 'S': case 'M': case 'K': case 'R': case 'Y':
                ret = true;
                break;
            default:
                return false;
            }
        }
    }

    return ret;
}

int check_compat_iupac(const Genotype &gt)
{
    for (auto &v : gt.allele) {
        for (auto &e : v) {
            if (e.size() != 1)
                return 1;
            auto a = e[0];
            if (a != 'A' && a != 'C' && a != 'G' && a != 'T')
                return 2;
        }
    }

    return 0;
}

int read_genotype_char(const std::string &filename, Genotype &gt)
{
    std::ifstream ifs(filename);
    if ( ! ifs ) {
        std::cerr << "ERROR: can't open file for reading: " << filename << "\n";
        return 1;
    }

    for (std::string line; std::getline(ifs, line); ) {
        if ( ! line.empty() && line.back() == '\r' )
            line.pop_back();

        std::vector<std::string> vs;
        split(line, " \t", vs);
        if ( vs.empty() )
            continue;

        if (vs.size() < 3) {
            std::cerr << "ERROR: expected at least 3 columns at the first line\n";
            return 1;
        }

        gt.ind.assign(vs.begin() + 3, vs.end());
        break;
    }

    size_t ploidy = 0;
    auto n = gt.ind.size();

    for (std::string line; std::getline(ifs, line); ) {
        if ( ! line.empty() && line.back() == '\r' )
            line.pop_back();

        std::vector<Token> vt;
        split(line, " \t/:", vt);
        if ( vt.empty() )
            continue;

        if (ploidy == 0) {
            ploidy = vt.size() > 3 + n ? (vt.size() - 3) / n : 1;
            if (ploidy > 2) {
                std::cerr << "ERROR: unsupported polyploidy (" << ploidy << ") genotype data: "
                          << vt[0].to_string() << "\n";
                return 1;
            }
        }

        if (vt.size() != 3 + ploidy * n) {
            std::cerr << "ERROR: column count doesn't match (" << vt.size() << " != "
                      << 3 + ploidy * n << "): " << vt[0].to_string() << "\n";
            return 1;
        }

        gt.loc.push_back(vt[0].to_string());
        gt.chr.push_back(vt[1].to_string());
        gt.pos.push_back(std::stoi(vt[2].to_string()));

        std::vector<allele_t> v;
        for (auto itr = vt.begin() + 3; itr != vt.end(); ++itr) {
            if (itr->size() == 1)
                v.push_back( static_cast<allele_t>( (*itr)[0] ) );
            else
                return -1;
        }

        for (auto &e : v) {
            if (e == 'N' || e == '-' || e == '.' || e == '?')
                e = 0;
        }

        gt.dat.push_back(v);
    }

    bool iupac = ploidy == 1 && is_iupac(gt.dat);

    for (auto &v : gt.dat) {
        if ( iupac ) {
            std::vector<allele_t> w;
            w.reserve(v.size() * 2);
            for (auto a : v) {
                auto p = decode_iupac( static_cast<char>(a) );
                w.push_back( static_cast<allele_t>( p.first ) );
                w.push_back( static_cast<allele_t>( p.second ) );
            }
            v.swap(w);
        }

        auto u = v;
        std::sort(u.begin(), u.end());
        u.erase(std::unique(u.begin(), std::remove(u.begin(), u.end(), allele_t(0))), u.end());

        std::vector<std::string> allele;
        for (auto a : u)
            allele.emplace_back(1, a);
        gt.allele.push_back(allele);

        for (auto &a : v)
            a = a == 0 ? 0 : static_cast<allele_t>( index(u,a) + 1 );
    }

    gt.ploidy = static_cast<int>(ploidy);

    return 0;
}

int read_genotype_string(const std::string &filename, Genotype &gt)
{
    std::ifstream ifs(filename);
    if ( ! ifs ) {
        std::cerr << "ERROR: can't open file for reading: " << filename << "\n";
        return 1;
    }

    for (std::string line; std::getline(ifs, line); ) {
        if ( ! line.empty() && line.back() == '\r' )
            line.pop_back();

        std::vector<std::string> vs;
        split(line, " \t", vs);
        if ( vs.empty() )
            continue;

        if (vs.size() < 3) {
            std::cerr << "ERROR: expected at least 3 columns at the first line\n";
            return 1;
        }

        gt.ind.assign(vs.begin() + 3, vs.end());
        break;
    }

    size_t ploidy = 0;
    auto n = gt.ind.size();
    const Token missing("?", 1);

    for (std::string line; std::getline(ifs, line); ) {
        if ( ! line.empty() && line.back() == '\r' )
            line.pop_back();

        std::vector<Token> vt;
        split(line, " \t/:", vt);
        if ( vt.empty() )
            continue;

        if (ploidy == 0) {
            ploidy = vt.size() > 3 + n ? (vt.size() - 3) / n : 1;
            if (ploidy > 2) {
                std::cerr << "ERROR: unsupported polyploidy (" << ploidy <<") genotype data: "
                          << vt[0].to_string() << "\n";
                return 1;
            }
        }

        if (vt.size() != 3 + ploidy * n) {
            std::cerr << "ERROR: column count doesn't match at line (" << vt.size() << " != "
                      << 3 + ploidy * n << "): " << vt[0].to_string() << "\n";
            return 1;
        }

        gt.loc.push_back(vt[0].to_string());
        gt.chr.push_back(vt[1].to_string());
        gt.pos.push_back(std::stoi(vt[2].to_string()));

        std::vector<Token> u(vt.begin() + 3, vt.end());
        std::sort(u.begin(), u.end());
        u.erase(std::unique(u.begin(), std::remove(u.begin(), u.end(), missing)), u.end());

        if (u.size() > std::numeric_limits<allele_t>::max()) {
            std::cerr << "ERROR: exceed the maximum number of alleles: " << u.size() << "\n";
            return 1;
        }

        std::vector<std::string> allele;
        for (auto &e : u)
            allele.push_back(e.to_string());
        gt.allele.push_back(allele);

        std::vector<allele_t> v;
        for (auto itr = vt.begin() + 3; itr != vt.end(); ++itr) {
            if (*itr == missing)
                v.push_back(0);
            else
                v.push_back( static_cast<allele_t>( index(u,*itr) + 1 ) );
        }

        gt.dat.push_back(v);
    }

    gt.ploidy = static_cast<int>(ploidy);

    return 0;
}

} // namespace


int read_geno(const std::string &filename, Genotype &gt)
{
    int info = read_genotype_char(filename, gt);

    if (info < 0) {
        gt.loc.clear();
        gt.chr.clear();
        gt.pos.clear();
        gt.dat.clear();
        gt.allele.clear();
        info = read_genotype_string(filename, gt);
    }

    return info;
}

int write_geno(const Genotype &gt, const std::string &filename)
{
    std::ofstream ofs(filename);
    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file for writing: " << filename << "\n";
        return 1;
    }

    auto m = gt.loc.size();
    auto n = gt.ind.size();
    bool haploid = gt.ploidy == 1;
    bool iupac = check_compat_iupac(gt) == 0;
    const std::string missing = iupac ? "N" : "?";

    std::string line;

    ofs << "Locus\tChromosome\tPosition";
    for (size_t i = 0; i < n; ++i)
        ofs << "\t" << gt.ind[i];
    ofs << "\n";

    for (size_t j = 0; j < m; ++j) {
        ofs << gt.loc[j] << "\t" << gt.chr[j] << "\t" << gt.pos[j];

        line.clear();

        for (size_t i = 0; i < n; ++i) {
            line.push_back('\t');
            if ( haploid ) {
                auto a = gt.dat[j][i];
                line.append(a == 0 ? missing : gt.allele[j][a-1]);
            }
            else {
                auto a = gt.dat[j][i*2];
                auto b = gt.dat[j][i*2+1];
                if ( iupac ) {
                    if (a == 0 || b == 0)
                        line.append(missing);
                    else
                        line.push_back(encode_iupac(gt.allele[j][a-1][0], gt.allele[j][b-1][0]));
                }
                else {
                    line.append(a == 0 ? missing : gt.allele[j][a-1]);
                    line.push_back('/');
                    line.append(b == 0 ? missing : gt.allele[j][b-1]);
                }
            }
        }

        ofs << line << "\n";
    }

    return 0;
}
