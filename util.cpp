#include "util.h"


bool starts_with(const std::string &s1, const std::string &s2)
{
    return s1.size() >= s2.size() && s1.compare(0, s2.size(), s2) == 0;
}

bool ends_with(const std::string &s1, const std::string &s2)
{
    return s1.size() >= s2.size() && s1.compare(s1.size() - s2.size(), s2.size(), s2) == 0;
}

std::string join(const std::vector<std::string> &vs, const std::string &sep)
{
    std::string s;

    auto itr = vs.begin();
    if (itr != vs.end()) {
        s += *itr;
        ++itr;
    }

    for (; itr != vs.end(); ++itr) {
        s += sep;
        s += *itr;
    }

    return s;
}
