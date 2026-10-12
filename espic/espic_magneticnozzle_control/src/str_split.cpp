#include "str_split.h"

void split(const std::string &s, char delim, std::vector<std::string> &elems)
{
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim))
    {
        if (!item.empty())
            elems.push_back(item);
    }
}

std::vector<std::string> split(const std::string &s, char delim)
{
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

void split(const std::string &s, char *delim, std::vector<std::string> &elems)
{
    char cstr[MaxLine];
    memcpy(cstr, s.data(), s.length());
    cstr[s.length()] = '\0';
    //  for (unsigned i = 0; i <= s.length(); ++i) cstr[i] = s[i];

    char *token = std::strtok(cstr, delim);
    while (token)
    {
        elems.push_back(std::string(token));
        token = std::strtok(NULL, delim);
    }
}

std::vector<std::string> split(const std::string &s, char *delim)
{
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}
