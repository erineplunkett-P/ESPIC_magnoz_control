#include "parse.h"
#include "cstdio"

bool parse_line(std::vector<std::string> &word, FILE *fp)
{
    if (!word.empty())
        word.clear();

    char delim[8] = " \t\n\r";
    std::stringstream ss;
    std::string s;

    while (1)
    {
        char line[MaxLine];
        if (0 == fgets(line, MaxLine, fp))
        {
            return false;
        } // file ends

        ss.clear();
        ss.str(std::string(line));
        std::getline(ss, s, '!');
        if (s.empty())
        {
            return true;
        } // this is a comment line

        std::vector<std::string> one = split(s, delim);
        if (one.empty())
        {
            return true;
        } // this is a blank line

        if (one.back() == "&")
        { // catenate next line
            word.insert(word.end(), one.begin(), one.end() - 1);
            //      for (std::vector<std::string>::const_iterator it = one.begin(); it != one.end()-1; ++it)
            //        word.push_back(*it);
            continue;
        }
        else
        {
            word.insert(word.end(), one.begin(), one.end());
            //      for (std::vector<std::string>::const_iterator it = one.begin(); it != one.end(); ++it)
            //        word.push_back(*it);
            break;
        }
    } // while (1)

    for (std::vector<std::string>::const_iterator it = word.begin(); it != word.end(); ++it)
    {
        printf("%s ", it->c_str());
    }
    printf("\n");

    return true;
}
