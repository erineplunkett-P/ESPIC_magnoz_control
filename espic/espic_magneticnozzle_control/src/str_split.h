#include <string>
#include <sstream>
#include <vector>
#include <cstring>

const int MaxLine = 1024;

void split(const std::string &, char, std::vector<std::string> &);
std::vector<std::string> split(const std::string &, char);

void split(const std::string &, char *, std::vector<std::string> &);
std::vector<std::string> split(const std::string &, char *);
