#include <iostream>
#include <string>
#include <sys/stat.h>
#include <cstdlib>

using namespace std;

void esPIC_exit(const string &erro_info)
{
    cerr << "Error: -> " << erro_info << ". <-\n(Program terminated ...)" << endl;
    exit(0);
}

bool is_file_exist(const string &filename)
{
    struct stat buffer;

    return (stat(filename.c_str(), &buffer) == 0);
}
