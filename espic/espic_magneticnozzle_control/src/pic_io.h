#ifndef ESPIC_IO_H
#define ESPIC_IO_H

#include <string>

extern void esPIC_exit(const std::string &error_info);
extern bool is_file_exist(const std::string &filename);

#endif
