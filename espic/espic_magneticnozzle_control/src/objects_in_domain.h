#ifndef OBJECTS_INDOMAIN_H
#define OBJECTS_INDOMAIN_H

#include <string>
#include <fstream>
#include "global.h"

struct Object
{
    double xmin[3];
    double xmax[3];
    bool realObj;
    bool uniform_pot;
    int lin_direction;
    double phi[2];

    double potential(const double x[3]) const;
};

class ObjectList
{
public:
    typedef int size_type;

    /* public methods */
    // default constructor
    ObjectList() : obj(NULL), nobjs(0) {}

    // copy constructor
    ObjectList(const ObjectList &);

    // constructor --  n objects
    ObjectList(const int n) : nobjs(n) { obj = new Object[n]; }

    // constructor -- read from file
    ObjectList(const std::string &ifilename);

    // destructor
    ~ObjectList();

    // return # of objects of this list
    size_type size() const { return nobjs; }

    // overload operator []
    Object &operator[](const int);
    const Object &operator[](const int) const;

private:
    Object *obj;
    size_type nobjs;
};

/* non-member function related */
// is one point inside the object?
bool is_inside_object(const double x[3], const Object &obj);

#endif
