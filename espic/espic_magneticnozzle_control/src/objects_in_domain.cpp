#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include "objects_in_domain.h"
#include "parse.h"
#include "pic_io.h"

using std::cout;
using std::endl;

// destructor
ObjectList::~ObjectList()
{
    if (nobjs > 0)
    {
        nobjs = 0;
        delete[] obj;
        obj = NULL;
    }
}

// copy constructor
ObjectList::ObjectList(const ObjectList &objlist) : nobjs(objlist.nobjs)
{
    obj = new Object[objlist.nobjs];

    for (size_type i = 0; i < objlist.nobjs; i++)
    {
        for (int a = 0; a < 3; a++)
        {
            obj[i].xmin[a] = objlist.obj[i].xmin[a];
            obj[i].xmax[a] = objlist.obj[i].xmax[a];
        }
        obj[i].realObj = objlist.obj[i].realObj;
        obj[i].phi[0] = objlist.obj[i].phi[0];
        obj[i].phi[1] = objlist.obj[i].phi[1];
        obj[i].uniform_pot = objlist.obj[i].uniform_pot;
    }
}

// constructor -- read from file
ObjectList::ObjectList(const std::string &ifilename)
{
    FILE *fp = fopen(ifilename.c_str(), "r");
    if (NULL == fp)
    {
        esPIC_exit("Cannot open file [" + ifilename + "]");
    }

    std::vector<std::string> word;

    // first line -- # of objects
    parse_line(word, fp);
    nobjs = (int)atoi(word[0].c_str());

    obj = new Object[nobjs];

    size_type iobj;
    for (iobj = 0; iobj < nobjs; iobj++)
    {

        parse_line(word, fp);
        for (int i = 0; i < 3; ++i)
        {
            obj[iobj].xmin[i] = (double)atof(word[2 * i].c_str());
            obj[iobj].xmax[i] = (double)atof(word[2 * i + 1].c_str());
        }

        assert(obj[iobj].xmin[0] <= obj[iobj].xmax[0]);
        assert(obj[iobj].xmin[1] <= obj[iobj].xmax[1]);
        assert(obj[iobj].xmin[2] <= obj[iobj].xmax[2]);

        parse_line(word, fp);
        obj[iobj].realObj = (word[0] == "1");

        // line defining phi
        parse_line(word, fp);

        if (word.size() == 1)
        { // uniform phi
            obj[iobj].phi[0] = obj[iobj].phi[1] = (double)atof(word[0].c_str());
            obj[iobj].uniform_pot = true;
        }
        else if (word.size() == 4)
        {
            if (word[0] == "linear")
            {
                obj[iobj].uniform_pot = false;
                if (word[1] == "x")
                {
                    obj[iobj].lin_direction = 0;
                }
                else if (word[1] == "y")
                {
                    obj[iobj].lin_direction = 1;
                }
                else
                {
                    cout << "--> Error: illegal definition of object potential. <--" << endl;
                    exit(-1);
                }
                obj[iobj].phi[0] = (double)atof(word[2].c_str());
                obj[iobj].phi[1] = (double)atof(word[3].c_str());
            }
            else
            {
                cout << "--> Error: illegal definition of object potential. <--" << endl;
                exit(-1);
            }
        }
        else
        {
            cout << "--> Error: illegal definition of object potential. <--" << endl;
            exit(-1);
        }
        // cout << "iobj = " << iobj << ": phi = (" << obj[iobj].phi[0] << ", " << obj[iobj].phi[1] << ")\n";
    }

    assert(iobj == nobjs);
}

// overload operator []
Object &ObjectList::operator[](const size_type iobj)
{

#ifdef DEBUG
    assert(iobj >= 0 && iobj < nobjs);
#endif

    return obj[iobj];
}

const Object &ObjectList::operator[](const size_type iobj) const
{

#ifdef DEBUG
    assert(iobj >= 0 && iobj < nobjs);
#endif

    return obj[iobj];
}

double Object::potential(const double x[3]) const
{
    if (uniform_pot)
        return phi[0];

    // return a linear interpolation of phi along the defined direction
    int a = lin_direction;
    return phi[0] + (phi[1] - phi[0]) * (x[a] - xmin[a]) / (xmax[a] - xmin[a]);
}

bool is_inside_object(const double x[3], const Object &obj)
{
    bool inside = true;
    for (int a = 0; a < 3; a++)
    {
        double xx = x[a] - domain[2 * a];

        inside = inside && (xx >= (obj.xmin[a] - domain[2 * a]));
        if (!inside)
            break;

        inside = inside && (xx <= (obj.xmax[a] - domain[2 * a]));
        if (!inside)
            break;
    }

    return inside;
}
