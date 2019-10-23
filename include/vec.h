#include <vector>
#include "solution.h"
class vec
{
    std::vector<Solution*> _items;
public:
    template<class T>
    vec(std::vector<T> v)
        : _items(v.begin(), v.end())
    {}
    // other stuff, see rule of 5
};
