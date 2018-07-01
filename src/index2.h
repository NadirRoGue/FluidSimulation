
#ifndef _INDEX2_HAS_BEEN_INCLUDED_
#define _INDEX2_HAS_BEEN_INCLUDED_

class Index2
{
public:
    unsigned int x, y;

    Index2()
        : x(0), y(0)
    {}

    Index2( unsigned int x, unsigned int y )
    {
        this->x = x;
        this->y = y;
    }

    unsigned int& operator[] ( unsigned int i )
    {
        return i == 0 ? x : y;
    }
    const unsigned int& operator[] ( unsigned int i ) const
    {
        return i == 0 ? x : y;
    }
};

bool operator== ( const Index2& lhs, const Index2& rhs );
bool operator!= ( const Index2& lhs, const Index2& rhs );

#endif
