
#ifndef _BBOX2_HAS_BEEN_INCLUDED_
#define _BBOX2_HAS_BEEN_INCLUDED_

#include "vec2.h"

class Bbox2
{
public:
    Vec2 minPosition;
    Vec2 maxPosition;

    Bbox2()
        : minPosition(0,0), maxPosition(0,0)
    {}

    Bbox2( const float& minx, const float& miny, const float& maxx, const float& maxy )
    {
        minPosition = Vec2( minx, miny );
        maxPosition = Vec2( maxx, maxy );
    }

    Bbox2( const Vec2& minpos, const Vec2& maxpos )
    {
        minPosition = minpos;
        maxPosition = maxpos;
    }
};

#endif
