
#ifndef _GRID2_HAS_BEEN_INCLUDED_
#define _GRID2_HAS_BEEN_INCLUDED_

#include "bbox2.h"
#include "index2.h"

class Grid2
{
public:
    Grid2()
        : domain( Vec2(-1,-1), Vec2(1,1) ),
          size( 100,100 )
    {}

    Grid2( const Bbox2& domain, const Index2& size )
    {
        init( domain, size );
    }

    void init( const Bbox2& domain, const Index2& size )
    {
        this->domain = domain;
        this->size = size;
    }

    Bbox2 getDomain( void ) const
    {
        return domain;
    }

    Index2 getSize() const
    {
        return size;
    }

    Vec2 getCellDx( void ) const
    {
        const Vec2 extent( domain.maxPosition - domain.minPosition );
        return Vec2( extent.x / float( size.x ), extent.y / float( size.y ) );
    }

    Index2 getSizeFacesX() const
    {
        return Index2( size.x + 1, size.y );
    }

    Index2 getSizeFacesY() const
    {
        return Index2( size.x, size.y + 1 );
    }

    Index2 getSizeNodes() const
    {
        return Index2( size.x + 1, size.y + 1 );
    }

    Vec2 getNodePos( const Index2& id ) const
    {
        const Vec2 pos(
            domain.minPosition + Vec2( float(id.x), float(id.y) ) * getCellDx() );
        return pos;
    }

    Vec2 getCellPos( const Index2& id ) const
    {
        const Vec2 pos(
            domain.minPosition + Vec2( float(id.x) + 0.5f, float(id.y) + 0.5f ) * getCellDx() );
        return pos;
    }

    Vec2 getFaceXPos( const Index2& id ) const
    {
        const Vec2 pos(
            domain.minPosition + Vec2( float(id.x), float(id.y) + 0.5f ) * getCellDx() );
        return pos;
    }

    Vec2 getFaceYPos( const Index2& id ) const
    {
        const Vec2 pos(
            domain.minPosition + Vec2( float(id.x) + 0.5f, float(id.y) ) * getCellDx() );
        return pos;
    }
 
    Vec2 getCellIndex( const Vec2& pos ) const
    {
        const Vec2 ispos( ( pos - domain.minPosition ) / getCellDx() - 0.5f );
        return ispos;
    }

    Vec2 getFaceIndex( const Vec2& pos, const unsigned int axis ) const
    {
        Vec2 ispos( ( pos - domain.minPosition ) / getCellDx() );
        if( axis == 0 ) ispos.y -= 0.5f;
        if( axis == 1 ) ispos.x -= 0.5f;
        return ispos;
    }

private:
    Bbox2 domain;
    Index2 size;
};

#endif
