
#ifndef _ARRAY2_HAS_BEEN_INCLUDED_
#define _ARRAY2_HAS_BEEN_INCLUDED_

#include "index2.h"

template< class T >
class Array2
{
public:
    Array2()
        : ptr( NULL )
    {}

    Array2( int x, int y )
        : ptr( NULL )
    {
        resize( Index2( x, y ) )
    }

    Array2( const Index2& size )
        : ptr( NULL )
    {
        resize( size );
    }

    Array2( const Array2< T >& array )
        : ptr( NULL )
    {
        copy( array );
    }

    ~Array2()
    {
        if( ptr ) delete[] ptr;
    }

    const Index2& getSize() const
    {
        return size;
    }

    const T* getData() const
    {
        return ptr;
    }

    unsigned int getLinearIndex( unsigned int i, unsigned int j ) const
    {
        return j * size.x + i;
    }

    bool resize( const Index2& size_ )
    {
        if( size != size_ )
        {
            if( ptr ) delete[] ptr;

            size = size_;
            ptr = new T[ size.x * size.y ];

            clear();

            return true;
        }
        return false;
    }

    void copy( const Array2< T >& src )
    {
        resize( src.size );

        for( unsigned int i = 0, n = size.x * size.y; i < n; ++i )
            ptr[ i ] = src[ i ];
    }

    void operator=( const Array2< T >& src )
    {
        copy( src );
    }

    const T& getValue( const Index2& id ) const
    {
        const unsigned int idx = getLinearIndex( id.x, id.y );
        return ptr[ idx ];
    }
    const T& getValue( const unsigned int i, const unsigned int j ) const
    {
        const unsigned int idx = getLinearIndex( i, j );
        return ptr[ idx ];
    }

    void setValue( const Index2& id, const T& value )
    {
        const unsigned int idx = getLinearIndex( id.x, id.y );
        ptr[ idx ] = value;
    }
    void setValue( const unsigned int i, const unsigned int j, const T& value )
    {
        const unsigned int idx = getLinearIndex( i, j );
        ptr[ idx ] = value;
    }

    const T& operator[] ( unsigned int i ) const
    {
        return ptr[ i ];
    }
    T& operator[] ( unsigned int i )
    {
        return ptr[ i ];
    }

    const T& operator[] ( const Index2& id ) const
    {
        const unsigned int idx = getLinearIndex( id.x, id.y );
        return ptr[ idx ];
    }
    T& operator[] ( const Index2& id )
    {
        const unsigned int idx = getLinearIndex( id.x, id.y );
        return ptr[ idx ];
    }

    void clear()
    {
        for( int i = 0, n = size.x * size.y; i < n; i++ )
            ptr[ i ] = T();
    }

private:
    Index2 size;
    T* ptr;
};

#endif