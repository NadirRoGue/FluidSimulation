
#ifndef _VEC2_HAS_BEEN_INCLUDED_
#define _VEC2_HAS_BEEN_INCLUDED_

#include <math.h>

class Vec2
{
public:
    float x, y;

public:
    Vec2( void )
        : x( 0.0f ), y( 0.0f )
    {}

    Vec2( const float x_, const float y_ )
        : x( x_ ), y( y_ )
    {}

    void operator+=( const Vec2& v )
    {
        x += v.x; y += v.y;
    }

    void operator-=( const Vec2& v )
    {
        x -= v.x; y -= v.y;
    }

    void operator*=( const Vec2& v )
    {
        x *= v.x; y *= v.y;
    }

    Vec2 operator+( const float& k ) const
    {
        return Vec2( x + k, y + k );
    }

    Vec2 operator-( const float& k ) const
    {
        return Vec2( x - k, y - k );
    }

    Vec2 operator*( const float& k ) const
    {
        return Vec2( x * k, y * k );
    }

    Vec2 operator+( const Vec2& v ) const
    {
        return Vec2( x + v.x, y + v.y );
    }

    Vec2 operator-( const Vec2& v ) const
    {
        return Vec2( x - v.x, y - v.y );
    }

    Vec2 operator*( const Vec2& v ) const
    {
        return Vec2( x * v.x, y * v.y );
    }

    Vec2 operator/( const Vec2& v ) const
    {
        return Vec2( x / v.x, y / v.y );
    }

    float dot( const Vec2& v ) const
    {
        return x * v.x + y * v.y;
    }

    float cross( const Vec2& v ) const
    {
        return x * v.y - y * v.x;
    }

    Vec2 rotate90( void ) const
    {
        return Vec2( -y, x );
    }

    float lengthSqr( void ) const
    {
        return x * x + y * y;
    }

    float length( void ) const
    {
        return sqrt( x * x + y * y );
    }

public:
    const static Vec2 ZERO;

    static Vec2 crossProd( const Vec2& v )
    {
        return Vec2( -v.y, v.x );
    }
};

Vec2 operator*( const float& k, const Vec2& v );


class Matrix2
{
public:
    float v[ 2 ][ 2 ];

    Matrix2( void )
    {
        v[ 0 ][ 0 ] = 0.0f; v[ 0 ][ 1 ] = 0.0f;
        v[ 1 ][ 0 ] = 0.0f; v[ 1 ][ 1 ] = 0.0f;
    }
    Matrix2( float a, float b, float c, float d )
    {
        v[ 0 ][ 0 ] = a; v[ 0 ][ 1 ] = b;
        v[ 1 ][ 0 ] = c; v[ 1 ][ 1 ] = d;
    }

    void operator+=( const Matrix2& m )
    {
        v[ 0 ][ 0 ] += m.v[ 0 ][ 0 ]; v[ 0 ][ 1 ] += m.v[ 0 ][ 1 ];
        v[ 1 ][ 0 ] += m.v[ 1 ][ 0 ]; v[ 1 ][ 1 ] += m.v[ 1 ][ 1 ];
    }
    void operator-=( const Matrix2& m )
    {
        v[ 0 ][ 0 ] -= m.v[ 0 ][ 0 ]; v[ 0 ][ 1 ] -= m.v[ 0 ][ 1 ];
        v[ 1 ][ 0 ] -= m.v[ 1 ][ 0 ]; v[ 1 ][ 1 ] -= m.v[ 1 ][ 1 ];
    }

    Matrix2 operator+( const Matrix2& m ) const
    {
        return Matrix2( v[ 0 ][ 0 ] + m.v[ 0 ][ 0 ], v[ 0 ][ 1 ] + m.v[ 0 ][ 1 ],
                        v[ 1 ][ 0 ] + m.v[ 1 ][ 0 ], v[ 1 ][ 1 ] + m.v[ 1 ][ 1 ] );
    }

    Matrix2 operator-( const Matrix2& m ) const
    {
        return Matrix2( v[ 0 ][ 0 ] - m.v[ 0 ][ 0 ], v[ 0 ][ 1 ] - m.v[ 0 ][ 1 ],
                        v[ 1 ][ 0 ] - m.v[ 1 ][ 0 ], v[ 1 ][ 1 ] - m.v[ 1 ][ 1 ] );
    }

    Vec2 operator*( const Vec2& m ) const
    {
        return Vec2( v[ 0 ][ 0 ] * m.x + v[ 0 ][ 1 ] * m.y,
                     v[ 1 ][ 0 ] * m.x + v[ 1 ][ 1 ] * m.y );
    }

    Matrix2 operator*( const float& k ) const
    {
        return Matrix2( v[ 0 ][ 0 ] * k, v[ 0 ][ 1 ] * k,
                        v[ 1 ][ 0 ] * k, v[ 1 ][ 1 ] * k );
    }

    Matrix2 operator*( const Matrix2& m ) const
    {
        return Matrix2( v[ 0 ][ 0 ] * m.v[ 0 ][ 0 ] + v[ 0 ][ 1 ] * m.v[ 1 ][ 0 ],
                        v[ 0 ][ 0 ] * m.v[ 0 ][ 1 ] + v[ 0 ][ 1 ] * m.v[ 1 ][ 1 ],
                        v[ 1 ][ 0 ] * m.v[ 0 ][ 0 ] + v[ 1 ][ 1 ] * m.v[ 1 ][ 0 ],
                        v[ 1 ][ 0 ] * m.v[ 0 ][ 1 ] + v[ 1 ][ 1 ] * m.v[ 1 ][ 1 ] );
    }

    const static Matrix2 ZERO;

    const static Matrix2 IDENTITY;

    friend Matrix2 operator*( float k, const Matrix2& m )
    {
        return Matrix2( k * m.v[ 0 ][ 0 ], k * m.v[ 0 ][ 1 ],
                        k * m.v[ 1 ][ 0 ], k * m.v[ 1 ][ 1 ] );
    }
};

Matrix2 operator*( const float& k, const Matrix2& m );

#endif
