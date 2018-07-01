
#include "Vec2.h"

const Vec2 Vec2::ZERO( 0.0, 0.0 );

Vec2 operator*( const float& k, const Vec2& v )
{
   return Vec2( k * v.x, k * v.y );
}

const Matrix2 Matrix2::ZERO( 0.0, 0.0, 0.0, 0.0 );

const Matrix2 Matrix2::IDENTITY( 1.0, 0.0, 0.0, 1.0 );

Matrix2 operator*( const float& k, const Matrix2& m )
{
    return Matrix2( m.v[ 0 ][ 0 ] * k, m.v[ 0 ][ 1 ] * k,
                    m.v[ 1 ][ 0 ] * k, m.v[ 1 ][ 1 ] * k );
}