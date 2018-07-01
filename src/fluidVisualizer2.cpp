
#include "fluidVisualizer2.h"
#include <vector>


FluidVisualizer2::FluidVisualizer2( const Fluid2& fluid_ )
    : fluid( fluid_ )
{
}

FluidVisualizer2::~FluidVisualizer2( void )
{
    glDeleteTextures( 3, &( texture[0] ) );
}

void FluidVisualizer2::init( void )
{
    glGenTextures( 3, &( texture[0] ) );
}

void FluidVisualizer2::drawBbox( const Bbox2& bbox )
{
    const Vec2& a = bbox.minPosition;
    const Vec2& b = bbox.maxPosition;

    glBegin( GL_LINE_LOOP );
        glColor3f( 1.0f, 1.0f,1.0f );
        glVertex2f( a.x, a.y );
        glVertex2f( a.x, b.y );
        glVertex2f( b.x, b.y );
        glVertex2f( b.x, a.y );
    glEnd();
}

void FluidVisualizer2::drawGrid( void )
{
    const Vec2 dx = fluid.getGrid().getCellDx();
    const Bbox2 gridDomain = fluid.getGrid().getDomain();

    const Vec2 cellCorner1( gridDomain.minPosition );
    const Vec2 cellCorner2( gridDomain.minPosition.x, gridDomain.maxPosition.y - dx.y );
    const Vec2 cellCorner4( gridDomain.maxPosition - dx );
    const Vec2 cellCorner3( gridDomain.maxPosition.x - dx.x, gridDomain.minPosition.y );

    drawBbox( Bbox2( gridDomain.minPosition, gridDomain.maxPosition ) );
    drawBbox( Bbox2( cellCorner1, cellCorner1 + dx ) );
    drawBbox( Bbox2( cellCorner2, cellCorner2 + dx ) );
    drawBbox( Bbox2( cellCorner3, cellCorner3 + dx ) );
    drawBbox( Bbox2( cellCorner4, cellCorner4 + dx ) );
}

void FluidVisualizer2::drawInkField( void )
{
    const Bbox2 domain = fluid.grid.getDomain();
    const Index2 size = fluid.ink.getSize();
    const unsigned int width = size.x;
    const unsigned int height = size.y;

    glBindTexture( GL_TEXTURE_2D, texture[ 0 ] );

    glTexImage2D( GL_TEXTURE_2D, 0, GL_RGB8,
        width, height, 0, GL_LUMINANCE, GL_FLOAT, ( GLvoid* ) fluid.ink.getData() );

    glTexEnvi( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );

    glEnable( GL_TEXTURE_2D );

    glBegin( GL_QUADS );
        glColor3f( 0.9f, 0.9f,0.9f );
        glTexCoord2f( 0.0f, 0.0f );
        glVertex2f( domain.minPosition.x, domain.minPosition.y );
        glTexCoord2f( 0.0f, 1.0f );
        glVertex2f( domain.minPosition.x, domain.maxPosition.y );
        glTexCoord2f( 1.0f, 1.0f );
        glVertex2f( domain.maxPosition.x, domain.maxPosition.y );
        glTexCoord2f( 1.0f, 0.0f );
        glVertex2f( domain.maxPosition.x, domain.minPosition.y );
    glEnd();

    glDisable( GL_TEXTURE_2D );
}

void FluidVisualizer2::drawVelocityField( void )
{
    const float VSCALE = 0.01f;
    const Grid2& grid = fluid.getGrid();
    const Vec2 dx = grid.getCellDx();

    {
        const Array2< float >& u = fluid.getVelocityX();
        const unsigned int numverts = u.getSize().x * u.getSize().y * 2;
        std::vector< Vec2 > verts( numverts );
        for( unsigned int i = 0; i < u.getSize().x; ++i )
        for( unsigned int j = 0; j < u.getSize().y; ++j )
        {
            const Index2 index( i, j );
            const unsigned int idx = u.getLinearIndex( i, j ) * 2;
            Vec2& head = verts[ idx ];
            Vec2& tail = verts[ idx + 1 ];
            head = grid.getFaceXPos( index );
            tail = Vec2( head.x + VSCALE * u[ index ], head.y );
        }

        glColor3f( 0.9f, 0.3f, 0.3f );
        glEnableClientState( GL_VERTEX_ARRAY );
        glVertexPointer( 2, GL_FLOAT, 0, &( verts[ 0 ] ) );
        glDrawArrays( GL_LINES, 0, numverts );
        glDisableClientState( GL_VERTEX_ARRAY );
    }
    {
        const Array2< float >& v = fluid.getVelocityY();
        const unsigned int numverts = v.getSize().x * v.getSize().y * 2;
        std::vector< Vec2 > verts( numverts );
        for( unsigned int i = 0; i < v.getSize().x; ++i )
        for( unsigned int j = 0; j < v.getSize().y; ++j )
        {
            const Index2 index( i, j );
            const unsigned int idx = v.getLinearIndex( i, j ) * 2;
            Vec2& head = verts[ idx ];
            Vec2& tail = verts[ idx + 1 ];
            head = grid.getFaceYPos( index );
            tail = Vec2( head.x, head.y + VSCALE * v[ index ] );
        }

        glColor3f( 0.3f, 0.9f, 0.3f );
        glEnableClientState( GL_VERTEX_ARRAY );
        glVertexPointer( 2, GL_FLOAT, 0, &( verts[ 0 ] ) );
        glDrawArrays( GL_LINES, 0, numverts );
        glDisableClientState( GL_VERTEX_ARRAY );
    }
}
