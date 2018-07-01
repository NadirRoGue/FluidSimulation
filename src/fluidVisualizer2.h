
#ifndef _FLUIDVISUALIZER2_HAS_BEEN_INCLUDED_
#define _FLUIDVISUALIZER2_HAS_BEEN_INCLUDED_

#include <iostream>
#include <GL/glut.h>

#include "fluid2.h"

class FluidVisualizer2
{
public:
    FluidVisualizer2( const Fluid2& fluid );
    ~FluidVisualizer2( void );

    void init( void );

    void drawGrid( void );
    void drawInkField( void );
    void drawVelocityField( void );
    void drawPressureField( void );

    void drawBbox( const Bbox2& bbox );

private:
    const Fluid2& fluid;

    GLuint texture[ 3 ];
};

#endif
