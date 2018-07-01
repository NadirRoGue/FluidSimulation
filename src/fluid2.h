
#ifndef _FLUID2_HAS_BEEN_INCLUDED_
#define _FLUID2_HAS_BEEN_INCLUDED_

#include "array2.h"
#include "grid2.h"

class Fluid2
{
public:
    Fluid2( const Grid2& grid_ )
        : grid( grid_ )
    {
    }

    void init( void )
    {
        velocityX.resize( grid.getSizeFacesX() );
        velocityY.resize( grid.getSizeFacesY() );
        ink.resize( grid.getSize() );
        pressure.resize( grid.getSize() );
    }

    const Grid2& getGrid( void ) const
    {
        return grid;
    }

    const Array2< float >& getVelocityX( void ) const
    {
        return velocityX;
    }

    const Array2< float >& getVelocityY( void ) const
    {
        return velocityY;
    }

    const Array2< float >& getPressure( void ) const
    {
        return pressure;
    }

    const Array2< float >& getInk( void ) const
    {
        return ink;
    }

    void advanceStep( const float dt )
    {
        fluidEmission();
        fluidAdvection( dt );
        fluidVolumeForces( dt );
        fluidViscosity( dt );
        fluidPressureProjection( dt );
    }

    // emission
    void fluidEmission();

    // advection
    void fluidAdvection( const float dt );

    // external forces
    void fluidVolumeForces( const float dt );

    // viscosity
    void fluidViscosity( const float dt );

    // pressure
    void fluidPressureProjection( const float dt );

public:
    Grid2 grid;

    Array2< float > velocityX;
    Array2< float > velocityY;
    Array2< float > pressure;
    Array2< float > ink;
};

#endif
