
#ifndef _SCENE_HAS_BEEN_INCLUDED_
#define _SCENE_HAS_BEEN_INCLUDED_

#include <vector>
#include "fluid2.h"
#include "fluidVisualizer2.h"

class Scene
{
public:
    enum
    {
        TEST_ADVECTION,
        SMOKE
    };

public:
    // settings
    static int testcase;
    static bool pauseFlag;

    static unsigned int nCellsX;
    static unsigned int nCellsY;
    static float step;
    static float kDensity;
    static float kGravity;
    static float kViscosity;

public:
    Scene( void );
    ~Scene( void );

    // initialization
    void init( void );
    void init( int argc, char* argv[] );
    void initAnimation( void );
    void printSettings( void );

    //Update
    void pause( void );
    void update( void );
    void animate( void );

    //Display
    void display( void );

private:
    Fluid2* fluid;
    FluidVisualizer2* fluidDisplay;

};

#endif
