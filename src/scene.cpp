
#include "scene.h"

#include <iostream>
#include <GL/glut.h>

using namespace std;

//int Scene::testcase = Scene::TEST_ADVECTION;
int Scene::testcase = Scene::SMOKE;
bool Scene::pauseFlag       = true;
unsigned int Scene::nCellsX = 100;
unsigned int Scene::nCellsY = 100;
float Scene::step           = 0.01f;
float Scene::kDensity       = 1.0f;
float Scene::kGravity       = -1.0f;
float Scene::kViscosity     = 0.001f;

Scene::Scene( void )
    : fluid( NULL ), fluidDisplay( NULL )
{
}

Scene::~Scene( void )
{
    if( fluidDisplay ) delete fluidDisplay;
    if( fluid ) delete fluid;
}

void Scene::printSettings( void )
{
    cerr << endl << "Current Settings:" << endl;
    cerr << "\t-gridsize " << nCellsX << " " << nCellsY << endl;
    cerr << "\t-step " << step << endl;
}

void Scene::init( int argc, char* argv[] )
{
    // TODO: read program parameters
    int arg = 1;
    while( arg < argc )
    {
        if( !strcmp( argv[ arg ], "-gridsize" ) )
        {
            arg++;
        }
        else if( !strcmp( argv[ arg ], "-step" ) )
        {
            arg++;
        }
        else
        {
            cerr << endl << "Unrecognized option " << argv[ arg ] << endl;
            cerr << "Usage: practica3.exe -[option1] [settings] -[option2] [settings] ..." << endl;
            cerr << "Options:" << endl;
            cerr << "\t-test [advection|smoke]" << endl;
            cerr << "\t-gridsize [gridcells x] [gridcells y]" << endl;
            cerr << "\t-step [step size in secs]" << endl;
            break;
        }
    }

    init();
}

void Scene::init( void )
{
    printSettings();

    const Index2 gridSize( nCellsX, nCellsY );
    const Bbox2  gridDomain( -2.0f, -2.0f, 2.0f,  2.0f );
    const Grid2  grid( gridDomain, gridSize );

    fluid = new Fluid2( grid );
    fluid->init();

    fluidDisplay = new FluidVisualizer2( *fluid );
    fluidDisplay->init();

    // initialize test
    initAnimation();
}

void Scene::initAnimation( void )
{
    switch( testcase )
    {
    case TEST_ADVECTION:
        {
            Array2< float >& ink = fluid->ink;
            for( unsigned int i = 2; i < ink.getSize().x / 4; ++i )
            for( unsigned int j = 2; j < ink.getSize().y / 4; ++j )
                ink[ Index2( i, j ) ] = 1.0f;

            Array2< float >& u = fluid->velocityX;
            for( unsigned int i = 0; i < u.getSize().x; ++i )
            for( unsigned int j = 0; j < u.getSize().y; ++j )
                u[ Index2( i, j ) ] = 2.0f;
            Array2< float >& v = fluid->velocityY;
            for( unsigned int i = 0; i < v.getSize().x; ++i )
            for( unsigned int j = 0; j < v.getSize().y; ++j )
                v[ Index2( i, j ) ] = 2.0f;

            break;
        }
    case SMOKE:
        {
            break;
        }
    }
}

void Scene::pause( void )
{
    pauseFlag = !pauseFlag;
}

void Scene::update( void )
{
    if( pauseFlag ) return;

    animate();
}

void Scene::animate( void )
{
    fluid->advanceStep( step );
}

void Scene::display( void )
{
    // draw ink
    fluidDisplay->drawInkField();

    // draw grid
    fluidDisplay->drawGrid();

    // draw velocity field
    fluidDisplay->drawVelocityField();

    // TODO: draw pressure on demand
}

