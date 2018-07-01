
#include "scene.h"
#include "pcg_solver.h"

#include <iostream>
#include <fstream>
#include <string>

typedef unsigned int uint;

namespace
{
	// Returns an Index2 modifying the given index by the given deltas
	inline Index2 modIndex(Index2 i, int dx = 0, int dy = 0)
	{
		return Index2(i.x + dx, i.y + dy);
	}

	// Returns safely a value from the domain 
	inline float getSafeV(Index2 & p, Array2<float> & v)
	{
		Index2 size = v.getSize();
		if (p.x >= 0 && p.x < size.x && p.y >= 0 && p.y < size.y)
			return v[p];

		return 0.0f;
	}

	// Sets a value in the domain enforcing in-range position
	inline void setSafeV(Index2 & p, Array2<float> &v, float value)
	{
		Index2 size = v.getSize();
		if (p.x >= 0 && p.x < size.x && p.y >= 0 && p.y < size.y)
			v[p] = value;
	}

	// Sets velocity X value in the domain enforcing boundary conditions
	inline void setSafeVx(Index2 & p, Array2<float> &v, float value)
	{
		Index2 size = v.getSize();
		if (p.x > 0 && p.x < size.x - 1 && p.y >= 0 && p.y < size.y)
		{
			v[p] = value;
		}
	}

	// Sets velocity Y value in the domain enforcing boundary conditions
	inline void setSafeVy(Index2 & p, Array2<float> & v, float value)
	{
		Index2 size = v.getSize();
		if (p.x >= 0 && p.x < size.x && p.y > 0 && p.y < size.y - 1)
		{
			v[p] = value;
		}
	}

	// Averages full velocity in a given cell
	inline Vec2 average(Index2 & point, Array2<float> & vx, Array2<float> & vy)
	{
		float vx1 = getSafeV(point, vx);
		float vx2 = getSafeV(modIndex(point, 1, 0), vx);
		float vy1 = getSafeV(point, vy);
		float vy2 = getSafeV(modIndex(point, 0, 1), vy);
		
		return Vec2((vx1 + vx2) * 0.5f, (vy1 + vy2) * 0.5f);
	}

	// Averages velocity X from the point of view of velocity Y within a cell
	inline float averageVX(Index2 & point, Array2<float> & v)
	{
		float a = getSafeV(point, v);
		float b = getSafeV(modIndex(point, 1, 0), v);
		float c = getSafeV(modIndex(point, 0, -1), v);
		float d = getSafeV(modIndex(point, 1, -1), v);

		float ab = (b + a) * 0.5f;
		float cd = (d + c) * 0.5f;

		return (ab + cd) * 0.5f;
	}

	// Averages velocity Y from the point of view of velocity X within a cell
	inline float averageVY(Index2 & point, Array2<float> & v)
	{
		float a = getSafeV(point, v);
		float b = getSafeV(modIndex(point, -1, 0), v);
		float c = getSafeV(modIndex(point, 0, 1), v);
		float d = getSafeV(modIndex(point, -1, 1), v);

		float ab = (b + a) * 0.5f;
		float cd = (d + c) * 0.5f;

		return (ab + cd) * 0.5f;
	}

	// Interpolate values based on its position in the current cell
	inline float interpolate(Vec2 & index, Array2<float> & vector)
	{
		float iflt = floor(index.x);
		float jflt = floor(index.y);
		
		float alphaX = index.x - iflt;
		float alphaY = index.y - jflt;

		uint x = uint(iflt);
		uint y = uint(jflt);

		Index2 p(x, y);
		Index2 p2 = modIndex(p, 1, 0);
		Index2 p3 = modIndex(p, 0, 1);
		Index2 p4 = modIndex(p, 1, 1);

		float v = getSafeV(p, vector);
		float v2 = getSafeV(p2, vector);
		float l1lerp = v * (1.0f - alphaX) + v2 * (alphaX);
		
		float v3 = getSafeV(p3, vector);
		float v4 = getSafeV(p4, vector);
		float l2lerp = v3 * (1.0f - alphaX) + v4 * (alphaX);

		return l1lerp * (1.0f - alphaY) + l2lerp * (alphaY);
	}

	// Adds coefficent to matrix if the given positions are within the matrix size
	inline void tryAddPressureCoeff(SparseMatrix<double> & m, unsigned int i, unsigned int j, float value)
	{
		if (i >= 0 && i < m.n && j >= 0 && j < m.n)
		{
			m.set_element(i, j, value);
		}
	}
}

// advection
void Fluid2::fluidAdvection( const float dt )
{
	Index2 gridSize = grid.getSize();
	Bbox2 box = grid.getDomain();
	Vec2 dmin = box.minPosition;
	Vec2 dmax = box.maxPosition;

    // ink advection
    {
		Array2<float> newInk(ink.getSize());
		for (uint i = 0; i < gridSize.x; i++)
		{
			for (uint j = 0; j < gridSize.y; j++)
			{
				// Get current position
				Index2 point(i, j);
				Vec2 domainPoint = grid.getCellPos(point);

				// Get current velocity
				Vec2 v = average(point, velocityX, velocityY);

				// Backwards integration
				Vec2 prevDomainPoint = domainPoint - dt * v;

				// Clamp out-of-domain positions
				prevDomainPoint.x = prevDomainPoint.x < dmin.x ? dmin.x : prevDomainPoint.x > dmax.x ? dmax.x : prevDomainPoint.x;
				prevDomainPoint.y = prevDomainPoint.y < dmin.y ? dmin.y : prevDomainPoint.y > dmax.y ? dmax.y : prevDomainPoint.y;

				// Get old ink
				Vec2 oldPointIndex = grid.getCellIndex(prevDomainPoint);
				float oldInk = interpolate(oldPointIndex, ink);
				// Write new ink
				newInk[point] = oldInk;
				
			}
		}
		// Update system
		ink = newInk;
    }

    // velocity advection
	if(Scene::testcase >= Scene::SMOKE)
    {	
		Index2 faceSizeX = grid.getSizeFacesX();
		Index2 faceSizeY = grid.getSizeFacesY();

		Array2<float> vx(velocityX.getSize());
		Array2<float> vy(velocityY.getSize());

		for (uint i = 0; i < faceSizeX.x; i++)
		{
			for (uint j = 0; j < faceSizeX.y; j++)
			{
				// Get current position
				Index2 point(i, j);
				Vec2 domainPoint = grid.getFaceXPos(point);
				
				// Get current velocity
				float v = getSafeV(point, velocityX);
				float u = averageVY(point, velocityY);
				Vec2 vel(v, u);

				// Backwards interpolation
				Vec2 prevDomainPoint = domainPoint - dt * vel;

				// Clamp out-of-domain positions
				prevDomainPoint.x = prevDomainPoint.x < dmin.x ? dmin.x : prevDomainPoint.x > dmax.x ? dmax.x : prevDomainPoint.x;
				prevDomainPoint.y = prevDomainPoint.y < dmin.y ? dmin.y : prevDomainPoint.y > dmax.y ? dmax.y : prevDomainPoint.y;

				// Lerp old velocity
				Vec2 faceIndexX = grid.getFaceIndex(prevDomainPoint, 0);

				float oldVx = interpolate(faceIndexX, velocityX);
				
				// Write new velocity
				setSafeVx(point, vx, oldVx);
			}
		}
	
		for (uint i = 0; i < faceSizeY.x; i++)
		{
			for (uint j = 0; j < faceSizeY.y; j++)
			{
				// Get current position
				Index2 point(i, j);
				Vec2 domainPoint = grid.getFaceYPos(point);

				// Get current velocity
				float v = averageVX(point, velocityX);
				float u = getSafeV(point, velocityY);
				Vec2 vel(v, u);

				// Backwards integration
				Vec2 prevDomainPoint = domainPoint - dt * vel;

				// Clamp out-of-domain positions
				prevDomainPoint.x = prevDomainPoint.x < dmin.x ? dmin.x : prevDomainPoint.x > dmax.x ? dmax.x : prevDomainPoint.x;
				prevDomainPoint.y = prevDomainPoint.y < dmin.y ? dmin.y : prevDomainPoint.y > dmax.y ? dmax.y : prevDomainPoint.y;
				
				// Lerp previous velocity
				Vec2 faceIndexY = grid.getFaceIndex(prevDomainPoint, 1);

				float oldVy = interpolate(faceIndexY, velocityY);

				// Write new velocity
				setSafeVy(point, vy, oldVy);
			}
		}
		
		// Update system
		velocityX = vx;
		velocityY = vy;	
    }
}

// emission
void Fluid2::fluidEmission()
{
    if( Scene::testcase >= Scene::SMOKE)
    {
		// Hardcoded transmitter ranges, my bad

        // emit source ink
        {
			for (uint i = 45; i < 55; i++)
			{
				for (uint j = 5; j < 10; j++)
				{
					Index2 emiPoint(i, j);
					ink[emiPoint] = 0.5f;
				}
			}
        }
        // emit source velocity
        {
			for (uint i = 45; i < 55; i++)
			{
				for (uint j = 5; j < 10; j++)
				{
					Index2 emiPoint(i, j);
					velocityY[emiPoint] = 5.0f;
				}
			}
        }
    }
}

// volume forces
void Fluid2::fluidVolumeForces( const float dt )
{
    if( Scene::testcase >= Scene::SMOKE )
    {
		// Update Y velocities with gravity
		Index2 faceYSize = grid.getSizeFacesY();
		for (uint i = 0; i < faceYSize.x; i++)
		{
			for (uint j = 0; j < faceYSize.y; j++)
			{
				Index2 p(i, j);
				float vy = getSafeV(p, velocityY);
				setSafeVy(p, velocityY, vy + dt * Scene::kGravity);
			}
		}
    }
}

// viscosity
void Fluid2::fluidViscosity( const float dt )
{
    if( Scene::testcase >= Scene::SMOKE )
    {
		float r = Scene::kDensity;
		float mu = Scene::kViscosity;
		
		Vec2 cellDx = grid.getCellDx();
		float dx2 = cellDx.x * cellDx.x;
		float dy2 = cellDx.y * cellDx.y;

		Index2 faceXSize = grid.getSizeFacesX();
		Index2 faceYSize = grid.getSizeFacesY();

		// Velocities X viscosity update
		Array2<float> vx(velocityX.getSize());
		for (uint i = 0; i < faceXSize.x; i++)
		{
			for (uint j = 0; j < faceXSize.y; j++)
			{
				Index2 p(i, j);

				float vijx = getSafeV(p, velocityX);
				float vip1j = getSafeV(modIndex(p, 1, 0), velocityX);
				float vim1j = getSafeV(modIndex(p, -1, 0), velocityX);
				float vijp1 = getSafeV(modIndex(p, 0, 1), velocityX);
				float vijm1 = getSafeV(modIndex(p, 0, -1), velocityX);
				
				float newVij = vijx + (dt / r) * mu * ( ((vip1j - 2.0f * vijx + vim1j) / dx2) + ((vijp1 - 2.0f * vijx + vijm1) / dy2) );
				setSafeVx(p, vx, newVij);
			}
		}

		Array2<float> vy(velocityY.getSize());
		for (uint i = 0; i < faceYSize.x; i++)
		{
			for (uint j = 0; j < faceYSize.y; j++)
			{
				Index2 p(i, j);

				float vij = getSafeV(p, velocityY);
				float vip1j = getSafeV(modIndex(p, 1, 0), velocityY);
				float vim1j = getSafeV(modIndex(p, -1, 0), velocityY);
				float vijp1 = getSafeV(modIndex(p, 0, 1), velocityY);
				float vijm1 = getSafeV(modIndex(p, 0, -1), velocityY);

				float newVij = vij + (dt / r) * mu * ( (vip1j - 2.0f * vij + vim1j) / dx2 + (vijp1 - 2.0f * vij + vijm1) / dy2 );
				setSafeVy(p, vy, newVij);
			}
		}

		velocityX = vx;
		velocityY = vy;
    }
}

// pressure
void Fluid2::fluidPressureProjection( const float dt )
{
    if( Scene::testcase >= Scene::SMOKE )
	{
		Index2 gridSize = grid.getSize();

		SparseMatrix<double> m (gridSize.x * gridSize.y, 5);

		PCGSolver<double> solver;
		solver.set_solver_parameters(1e-2, 1000);

		std::vector<double> rhs;
		rhs.resize(gridSize.x * gridSize.y);

		std::vector<double> result(gridSize.x * gridSize.y);

		Vec2 cellDx = grid.getCellDx();
		float dx2 = cellDx.x * cellDx.x;
		float dy2 = cellDx.y * cellDx.y;

		float oneByDx2 = -1.0f / dx2;
		float oneByDy2 = -1.0f / dy2;
		float rdt = -Scene::kDensity / dt;

		Index2 faceXSize = grid.getSizeFacesX();
		Index2 faceYSize = grid.getSizeFacesY();

		// SOLVE PRESSURES
		for (uint i = 0; i < gridSize.x; i++)
		{
			for (uint j = 0; j < gridSize.y; j++)
			{
				Index2 pij(i, j);

				unsigned int row = pressure.getLinearIndex(i, j);
				
				float xCoeff = 2.0f / dx2;
				float yCoeff = 2.0f / dy2;
				if (i == 0 || i == gridSize.x - 1)
					xCoeff = 1.0f / dx2;
				if (j == 0 || j == gridSize.y - 1)
					yCoeff = 1.0f / dy2;
				
				tryAddPressureCoeff(m, row, row, xCoeff + yCoeff);

				Index2 pim1j(i - 1, j);
				if(pim1j.x < gridSize.x)
					tryAddPressureCoeff(m, row, pressure.getLinearIndex(pim1j.x, pim1j.y), oneByDx2);

				Index2 pip1j(i + 1, j);
				if(pip1j.x < gridSize.x)
					tryAddPressureCoeff(m, row, pressure.getLinearIndex(pip1j.x, pip1j.y), oneByDx2);

				Index2 pijm1(i, j - 1);
				if(pijm1.y < gridSize.y)
					tryAddPressureCoeff(m, row, pressure.getLinearIndex(pijm1.x, pijm1.y), oneByDy2);

				Index2 pijp1(i, j + 1);
				if(pijp1.y < gridSize.y)
					tryAddPressureCoeff(m, row, pressure.getLinearIndex(pijp1.x, pijp1.y), oneByDy2);

				float vxm1 = getSafeV(pij, velocityX);
				float vym1 = getSafeV(pij, velocityY);
				float vxp1 = getSafeV(pip1j, velocityX);
				float vyp1 = getSafeV(pijp1, velocityY);

				float deltaVx = (vxp1 - vxm1) / cellDx.x;
				float deltaVy = (vyp1 - vym1) / cellDx.y;

				float rhsValue = rdt * (deltaVx + deltaVy);

				rhs[row] = rhsValue;
			}
		}
		
		double residual;
		int iterations;
		solver.solve(m, rhs, result, residual, iterations);

		// UPDATE ARRAY
		for (uint i = 0; i < gridSize.x; i++)
		{
			for (uint j = 0; j < gridSize.y; j++)
			{
				uint li = pressure.getLinearIndex(i, j);
				pressure[Index2(i, j)] = result[li];
			}
		}

		//std::cout << "Solver: Residual value: " << residual << std::endl << "Solver: Iterations: " << iterations << std::endl;

		// UPDATE VELOCITIES
		
		float dtr = dt / Scene::kDensity;

		Array2<float> vxtemp(velocityX.getSize());
		Array2<float> vytemp(velocityY.getSize());
		for (uint i = 0; i < gridSize.x; i++)
		{
			for (uint j = 0; j < gridSize.y; j++)
			{
				Index2 p(i, j);

				float press = getSafeV(p, pressure);

				Index2 xp = modIndex(p, 1, 0);
				float vx = getSafeV(xp, velocityX);
				float px1 = getSafeV(xp, pressure);
				float vxn1 = vx - dtr * (px1 - press) / cellDx.x;
				setSafeV(xp, vxtemp, vxn1);

				Index2 yp = modIndex(p, 0, 1);
				float vy = getSafeV(yp, velocityY);
				float py1 = getSafeV(yp, pressure);
				float vyn1 = vy - dtr * (py1 - press) / cellDx.y;
				setSafeV(yp, vytemp, vyn1);
			}
		}

		velocityX = vxtemp;
		velocityY = vytemp;
    }
}
