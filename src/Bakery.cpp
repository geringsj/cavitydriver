
#include "Bakery.hpp"
#include "GridFunction.hpp"
#include "Boundary.hpp"

#define _USE_MATH_DEFINES
#include <math.h>

#include <list>

namespace Bakery {

	namespace {
		struct vec2 {
			Real x;
			Real y;

			vec2() : x(0), y(0) {}
			vec2(Real x, Real y) : x(x), y(y) {}
			vec2 normal() const { return vec2(-y,x); }
		};

		void setOuterBoundaries(
			SimulationParameters& simpams,
			Range inner,
			Boundary::Grid grid,
			Boundary::Condition condU,
			Real valueU,
			Boundary::Condition condD,
			Real valueD,
			Boundary::Condition condL,
			Real valueL,
			Boundary::Condition condR,
			Real valueR)
		{
		simpams.boundary_conditions.push_back(
			Boundary::BoundaryPiece(
				Boundary::Direction::Up, condU,
				grid, valueU, Range(
				Index(inner.begin.i, inner.end.j),
				Index(inner.end.i, inner.end.j) ) ));
		simpams.boundary_conditions.push_back(
			Boundary::BoundaryPiece(
				Boundary::Direction::Down, condD,
				grid, valueD, Range(
				Index(inner.begin.i, inner.begin.j),
				Index(inner.end.i, inner.begin.j) ) ));
		simpams.boundary_conditions.push_back(
			Boundary::BoundaryPiece(
				Boundary::Direction::Left, condL,
				grid, valueL, Range(
				Index(inner.begin.i, inner.begin.j),
				Index(inner.begin.i, inner.end.j) ) ));
		simpams.boundary_conditions.push_back(
			Boundary::BoundaryPiece(
				Boundary::Direction::Right, condR,
				grid, valueR, Range(
				Index(inner.end.i, inner.begin.j),
				Index(inner.end.i, inner.end.j) ) ));
		}

		void setOuterBoundariesWithObstacle(
			SimulationParameters& simpams,
			Range inner,
			GridFunction& field,
			Boundary::Grid grid,
			Boundary::Condition condU,
			Real valueU,
			Boundary::Condition condD,
			Real valueD,
			Boundary::Condition condL,
			Real valueL,
			Boundary::Condition condR,
			Real valueR)
		{
		/* up */
		for_range(i,j,
				Range(
				Index(inner.begin.i, inner.end.j),
				Index(inner.end.i, inner.end.j) )
		)
		if(field(i,j) == 0.0)
		simpams.boundary_conditions.push_back(
			Boundary::BoundaryPiece(
				Boundary::Direction::Up, condU,
				grid, valueU, Range(Index(i,j),Index(i,j)) ));
		/* down */
		for_range(i,j,
				Range(
				Index(inner.begin.i, inner.begin.j),
				Index(inner.end.i, inner.begin.j) )
		)
		if(field(i,j) == 0.0)
		simpams.boundary_conditions.push_back(
			Boundary::BoundaryPiece(
				Boundary::Direction::Down, condD,
				grid, valueD, Range(Index(i,j),Index(i,j)) ));
		/* left */
		for_range(i,j,
				Range(
				Index(inner.begin.i, inner.begin.j),
				Index(inner.begin.i, inner.end.j) )
		)
		if(field(i,j) == 0.0)
		simpams.boundary_conditions.push_back(
			Boundary::BoundaryPiece(
				Boundary::Direction::Left, condL,
				grid, valueL, Range(Index(i,j),Index(i,j)) ));
		/* right */
		for_range(i,j,
				Range(
				Index(inner.end.i, inner.begin.j),
				Index(inner.end.i, inner.end.j) ) 
		)
		if(field(i,j) == 0.0)
		simpams.boundary_conditions.push_back(
			Boundary::BoundaryPiece(
				Boundary::Direction::Right, condR,
				grid, valueR, Range(Index(i,j),Index(i,j)) ));
		}

		void getDrivenCavity(SimulationParameters& simpams, Range inner)
		{
		/* set P */
		setOuterBoundaries(simpams, inner, Boundary::Grid::P,
			Boundary::Condition::OUTFLOW, 0.0, /* UP */
			Boundary::Condition::OUTFLOW, 0.0, /* DOWN */
			Boundary::Condition::OUTFLOW, 0.0, /* LEFT */
			Boundary::Condition::OUTFLOW, 0.0); /* RIGHT */
		/* set U */
		setOuterBoundaries(simpams, inner, Boundary::Grid::U,
			Boundary::Condition::INFLOW, 1.0,
			Boundary::Condition::NOSLIP, 0.0,
			Boundary::Condition::NOSLIP, 0.0,
			Boundary::Condition::NOSLIP, 0.0);
		/* set V */
		setOuterBoundaries(simpams, inner, Boundary::Grid::V,
			Boundary::Condition::INFLOW, 0.0,
			Boundary::Condition::NOSLIP, 0.0,
			Boundary::Condition::NOSLIP, 0.0,
			Boundary::Condition::NOSLIP, 0.0);
		}

		void getChannelFlow(SimulationParameters& simpams, Range inner)
		{
		/* set P */
		setOuterBoundaries(simpams, inner, Boundary::Grid::P,
			Boundary::Condition::OUTFLOW, 0.0, /* UP */
			Boundary::Condition::OUTFLOW, 0.0, /* DOWN */
			Boundary::Condition::INFLOW, 1.0, /* LEFT */
			Boundary::Condition::OUTFLOW, 0.0); /* RIGHT */
		/* set U */
		setOuterBoundaries(simpams, inner, Boundary::Grid::U,
			Boundary::Condition::NOSLIP, 0.0,
			Boundary::Condition::NOSLIP, 0.0,
			Boundary::Condition::OUTFLOW, 0.0,
			Boundary::Condition::OUTFLOW, 0.0);
		/* set V */
		setOuterBoundaries(simpams, inner, Boundary::Grid::V,
			Boundary::Condition::NOSLIP, 0.0,
			Boundary::Condition::NOSLIP, 0.0,
			Boundary::Condition::OUTFLOW, 0.0,
			Boundary::Condition::OUTFLOW, 0.0);
		}

		void getChannelFlowUpperHalf(SimulationParameters& simpams, Range inner)
		{
		/* set P */
		setOuterBoundaries(simpams, inner, Boundary::Grid::P,
			Boundary::Condition::OUTFLOW, 0.0, /* UP */
			Boundary::Condition::OUTFLOW, 0.0, /* DOWN */
			Boundary::Condition::INFLOW, 1.0, /* LEFT */
			Boundary::Condition::OUTFLOW, 0.0); /* RIGHT */
		/* set U */
		setOuterBoundaries(simpams, inner, Boundary::Grid::U,
			Boundary::Condition::NOSLIP, 0.0,
			Boundary::Condition::SLIP, 0.0,
			Boundary::Condition::OUTFLOW, 0.0,
			Boundary::Condition::OUTFLOW, 0.0);
		/* set V */
		setOuterBoundaries(simpams, inner, Boundary::Grid::V,
			Boundary::Condition::NOSLIP, 0.0,
			Boundary::Condition::SLIP, 0.0,
			Boundary::Condition::OUTFLOW, 0.0,
			Boundary::Condition::OUTFLOW, 0.0);
		}

		void getChannelFlowObstacle(SimulationParameters& simpams, Range inner, GridFunction& field)
		{
		/* set P */
		setOuterBoundariesWithObstacle(simpams, inner, field, Boundary::Grid::P,
			Boundary::Condition::OUTFLOW, 0.0, /* UP */
			Boundary::Condition::OUTFLOW, 0.0, /* DOWN */
			Boundary::Condition::INFLOW, 1.0, /* LEFT */
			Boundary::Condition::OUTFLOW, 0.0); /* RIGHT */
		/* set U */
		setOuterBoundariesWithObstacle(simpams, inner, field, Boundary::Grid::U,
			Boundary::Condition::NOSLIP, 0.0,
			Boundary::Condition::NOSLIP, 0.0,
			Boundary::Condition::OUTFLOW, 0.0,
			Boundary::Condition::OUTFLOW, 0.0);
		/* set V */
		setOuterBoundariesWithObstacle(simpams, inner, field, Boundary::Grid::V,
			Boundary::Condition::NOSLIP, 0.0,
			Boundary::Condition::NOSLIP, 0.0,
			Boundary::Condition::OUTFLOW, 0.0,
			Boundary::Condition::OUTFLOW, 0.0);
		}

		Real dot(const vec2 a, const vec2 b)
		{
			return a.x*b.x + a.y*b.y;
		}
		vec2 operator-(const vec2& a, const vec2& b)
		{
			return vec2(a.x-b.x, a.y-b.y);
		}
		bool pointInTri(const vec2 point, const vec2* tri)
		{
			return 
				(  dot( point-tri[0], (tri[1]-tri[0]).normal() )
				 *dot( tri[2]-tri[0], (tri[1]-tri[0]).normal() ) >= 0.0)  &&

				(  dot( point-tri[1], (tri[2]-tri[1]).normal() )
				 *dot( tri[0]-tri[1], (tri[2]-tri[1]).normal() ) >= 0.0)  &&

				(  dot( point-tri[2], (tri[0]-tri[2]).normal() )
				 *dot( tri[1]-tri[2], (tri[0]-tri[2]).normal() ) >= 0.0);
		}
		bool pointInRect(vec2 point, vec2* rect)
		{
			return pointInTri(point, rect) || pointInTri(point, rect+2);
		}

		void addUVP_NOSLIP(Boundary::Direction dir, Index ind, SimulationParameters& simpams)
		{
			simpams.boundary_conditions.push_back(
					Boundary::BoundaryPiece(
						dir, Boundary::Condition::NOSLIP, Boundary::Grid::U,
						0.0, Range(ind,ind) ) );
			simpams.boundary_conditions.push_back(
					Boundary::BoundaryPiece(
						dir, Boundary::Condition::NOSLIP, Boundary::Grid::V,
						0.0, Range(ind,ind) ) );
			simpams.boundary_conditions.push_back(
					Boundary::BoundaryPiece(
						dir, Boundary::Condition::NOSLIP, Boundary::Grid::P,
						0.0, Range(ind,ind) ) );
		}
	}

	SimulationParameters get(const Setting setting, SimulationParameters simpams)
	{
		simpams.boundary_conditions.clear();

		/* update simpams according to setting */
		Range inner(Index(1,1), Index(simpams.iMax, simpams.jMax));
		switch(setting)
		{
			case Setting::DrivenCavity:
				getDrivenCavity(simpams, inner);
				break;
			case Setting::ChannelFlow:
				if(simpams.xLength < 5.0*simpams.yLength)
					simpams.xLength = 5.0 * simpams.yLength;
				getChannelFlow(simpams, inner);
				break;
			case Setting::ChannelFlowUpperHalf:
				if(simpams.xLength < 2.0*5.0*simpams.yLength)
					simpams.xLength = 2.0*5.0 * simpams.yLength;
				getChannelFlowUpperHalf(simpams, inner);
				break;
			case Setting::StepFlow:
				if(simpams.xLength < 5.0*simpams.yLength)
					simpams.xLength = 5.0 * simpams.yLength;
				/* given in assignment 4 */
				simpams.re = 100;
				/* boundaries are set below, after rasterizing step object */
				break;
			case Setting::ObstacleChannelFlow:
				if(simpams.xLength < 5.0*simpams.yLength)
					simpams.xLength = 5.0 * simpams.yLength;
				/* boundaries are set below, after rasterizing floating object */
				break;
		}

		/* set, transform and mark object */
		if(setting==Setting::StepFlow || setting==Setting::ObstacleChannelFlow)
		{
			vec2 object_corners[5];
			switch(setting)
			{
				case Setting::StepFlow:
					object_corners[0].x = 0.0; object_corners[0].y = 0.0;
					object_corners[1].x = simpams.yLength/2.; object_corners[1].y = 0.0;
					object_corners[2].x = simpams.yLength/2.; object_corners[2].y = simpams.yLength/2.;
					object_corners[3].x = 0.0; object_corners[3].y = simpams.yLength/2.;
					break;
				case Setting::ObstacleChannelFlow:
					/* set object corners using KarmanObjectWidth*/
					object_corners[0].x = 0.0; object_corners[0].y = 0.0;
					object_corners[1].x = simpams.KarmanObjectWidth; object_corners[1].y = 0.0;
					object_corners[2].x = simpams.KarmanObjectWidth; object_corners[2].y = simpams.yLength/2.;
					object_corners[3].x = 0.0; object_corners[3].y = simpams.yLength/2.;

					/* rotate object according to KarmanAngle */
					for(int i : {0,1,2,3})
					{
						/* move center to 0 */
						object_corners[i].x -= simpams.KarmanObjectWidth/2.;
						object_corners[i].y -= simpams.yLength/4.;

						/* rotate by KarmanAngle */
						Real minus90degrees = - 2.0 * M_PI / 4.0;
						Real s = sin(minus90degrees + simpams.KarmanAngle);
						Real c = cos(minus90degrees + simpams.KarmanAngle);
						Real x = object_corners[i].x;
						Real y = object_corners[i].y;
						object_corners[i].x = x*c - y*s;
						object_corners[i].y = x*s + y*c;

						/* move center of obj to final position */
						object_corners[i].x += simpams.yLength/2.;
						object_corners[i].y += simpams.yLength/2.;
					}
					break;
				default:
					break;
			}
			object_corners[4].x = object_corners[0].x;
			object_corners[4].y = object_corners[0].y;

			/* rasterize object / mark blacklisted cells under object in field */
			GridFunction field(Dimension(simpams.iMax+2,simpams.jMax+2));
			Range inner(Index(1,1),Index(simpams.iMax, simpams.jMax));
			Delta d(simpams.xLength / simpams.iMax, simpams.yLength / simpams.jMax);
			/* field has initially 0.0 everywhere */
			for_range(i,j,inner)
			if( pointInRect(vec2(d.x*i,d.y*j)-vec2(-d.x/2.,-d.y/2.), object_corners) )
			{
				field(i,j) = 1.0;
				/* fulfill boundary conditions on object:
				 * boundaries must at least be two cells thick */
				field(i+1,j) = 1.0;
				field(i,j-1) = 1.0;
				field(i+1,j-1) = 1.0;
			}
			//field.printSTDOUT();

			/* collect object boundaries in inner */
			for_range(i,j,inner)
			if(field(i,j) == 0.0)
			{
				/* left */
				if(field(i-1,j) == 1.0)
					addUVP_NOSLIP(Boundary::Direction::Left, Index(i,j), simpams);
				/* right */
				if(field(i+1,j) == 1.0)
					addUVP_NOSLIP(Boundary::Direction::Right, Index(i,j), simpams);
				/* up */
				if(field(i,j+1) == 1.0)
					addUVP_NOSLIP(Boundary::Direction::Up, Index(i,j), simpams);
				/* down */
				if(field(i,j-1) == 1.0)
					addUVP_NOSLIP(Boundary::Direction::Down, Index(i,j), simpams);
			}

			/* add missing domain boundaries depending on case (maybe do somewhere before?) */
			getChannelFlowObstacle(simpams, inner, field);
		}

		return simpams;
	}
};

