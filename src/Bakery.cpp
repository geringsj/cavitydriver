
#include "Bakery.hpp"
#include "GridFunction.hpp"
#include "Boundary.hpp"

#include <cmath>

namespace Bakery {

	namespace {
		struct vec2 {
			Real x;
			Real y;
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
				Boundary::Direction::Up, 
				condU, 
				grid, valueU, Range(
				Index(inner.begin.i, inner.end.j),
				Index(inner.end.i, inner.end.j) ) ));
		simpams.boundary_conditions.push_back(
			Boundary::BoundaryPiece(
				Boundary::Direction::Down, 
				condD, 
				grid, valueD, Range(
				Index(inner.begin.i, inner.begin.j),
				Index(inner.end.i, inner.begin.j) ) ));
		simpams.boundary_conditions.push_back(
			Boundary::BoundaryPiece(
				Boundary::Direction::Left, 
				condL, 
				grid, valueL, Range(
				Index(inner.begin.i, inner.begin.j),
				Index(inner.begin.i, inner.end.j) ) ));
		simpams.boundary_conditions.push_back(
			Boundary::BoundaryPiece(
				Boundary::Direction::Right, 
				condR,
				grid, valueR, Range(
				Index(inner.end.i, inner.begin.j),
				Index(inner.end.i, inner.end.j) ) ));
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
				/* given in assignment 4 */
				simpams.re = 100;
				/* set boundaries ? */
				break;
			case Setting::ObstacleChannelFlow:
				/* set boundaries ?:
				 * set velocities flow profile via analytic solution to channel flow */
				break;
		}

		/* set, transform and mark object */
		if(setting==Setting::StepFlow || setting==Setting::ObstacleChannelFlow)
		{
			vec2 object_corners[4];
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
						/* move to 0 */
						object_corners[i].x -= simpams.KarmanObjectWidth/2.;
						object_corners[i].y -= simpams.yLength/2.;

						/* rotate by KarmanAngle */
						Real minus45degrees = - 2.0 * M_PI / 8.0;
						Real s = sin(minus45degrees + simpams.KarmanAngle);
						Real c = cos(minus45degrees + simpams.KarmanAngle);
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

			/* rasterize object / mark blacklisted cells under object in field */
			GridFunction field(Dimension(simpams.iMax+2,simpams.jMax+2));

			/* fulfill boundary conditions on object:
			 * boundaries must at least be two cells thick */

			/* collect object boundaries in inner */

			/* add missing domain boundaries depending on case (maybe do somewhere before?) */
		}

		return simpams;
	}
};

