
#ifndef SimulationParameters_hpp
#define SimulationParameters_hpp

#include "Structs.hpp"
#include "Boundary.hpp"

#include <string>

/** Struct that holds all simulation parameters the program will use.
 * The struct holds methods to paese a struct file and construct a 
 * SimulationParameters object.
 */
struct SimulationParameters
{
public:
	/* domain parameters */
	Real xLength;
	Real yLength;
	union{
		int xCells;
		int iMax;
	};
	union{
		int yCells;
		int jMax;
	};

	/* time parameters */
	Real tau;
	Real tEnd;
	Real deltaT;
	union{
		Real tDeltaWrite;
		Real deltaVec;
	};

	/* simulation parameters */
	Real eps;
	Real omg;
	Real alpha;
	int iterMax;

	/* forces and constants */
	Real re;
	Real gx;
	Real gy;

	/* field starting values */
	Real ui;
	Real vi;
	Real pi;

	/* extended simulation parameters */
	std::string name;
	Real KarmanAngle;
	Real KarmanObjectWidth;

	/* From the bakery we take the following semantics for 'ComplexGeometry':
	 * DrivenCavity = 0
	 * ChannelFlow = 1
	 * ChannelFlowUpperHalf = 2
	 * StepFlow = 3
	 * ObstacleChannelFlow = 4
	 */
	int useComplexGeometry;

	/* inner ranges of boundaries (including objects).
	 * note: the ranges are from the inner-part of the domain and the direction
	 * included in the piece struct tells in which direction from the cell the 
	 * boundary lies. */
	std::vector<Boundary::BoundaryPiece> boundary_conditions;

public:

	/* TODO set defaults in cpp */
	SimulationParameters();

	SimulationParameters(std::string settingsfile);

	/** Write the current state of the simulation parameters to stdout. 
	*/
	void writeToSTDOUT();

	void writeSettingsFile(std::string settingsfile);

	bool readInputfile(std::string settingsfile);

private:
	void init();
};


#endif
