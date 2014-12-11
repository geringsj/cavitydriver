
#ifndef SimulationParameters_hpp
#define SimulationParameters_hpp

#include "Structs.hpp"

#include <iostream>

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
		Real tDeltaWriteVTK;
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

public:

	SimulationParameters(std::string settingsfile);

	/** Write the current state of the simulation parameters to stdout. 
	*/
	void writeToSTDOUT();

	bool readInputfile(std::string settingsfile);

};


#endif
