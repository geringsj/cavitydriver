//! Here we define our sane, well-designed and useful types and structs!
/*!
 * @author becherml, friesfn, geringsj
 * @date 2014
 */

#ifndef Structs_hpp
#define Structs_hpp

#define NULL nullptr

#define DIMENSIONS 2

typedef unsigned int uint;

typedef double Real;

struct Index
{
	Index() : i(0), j(0), k(1) {}
	Index(int i, int j, int k = 1)
		: i(i), j(j), k(k) {} 

	int i,j,k;

	int& operator[](const uint index)
	{
		switch (index)
		{
		case 0:
			return i;
		case 1:
			return j;
		case 2:
			return k;
		default:
			return i;
		}
	}
	int operator[](const uint index) const
	{
		switch (index)
		{
		case 0:
			return i;
		case 1:
			return j;
		case 2:
			return k;
		default:
			return i;
		}
	}
};
typedef Index Dimension;

struct Point 
{
	Point() : x(0), y(0), z(0) {}
	Point(Real x, Real y, Real z = 0)
		: x(x), y(y), z(z) {} 

	Real x,y,z;

	Real& operator[](const uint index)
	{
		switch (index)
		{
		case 0:
			return x;
		case 1:
			return y;
		case 2:
			return z;
		default:
			return x;
		}
	}
	Real operator[](const uint index) const
	{
		switch (index)
		{
		case 0:
			return x;
		case 1:
			return y;
		case 2:
			return z;
		default:
			return x;
		}
	}
};
typedef Point Delta;

//! Struct that holds the simulation parameters.
#include <iostream>
struct Simparam 
{
	Real xLength;
	Real yLength;
	int iMax;
	int jMax;
	Real tEnd;
	Real deltaT;
	Real tau;
	Real deltaVec;
	int iterMax;
	Real eps;
	Real omg;
	Real alpha;
	Real re;
	Real gx;
	Real gy;
	Real ui;
	Real vi;
	Real pi;
	int xProcs;
	int yProcs;
	/**
	* @brief Write the current state of the simulation parameters to stdout.
	*/
	void writeSimParamToSTDOUT()
	{
		std::cout << "SimParam: " << std::endl <<
			"xLength=" <<
			xLength << std::endl <<
			"yLength=" <<
			yLength << std::endl <<
			"iMax=" <<
			iMax << std::endl <<
			"jMax=" <<
			jMax << std::endl <<
			"tEnd=" <<
			tEnd << std::endl <<
			"deltaT=" <<
			deltaT << std::endl <<
			"tau=" <<
			tau << std::endl <<
			"deltaVec=" <<
			deltaVec << std::endl <<
			"iterMax=" <<
			iterMax << std::endl <<
			"eps=" <<
			eps << std::endl <<
			"omg=" <<
			omg << std::endl <<
			"alpha=" <<
			alpha << std::endl <<
			"re=" <<
			re << std::endl <<
			"gx=" <<
			gx << std::endl <<
			"gy=" <<
			gy << std::endl <<
			"ui=" <<
			ui << std::endl <<
			"vi=" <<
			vi << std::endl <<
			"pi=" <<
			pi << std::endl;
		/* TODO remove comment when implementing MPI stuff */
		// <<
		// "xProcs=" <<
		// xProcs << std::endl <<
		// "yProcs=" <<
		// yProcs << std::endl;
	}
};

#endif
