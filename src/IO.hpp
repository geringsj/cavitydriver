//! The class implements the IO 
/*!
 * @author diehlpk
 * @date 2012
 */
/* we haven't really changed the IO, so we'll let diehlpk have this one... */

#ifndef IO_HPP_
#define IO_HPP_

#include "Structs.hpp"
#include "Debug.hpp"

#include "GridFunction.hpp"

#include <iostream>
#include <fstream>

#ifdef _MSC_VER
#include <windows.h>
#endif

class IO
{
public:

  /*! Construktor
   * @param input Path to the file with the simulation parameters.
   * @param outout Path to the directory for the vtk files.
   */
  IO (const char *input, const char *output);

  IO(int argc, char** argv);

  //! Destructor
   ~IO ();

  //! Method writes the GridFunctions u,v,p in the vtk data format to the hard disk. 
  //! The files are named in the following convention: field_(step).vts.
  /*!
   * \param griddimension The dimension of the gridfunctions.
   * \param u The gridfunction u.
   * \param v The gridfunction v.
   * \param p The gridfunction p.
   * \param delta The mesh width in x direction and y direction
   * \param step The number of the timestep.
   */
  void writeVTKFile (const Index & griddimension,
		     GridFunction & u, GridFunction & v,
		     GridFunction & p, const Point & delta, int step);


  /**
	* @brief Write the current state of the simulation parameters to stdout.
	*/
  void writeSimParamToSTDOUT();

  Real getXLength() const { return simparam.xLength; };
  Real getYLength() const { return simparam.yLength; };
  int getIMax() const { return simparam.iMax; };
  int getJMax() const { return simparam.jMax; };
  Real getTEnd() const { return simparam.tEnd; };
  Real getDeltaT() const { return simparam.deltaT; };
  Real getTau() const { return simparam.tau; };
  Real getDeltaVec() const { return simparam.deltaVec; };
  int getIterMax() const { return simparam.iterMax; };
  Real getEps() const { return simparam.eps; };
  Real getOmg() const { return simparam.omg; };
  Real getAlpha() const { return simparam.alpha; };
  Real getRe() const { return simparam.re; };
  Real getGx() const { return simparam.gx; };
  Real getGy() const { return simparam.gy; };
  Real getUi() const { return simparam.ui; };
  Real getVi() const { return simparam.vi; };
  Real getPi() const { return simparam.pi; };
  int getXProcs() const { return simparam.xProcs; };
  int getYProcs() const { return simparam.yProcs; };

  /*!
  * Methods reads the simulation parameters from the specified input file.
  */
  Simparam readInputfile();

private:
  Simparam simparam;

  // default values for the output directory and the settings file
  const char *output = "out";
  const char *settings = "inputvals";

  /**
   * Print the help text in order to show what
   * arguments can be processed.
   */
  void dieSynopsis();

  /**
	* Parse the input arguments. If there are no arguments
	* the default values are used.
	* @param the length of the input
	* @param the arguments
	*/
  void parseArguments(int argc, char** argv);

  //! Method interpolates the velocity for u in the staggered grid.
  /*!
   * \param x Value of the x coordinate.
   * \param y Value of the y coordinate.
   * \param u The gridfunction u.
   * \param delta The mesh width in x direction and y direction.
   */
  Real interpolateVelocityU (Real x, Real y, GridFunction & u,
				 const Point & delta);

  //! Method interpolates the velocity for v in the staggered grid.
  /*!
   * \param x Value of the x coordinate.
   * \param y Value of the y coordinate.
   * \param u The gridfunction v.
   * \param delta The mesh width in x direction and y direction.
   */
  Real interpolateVelocityV (Real x, Real y, GridFunction & v,
				 const Point & delta);

};

#endif
