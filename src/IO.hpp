//! The class implements the IO 
/*!
 * @author diehlpk
 * @date 2012
 */

#include "typedef.hpp"

#ifndef IO_HPP_
#define IO_HPP_

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
  IO (char *input, char *output);

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
  void writeVTKFile (const MultiIndexType & griddimension,
		     GridFunctionType & u, GridFunctionType & v,
		     GridFunctionType & p, const PointType & delta, int step);


  /**
	* @brief Write the current state of the simulation parameters to stdout.
	*/
  void writeSimParamToSTDOUT();

  RealType getXLength() const { return simparam.xLength; };
  RealType getYLength() const { return simparam.yLength; };
  int getIMax() const { return simparam.iMax; };
  int getJMax() const { return simparam.jMax; };
  RealType getTEnd() const { return simparam.tEnd; };
  RealType getDeltaT() const { return simparam.deltaT; };
  RealType getTau() const { return simparam.tau; };
  RealType getDeltaVec() const { return simparam.deltaVec; };
  int getIterMax() const { return simparam.iterMax; };
  RealType getEps() const { return simparam.eps; };
  RealType getOmg() const { return simparam.omg; };
  RealType getAlpha() const { return simparam.alpha; };
  RealType getRe() const { return simparam.re; };
  RealType getGx() const { return simparam.gx; };
  RealType getGy() const { return simparam.gy; };
  RealType getUi() const { return simparam.ui; };
  RealType getVi() const { return simparam.vi; };
  RealType getPi() const { return simparam.pi; };
  int getXProcs() const { return simparam.xProcs; };
  int getYProcs() const { return simparam.yProcs; };

private:

//! Struct that holds the simulation parameters.
	struct {
		RealType xLength;
		RealType yLength;
		int iMax;
		int jMax;
		RealType tEnd;
		RealType deltaT;
		RealType tau;
		RealType deltaVec;
		int iterMax;
		RealType eps;
		RealType omg;
		RealType alpha;
		RealType re;
		RealType gx;
		RealType gy;
		RealType ui;
		RealType vi;
		RealType pi;
		int xProcs;
		int yProcs;
	} simparam;

//! Path where to write the vtk files.
  char *output;

/*!
   * Methods reads the simulation parameters from the specified input file.
   * 
   * @param filename The name of the file with the simulations paremters
   */
  void readInputfile (char *filename);


  //! Method interpolates the velocity for u in the staggered grid.
  /*!
   * \param x Value of the x coordinate.
   * \param y Value of the y coordinate.
   * \param u The gridfunction u.
   * \param delta The mesh width in x direction and y direction.
   */
  RealType interpolateVelocityU (RealType x, RealType y, GridFunctionType & u,
				 const PointType & delta);

  //! Method interpolates the velocity for v in the staggered grid.
  /*!
   * \param x Value of the x coordinate.
   * \param y Value of the y coordinate.
   * \param u The gridfunction v.
   * \param delta The mesh width in x direction and y direction.
   */
  RealType interpolateVelocityV (RealType x, RealType y, GridFunctionType & v,
				 const PointType & delta);

};

#endif
