
#ifndef Bakery_hpp
#define Bakery_hpp

#include "GridFunction.hpp"
#include "SimulationParameters.hpp"

#include <iostream>

/* the bakery gets input and computes stuff, which is saved in a 
 * SimulationParameters object */
class Bakery {
public:
	enum class Setting {
		DrivenCavity,
		ChannelFlow,
		ChannelFlowUpperHalf,
		StepFlow,
		ObstacleChannelFlow
	};

private:
	struct Vec2 {
		Real x;
		Real y;
	};
	Setting m_setting;

	SimulationParameters m_simparams;

	GridFunction m_field;

	Vec2 m_object_corners[4];


	void initSimParamsDrivenCavityDefault();

public:
	/* construct:
	 * - Setting-struct + needed infos (make up missing, non-important parameters)
	 * - get Setting-struct + SimParams object (manipulate SimParams to respect setting)
	 */

	/* return computed setting simparams */
	SimulationParameters getSettingSimulationParameters();

};


#endif
