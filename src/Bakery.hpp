
#ifndef Bakery_hpp
#define Bakery_hpp

#include "SimulationParameters.hpp"

namespace Bakery {

	enum class Setting {
		DrivenCavity = 0,
		ChannelFlow = 1,
		ChannelFlowUpperHalf = 2,
		StepFlow = 3,
		ObstacleChannelFlow = 4
	};

	SimulationParameters get(Setting setting, SimulationParameters simpams);
};


#endif
