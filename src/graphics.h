#pragma once

#include <string>

#include "matlab_engine.h"

class graphics
{
public:

	graphics(matlab *engine);
	~graphics();

	int  plot_patch(std::string vertices, std::string faces, std::string face_data);

private:
	matlab* m_engine;
	
};