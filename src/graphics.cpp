#include "graphics.h"

//returns the handle from matlab
//vertices should be in the form [x(:) y(:) z(:)]
//faces should be in the form  tri.Triangulation
int graphics::plot_patch( std::string vertices, std::string faces, std::string face_data )
{
	std::string command = std::string("path_handle = patch('Vertices',") + 	vertices + 
		std::string(",'Faces',") + 	faces +
		std::string(",'facevertexcdata',") + face_data +
		std::string(",'facecolor','flat', 'edgecolor','black');");
	m_engine->evaluate(command.c_str());

	mxArray* handle =  m_engine->get("path_handle");
	return *mxGetPr(handle);
}

graphics::graphics( matlab *engine )
{
	m_engine = engine;
}
