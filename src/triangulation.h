#pragma once

#include <vector>

//matlab engine
#include <engine.h>


#include "matlab_engine.h"
#include "triangle.h"
class triangulation
{

public:
	triangulation(matlab* engine);
	~triangulation();
	void create_delaunay(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z );
	int get_size();
	std::vector<int> get_tri(int t);
private:
	//triangulation, stored as a list of indexes into the x,y,z data. This is like Matlab
//	std::vector<std::vector<double> > m_tri;  //size * 3


	int m_size;//number of triangulations

	std::vector<triangle*> m_triangles;

	//ptr to the matlab engine
	matlab* m_engine;

};




