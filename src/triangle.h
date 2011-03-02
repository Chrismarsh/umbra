#pragma once

#include <vector>

#include "vertex.h"
class triangle
{
public:
	//use xyz triples in the vector
	//store the index like matlab
	triangle(int vertex1, int vertex2, int vertex3);

	//return vth vertex
	int get_vertex(int v);
private:
	//store index like matlab
	int vertex[3];

	std::vector<int> m_centre;

};