#include "triangle.h"

triangle::triangle( int vertex1, int vertex2, int vertex3 )
{
	vertex[0] = vertex1;
	vertex[0] = vertex2;
	vertex[0] = vertex3;
}

int triangle::get_vertex( int v )
{
	//add boundary check here
	return vertex[v];
}
