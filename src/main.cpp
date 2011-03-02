#include <iostream>
#include <string.h>
#include <vector>

#include "matlab_engine.h"
#include "triangulation.h"

int main()

{
	try
	{
	
		matlab* engine = new matlab();
		engine->start();

		//loads the data via matlab
		engine->evaluate("load square_nodes_5m.csv");

		mxArray* square_nodes_5m = engine->get("square_nodes_5m");
		const mwSize* num_nodes = mxGetDimensions(square_nodes_5m);

		//store our spatial data
		std::vector<double> x;
		std::vector<double> y;
		std::vector<double> z;

		double* ptr = mxGetPr(square_nodes_5m);

		//allocate data
		x.resize(num_nodes[0]);
		y.resize(num_nodes[0]);
		z.resize(num_nodes[0]);

		for (int i =0; i<num_nodes[0]; i++)
		{
			x[i] = ptr[i*3+0];
			y[i] = ptr[i*3+1];
			z[i] = ptr[i*3+2];
		}
		
		//perform the triangulation
		triangulation* tri = new triangulation(engine);
		tri->create_delaunay(x,y,z);

		

		engine->stop();
	}
	catch(std::exception e)
	{
		std::cout << "Failure with error: " << e.what() << std::endl;
	}
	
	return 0;
}








