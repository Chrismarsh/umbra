// % Copyright (c) 2011, Chris Marsh
// 	% All rights reserved.
// 	% 
// 	% Redistribution and use in source and binary forms, with or without 
// 	% modification, are permitted provided that the following conditions are 
// 	% met:
// % 
// 	%     * Redistributions of source code must retain the above copyright 
// 	%       notice, this list of conditions and the following disclaimer.
// 	%     * Redistributions in binary form must reproduce the above copyright 
// 	%       notice, this list of conditions and the following disclaimer in 
// 	%       the documentation and/or other materials provided with the distribution
// 	%     * Neither the name of the University of Saskatchewan nor the names 
// 	%       of its contributors may be used to endorse or promote products derived 
// 	%       from this software without specific prior written permission.
// 	%       
// 	% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// 	% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// 	% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
// 	% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
// 	% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
// 	% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
// 	% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
// 	% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
// 	% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
// 	% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
// 	% POSSIBILITY OF SUCH DAMAGE.


#include <iostream>
#include <string.h>
#include <vector>

#include "matlab_engine.h"
#include "graphics.h"
#include "triangulation.h"

int main()

{
	try
	{
	
		matlab* engine = new matlab();
		graphics* gfx = new graphics(engine);
		engine->start();
		engine->set_working_dir();

		std::cout << "Matlab engine created" << std::endl;

		//loads the data via matlab
		std::cout << "Loading data" << std::endl;
		engine->evaluate("load square_nodes_5m.csv");

		mxArray* square_nodes_5m = engine->get("square_nodes_5m");

		if(!square_nodes_5m)
			throw std::exception("Variable not found in Matlab workspace");

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
		std::cout << "Creating triangulation..." <<std::endl;
		triangulation* tri = new triangulation(engine);
		tri->create_delaunay(x,y,z);
		std::cout <<"Finished!" <<std::endl;

		std::cout << "Sending back to matlab..." <<std::endl;
		mxArray* domain = mxCreateDoubleMatrix(tri->get_size(),3,mxREAL);
		double* domain_ptr = mxGetPr(domain);

		//send back to matlab for testing
		for (int i = 0; i<tri->get_size();i++)
		{
			std::vector<int> v = tri->get_tri(i);
			domain_ptr[i*3+0] = v[0];
			domain_ptr[i*3+1] = v[1];
			domain_ptr[i*3+2] = v[2];
		}

		engine->put("new_tri",domain);
		std::cout << "Finished!" << std::endl;



		int handle = gfx->plot_patch("[square_nodes_5m(:,1) square_nodes_5m(:,2) square_nodes_5m(:,3)]","new_tri","square_nodes_5m(:,3)");
		std::cout << "Plotted with handle " << handle << std::endl;

		std::string lol;
		std::cout << "Finished plotting...hopefully. Press [enter]." <<std::endl;
		std::cin >> lol;

		engine->stop();
	}
	catch(std::exception e)
	{
		std::cout << "Failure with error: " << e.what() << std::endl;
	}
	
	return 0;
}








