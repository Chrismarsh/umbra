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
#include <armadillo>

#define _USE_MATH_DEFINES
#include <math.h>

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
		int size = num_nodes[0];
		//create matrix to hold our data
		arma::mat xyz(size,3);

		//store the rotated matrix data
		arma::mat rot_domain(size,3);

		//get the data we loaded
		double* ptr = mxGetPr(square_nodes_5m);
		for (int i =0; i<size; i++)
		{
			//col major in ptr!!!
			xyz(i,0) = ptr[i+0*size];
			xyz(i,1) = ptr[i+1*size];
			xyz(i,2) = ptr[i+2*size];
		}
		

		//perform the triangulation
		std::cout << "Creating triangulation..." <<std::endl;
		triangulation* tri = new triangulation(engine);
		tri->create_delaunay(xyz.unsafe_col(0),xyz.unsafe_col(1));
		std::cout <<"Finished!" <<std::endl;

		//pretend this is the start of the time loop....

		//get the solar position
		//need to add UTC 7 to this
		engine->evaluate("[Az El] = SolarAzEl('2010/10/01 20:00:00',50.960873, -115.187890,0);");
		double Az = *(mxGetPr(engine->get("Az")));
		double El = *(mxGetPr(engine->get("El")));

		//euler rotation matrix K
		arma::mat K; 

		//convert to rads
		double	A  = Az * M_PI/180;
		double E  = El * M_PI/180;

		//eqns (6) & (7) in Montero
		double	z0 = M_PI-A;
		double q0 = M_PI/2 - E;

		K << cos(z0) << sin(z0) << 0 <<arma::endr
			<< -cos(q0)*sin(z0) << cos(q0)*cos(z0) << sin(q0) << arma::endr
			<< sin(q0)*sin(z0) << -cos(z0)*sin(q0) << cos(q0);

		//perform the euler rotation
		for(int i = 0; i<size;i++)
		{
			arma::vec coord(3);
			coord(0) = xyz(i,0);
			coord(1) = xyz(i,1);
			coord(2) = xyz(i,2);

			coord = K*coord;
			rot_domain(i,0)=coord(0);
			rot_domain(i,1)=coord(1);
			rot_domain(i,2)=coord(2);

		}
		
		std::cout << "Sending back to matlab..." <<std::endl;
		mxArray*  mxTri = mxCreateDoubleMatrix(tri->get_size(),3,mxREAL);
		double* mxTri_ptr = mxGetPr(mxTri);

		//send back to matlab for testing
		for (int i = 0; i<tri->get_size();i++)
		{
			arma::uvec v = tri->get_tri(i);

			int v1 = v[0];
			int v2 = v[1];
			int v3 = v[3];

			mxTri_ptr[i+0*tri->get_size()] = v[0];
			mxTri_ptr[i+1*tri->get_size()] = v[1];
			mxTri_ptr[i+2*tri->get_size()] = v[2];
		}

		engine->put("tri",mxTri);

		mxArray*  mxDomain = mxCreateDoubleMatrix(xyz.n_rows,3,mxREAL);
		double*   mxDomain_ptr = mxGetPr(mxDomain);
		//send back to matlab for testing
		for (int i = 0; i<xyz.n_rows;i++)
		{
			mxDomain_ptr[i+0*xyz.n_rows] = xyz(i,0);
			mxDomain_ptr[i+1*xyz.n_rows] = xyz(i,1);
			mxDomain_ptr[i+2*xyz.n_rows] = xyz(i,2);
		}
		engine->put("mxDomain",mxDomain);
		std::cout << "Finished!" << std::endl;



		int handle = gfx->plot_patch("[mxDomain(:,1) mxDomain(:,2) mxDomain(:,3)]","tri","mxDomain(:,3)");
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








