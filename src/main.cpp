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
#include <sstream>

#include <conio.h>

#include <armadillo>

#include <windows.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>



#include "matlab_engine.h"
#include "graphics.h"
#include "triangulation.h"

using namespace boost;
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

		int num_nodes = mxGetM(square_nodes_5m);
		//create matrix to hold our data
		arma::mat xyz(num_nodes,3);

		//store the rotated matrix data
		arma::mat rot_domain(num_nodes,3);

		//get the data we loaded
		double* ptr = mxGetPr(square_nodes_5m);
		for (int i =0; i<num_nodes; i++)
		{
			//col major in ptr!!!
			xyz(i,0) = ptr[i+0*num_nodes];
			xyz(i,1) = ptr[i+1*num_nodes];
			xyz(i,2) = ptr[i+2*num_nodes];
		}
		

		//perform the triangulation
		std::cout << "Creating triangulation..." <<std::endl;
		triangulation* tri = new triangulation(engine);
		tri->create_delaunay(xyz.unsafe_col(0),xyz.unsafe_col(1));
		std::cout <<"Finished!" <<std::endl;

		std::cout << "Sending back to matlab..." <<std::endl;
		mxArray*  mxTri = mxCreateDoubleMatrix(tri->get_size(),3,mxREAL);
		double* mxTri_ptr = mxGetPr(mxTri);

		//send back to matlab for testing

		//create the triangulation array
		for (int i = 0; i<tri->get_size();i++)
		{
			arma::uvec v = tri->get_tri(i);

			mxTri_ptr[i+0*tri->get_size()] = v[0];
			mxTri_ptr[i+1*tri->get_size()] = v[1];
			mxTri_ptr[i+2*tri->get_size()] = v[2];
		}

		engine->put("tri",mxTri);
		engine->copy_doublematrix_to("mxDomain",xyz);

		//start up time
		posix_time::ptime time (gregorian::date(2011,gregorian::Mar,6), 
							posix_time::hours(6)); //start at 6am
		
		posix_time::ptime end_time (gregorian::date(2011,gregorian::Mar,6), 
			posix_time::hours(18)); //end at 6pm

		//time step
		posix_time::time_duration dt = posix_time::minutes(15);

		//UTC offset. Don't know how to use datetime's UTC converter yet....
		posix_time::time_duration UTC_offset = posix_time::hours(7);
		
		posix_time::time_facet* facet =new posix_time::time_facet("%Y/%m/%d %H:%M:%S");
		std::stringstream ss;

		ss.imbue(std::locale(ss.getloc(),facet));
		
		//setup the plot
		engine->evaluate("ff=figure; set(gcf,'units','normalized','outerposition',[0 0 1 1]);");
		engine->evaluate("set(ff,'Renderer','OpenGL')");

		std::string viewpoint = "sun";
		//for as basin
		if(viewpoint == "basin")
			engine->evaluate(" campos(  1.0e+006 .*[ 0.6651    5.6380    0.0080] )");
		

		//plot handle
		double handle=-1.0;
		//title handle
		double ht = -1.0;

		while (time != end_time)
		{
			ss.str("");
			ss << time;
			std::cout << "Time: " << ss.str() << std::endl;

			posix_time::ptime temp_time = time + UTC_offset;
			ss.str("");
			ss << temp_time;

			//get the solar position
			//need to add UTC 7 to this
			engine->evaluate(std::string("[Az El] = SolarAzEl('") +
										ss.str() + std::string("'") +
										std::string(",50.960873, -115.187890,0);"));

			//holy LOLs. This is dangerous.
			double Az = *(mxGetPr(engine->get("Az")));
			double El = *(mxGetPr(engine->get("El")));

			//euler rotation matrix K
			arma::mat K; 

			//convert to rads
			double	A  = Az * M_PI/180;
			double E  = El * M_PI/180;

			if(E < 0)
				goto skip;

			//eqns (6) & (7) in Montero
			double	z0 = M_PI-A;
			double q0 = M_PI/2 - E;

			K   << cos(z0)          << sin(z0)          << 0       <<arma::endr
				<< -cos(q0)*sin(z0) << cos(q0)*cos(z0)  << sin(q0) << arma::endr
				<< sin(q0)*sin(z0)  << -cos(z0)*sin(q0) << cos(q0);

			//perform the euler rotation
			for(int i = 0; i<num_nodes;i++)
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

			engine->copy_doublematrix_to("mxRot",rot_domain);


			if(viewpoint=="basin")
			{
				//for basin view
				engine->copy_doublematrix_to("mxRot",rot_domain);
				if(handle == -1)
					handle = gfx->plot_patch("[mxDomain(:,1) mxDomain(:,2) mxDomain(:,3)]","tri","mxRot(:,3)");
				else
					handle = gfx->update_patch(handle,"[mxDomain(:,1) mxDomain(:,2) mxDomain(:,3)]","mxRot(:,3)");
			}
			else
			{
				//as sun
				if(handle == -1)
					handle = gfx->plot_patch("[mxRot(:,1) mxRot(:,2)]","tri","mxRot(:,3)");
				else
					handle = gfx->update_patch(handle,"[mxRot(:,1) mxRot(:,2)]","mxRot(:,3)");
				engine->evaluate("axis tight");

			}
					

			//update time w/o UTC offset.
			ss.str("");
			ss << time;			
		 			
			ht = gfx->add_title(ss.str());
			
			

			//_getch();
			Sleep(50);

skip:
			time = time + dt;
		}
	

		_getch();
		engine->stop();
	}
	catch(std::exception e)
	{
		std::cout << "Failure with error: " << e.what() << std::endl;
	}
	
	
	return 0;
}








