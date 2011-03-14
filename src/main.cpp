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
#include <cmath>
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
#include "bounding_rect.h"

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

		arma::mat* xyz = engine->get_double_matrix("square_nodes_5m");

		if(!xyz)
			throw std::exception(engine->get_last_error().c_str());

		engine->evaluate("clear square_nodes_5m");


		int num_nodes = xyz->n_rows;
		arma::mat rot_domain(num_nodes,3);

		//perform the triangulation
		std::cout << "Creating triangulation..." <<std::endl;
		triangulation* tri = new triangulation(engine);
		tri->create_delaunay(xyz->unsafe_col(0),xyz->unsafe_col(1));
		std::cout <<"Finished!" <<std::endl;

		//assign the "real values" too
		for (int i = 0; i<tri->get_size();i++)
		{
			triangle& t = *(tri->get_ptr(i));

			int v1 = t.get_vertex(0)-1;
			int v2 = t.get_vertex(1)-1;
			int v3 = t.get_vertex(2)-1;
			
			ptr_point p1;
			p1.x = &rot_domain(v1,0);
			p1.y = &rot_domain(v1,1);
			p1.z = &rot_domain(v1,2);

			ptr_point p2;
			p2.x = &rot_domain(v2,0);
			p2.y = &rot_domain(v2,1);
			p2.z = &rot_domain(v2,2);

			ptr_point p3;
			p3.x = &rot_domain(v3,0);
			p3.y = &rot_domain(v3,1);
			p3.z = &rot_domain(v3,2);

			t.set_vertex_values(p1,p2,p3);
		}

//send back to matlab for testing
		std::cout << "Sending back to matlab..." <<std::endl;
		mxArray*  mxTri = mxCreateDoubleMatrix(tri->get_size(),3,mxREAL);
		double* mxTri_ptr = mxGetPr(mxTri);

		
		for (int i = 0; i<tri->get_size();i++)
		{
			arma::uvec v = tri->get_tri(i);

			mxTri_ptr[i+0*tri->get_size()] = v[0];
			mxTri_ptr[i+1*tri->get_size()] = v[1];
			mxTri_ptr[i+2*tri->get_size()] = v[2];
		}

		engine->put("tri",mxTri);
		engine->put_double_matrix("mxDomain",xyz);
		mxDestroyArray(mxTri);
		mxTri_ptr = NULL;

		//start up time
		posix_time::ptime time (gregorian::date(2011,gregorian::Mar,6), 
							posix_time::hours(6)); //start at 6am
		
		posix_time::ptime end_time (gregorian::date(2011,gregorian::Mar,6), 
			posix_time::hours(18)); //end at 6pm

		//time step
		posix_time::time_duration dt = posix_time::minutes(15);

		//UTC offset. Don't know how to use datetime's UTC converter yet....
		posix_time::time_duration UTC_offset = posix_time::hours(7);
		
		//set the format
		posix_time::time_facet* facet =new posix_time::time_facet("%Y/%m/%d %H:%M:%S");
		std::stringstream ss;
		ss.imbue(std::locale(ss.getloc(),facet));
		
		//setup the plot
		engine->evaluate("ff=figure; set(gcf,'units','normalized','outerposition',[0 0 1 1]);");
		engine->evaluate("set(ff,'Renderer','OpenGL')");

		std::string viewpoint = "basin";
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

			//check negative solar elevation 
			if(E > 0)
			{	

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
					coord(0) = (*xyz)(i,0);
					coord(1) = (*xyz)(i,1);
					coord(2) = (*xyz)(i,2);

					coord = K*coord;
					rot_domain(i,0)=coord(0);
					rot_domain(i,1)=coord(1);
					rot_domain(i,2)=coord(2);

				}

				engine->put_double_matrix("mxRot",&rot_domain);

				//build bounding rect
				bounding_rect* BBR = new bounding_rect(engine);
				BBR->make(&(rot_domain.unsafe_col(0)),&(rot_domain.unsafe_col(1)),20);
			
				//not a great way of doing this, but it's compatible with matlab's plotting
				arma::vec shadows(xyz->n_rows);
				shadows.ones();

				//for each triangle
				for(int i = 0; i< tri->get_size();i++)
				{
					arma::mat c(3,3);
					triangle* t = tri->get_ptr(i);
					c << rot_domain(t->get_vertex(0)-1,0) << rot_domain(t->get_vertex(0)-1,1) << rot_domain(t->get_vertex(0)-1,2) << arma::endr
					  << rot_domain(t->get_vertex(1)-1,0) << rot_domain(t->get_vertex(1)-1,1) << rot_domain(t->get_vertex(1)-1,2) << arma::endr
					  << rot_domain(t->get_vertex(2)-1,0) << rot_domain(t->get_vertex(2)-1,1) << rot_domain(t->get_vertex(2)-1,2) << arma::endr;

					engine->put_double_matrix("t_set",&c);
					engine->evaluate("c=tri_center( [t_set(1,1) t_set(1,2) t_set(1,3)],[t_set(2,1) t_set(2,2) t_set(2,3)],[t_set(3,1) t_set(3,2) t_set(3,3)],'circumcenter')");
					arma::vec* pos = engine->get_double_vector("c");

					tri->get_ptr(i)->center.x = (*pos)(0);
					tri->get_ptr(i)->center.y = (*pos)(1);
					tri->get_ptr(i)->center.z = (*pos)(2);

					ptr_point v1;
					v1.x = &rot_domain(t->get_vertex(0)-1,0);
					v1.y = &rot_domain(t->get_vertex(0)-1,1);
					v1.z = &rot_domain(t->get_vertex(0)-1,2);

					ptr_point v2;
					v2.x = &rot_domain(t->get_vertex(1)-1,0);
					v2.y = &rot_domain(t->get_vertex(1)-1,1);
					v2.z = &rot_domain(t->get_vertex(1)-1,2);

					ptr_point v3;
					v3.x = &rot_domain(t->get_vertex(2)-1,0);
					v3.y = &rot_domain(t->get_vertex(2)-1,1);
					v3.z = &rot_domain(t->get_vertex(2)-1,2);


					tri->get_ptr(i)->set_vertex_values(v1,v2,v3);
					tri->get_ptr(i)->update_subtri();

					for(int j=0; j < BBR->n_segments;j++)
					{
						arma::uvec v = tri->get_tri(i);
						if (  BBR->pt_in_rect(rot_domain(v(0)-1,0),rot_domain(v(0)-1,1), BBR->get_rect(j)) || //pt1
							  BBR->pt_in_rect(rot_domain(v(1)-1,0),rot_domain(v(1)-1,1), BBR->get_rect(j)) || //pt2
							  BBR->pt_in_rect(rot_domain(v(2)-1,0),rot_domain(v(2)-1,1), BBR->get_rect(j)) )  //pt3
						{
							BBR->get_rect(j)->triangles.push_back(tri->get_ptr(i));
							//shadows(i)=j;
						}
					}
				}

				std::cout << "Generating shadows" << std::endl;


				

				//for each rect
				for(int i = 0; i<BBR->n_segments; i++)
				{
					std::cout << i <<std::endl;
					int num_tri = BBR->get_rect(i)->triangles.size();

					//for each triangle in each rect
					for(int j=0; j<num_tri;j++)
					{

						//jth triangle
						triangle* tj = (BBR->get_rect(i)->triangles.at(j));

						//compare to other triangles
						for(int k=0; k<num_tri;k++)
						{
							//out current kth triangle
							triangle* tk = (BBR->get_rect(i)->triangles.at(k));

							if(j != k)// && tj->get_vertex_value(0).z > tk->get_vertex_value(0).z)
							{
								//check each vertex
								if (tk->contains(tj->get_vertex_value(0).x,tj->get_vertex_value(0).y) ||
									tk->contains(tj->get_vertex_value(1).x,tj->get_vertex_value(1).y) ||
									tk->contains(tj->get_vertex_value(2).x,tj->get_vertex_value(2).y) ||
									tk->contains(tj->center.x,tj->center.y))
								{
										tk->set_shadow(true);
										shadows(tk->get_vertex(0)-1) = 1.0;

								}
								else
								{
										shadows(tk->get_vertex(0)-1) = 0.0;
								}
							}

						}

					}

				}

				double  lol = shadows(9999);
				engine->put_double_vector("shadows",&shadows);

				//plot BBR
// 				gfx->hold_on();
// 				gfx->plot_line(BBR->bbx,BBR->bby,"'color','red'");
// 				for (int i = 0; i<BBR->n_segments;i++)
// 				{
// 					gfx->plot_line(&(BBR->get_rect(i)->coord->unsafe_col(0)),&(BBR->get_rect(i)->coord->unsafe_col(1)));
// 				}
// 				gfx->hold_off();

			
				//plot
				if(viewpoint=="basin")
				{
					//for basin view
					engine->put_double_matrix("mxRot",&rot_domain);
					if(handle == -1)
						handle = gfx->plot_patch("[mxDomain(:,1) mxDomain(:,2) mxDomain(:,3)]","tri","shadows");
					else
						handle = gfx->update_patch(handle,"[mxDomain(:,1) mxDomain(:,2) mxDomain(:,3)]","shadows");
				}
				else
				{
					//as sun
	// 				if(handle == -1)
	// 					handle = gfx->plot_patch("[mxRot(:,1) mxRot(:,2)]","tri","mxRot(:,3)");
	// 				else
	// 					handle = gfx->update_patch(handle,"[mxRot(:,1) mxRot(:,2)]","mxRot(:,3)");
	// 				engine->evaluate("axis tight");

					if(handle == -1)
						handle = gfx->plot_patch("[mxRot(:,1) mxRot(:,2)]","tri","shadows(:)");
					else
						handle = gfx->update_patch(handle,"[mxRot(:,1) mxRot(:,2)]","shadows(:)");
					engine->evaluate("axis tight");

				}
					

				//update time w/o UTC offset.
				ss.str("");
				ss << time;			
		 			
				ht = gfx->add_title(ss.str());
			
			
				std::cout << "paused" << std::endl;
				_getch();
				//Sleep(50);

			}
			time = time + dt;
		}
	

		_getch();
		engine->stop();
	}
	catch(std::exception e)
		{
			std::cout << e.what() << std::endl;
			_getch();
		}
		
	
	return 0;
}








