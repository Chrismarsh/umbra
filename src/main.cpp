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
//#include <omp.h>
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

		std::cout.precision(30);
// 		std::cout << "Wait for debug attach" <<std::endl;
// 		_getch();
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

		//perform the triangulation
		std::cout << "Creating triangulation..." <<std::endl;
		triangulation* tri = new triangulation(engine);
		tri->create_delaunay(xyz->unsafe_col(0),xyz->unsafe_col(1),xyz->unsafe_col(2));

		std::cout << "Creating face normals..." <<std::endl;
		tri->compute_face_normals();


		//send 
		std::cout << "Sending triangulation to matlab..." <<std::endl;
		engine->put_double_matrix("tri",&(arma::conv_to<arma::mat>::from(tri->get_tri_index())));
		engine->put_double_matrix("mxDomain",xyz);


		int num_nodes = xyz->n_rows;
		arma::mat rot_domain(num_nodes,3);
 
		//start up time
		posix_time::ptime time (gregorian::date(2010,gregorian::Sep,15), 
							posix_time::hours(15)+posix_time::minutes(30)); //start at 6am
		
		posix_time::ptime end_time (gregorian::date(2010,gregorian::Sep,15), 
			posix_time::hours(15)+posix_time::minutes(45)); //end at 6pm

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


			double Az = engine->get_scaler("Az");
			double El = engine->get_scaler("El");

			//euler rotation matrix K
			arma::mat K; 

			//convert to radians
			double	A  = Az * M_PI/180;
			double E  = El * M_PI/180;
			arma::vec S;
			S << cos(El) * sin(Az) << arma::endr
			  << cos(El) * cos(Az) << arma::endr
			  << sin(El) << arma::endr;

			//check negative solar elevation 
			if(E > 0)
			{	

				//eqns (6) & (7) in Montero
				double	z0 = M_PI-A;
				double q0 = M_PI/2 - E;

				K   << cos(z0)          << sin(z0)          << 0       <<arma::endr
					<< -cos(q0)*sin(z0) << cos(q0)*cos(z0)  << sin(q0) << arma::endr
					<< sin(q0)*sin(z0)  << -cos(z0)*sin(q0) << cos(q0) << arma::endr;

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
				//need to update what the triangle data points to, as it will be by default initialized to that of the triangulation data
				tri->set_vertex_data(rot_domain);
				engine->put_double_matrix("mxRot",&rot_domain);


		
			
				//not a great way of doing this, but it's compatible with matlab's plotting
				arma::vec shadows(tri->get_num_tri());
				shadows.zeros();
				arma::vec radiation(tri->get_num_tri());
				radiation.zeros();


				std::cout <<"Building BBR..." <<std::endl;

				//build bounding rect
				bounding_rect* BBR = new bounding_rect(engine);
				BBR->make(&(rot_domain.unsafe_col(0)),&(rot_domain.unsafe_col(1)),20,20);
				
				//for each triangle
				for(size_t i = 0; i< tri->get_num_tri();i++)
				{
					//for each bounding box segment
					#pragma omp parallel for
					for(int j=0; j < BBR->n_rows;j++)
					{
						for(int k = 0; k < BBR->n_cols; k++)
						{
							triangle& t = tri->operator()(i);

							if (  BBR->pt_in_rect(t(0).x, t(0).y, BBR->get_rect(j,k)) || //pt1
								  BBR->pt_in_rect(t(1).x, t(1).y, BBR->get_rect(j,k)) || //pt2
								  BBR->pt_in_rect(t(2).x, t(2).y, BBR->get_rect(j,k)) )  //pt3
							{
								BBR->get_rect(j,k)->triangles.push_back(&tri->operator()(i));
								BBR->get_rect(j,k)->m_globalID.push_back(i);
								//shadows(i)=j*k; //uncommenting this will color the triangles based on BBR segment. good for debugging


							}
						}
					}

					//radiation data
					//solar vector
					arma::vec s_t;
					s_t << sin(S(1))*cos(S(2)) << arma::endr
						<< sin(S(1))*sin(S(2)) << arma::endr
						<< cos (S(1)) << arma::endr;
					//face normal
					arma::vec normal = tri->operator()(i).get_facenormal();
					triangle lol = tri->operator()(i);
					double angle = acos(arma::norm_dot(s_t,tri->operator()(i).get_facenormal()));
					double rad =  1370.0/1.0344 *  cos(angle) *0.75; //use a "default" transmittance
					rad = rad <0 ? 0: rad;
					radiation(i) = rad;
				}

				
				std::cout << "Generating shadows" << std::endl;

				//for each rect
				for(int i = 0; i<BBR->n_rows; i++)
				{
					for(int ii = 0; ii<BBR->n_cols;ii++)
					{
						int num_tri = BBR->get_rect(i,ii)->triangles.size();

						//for each triangle in each rect
						for(int j=0; j<num_tri;j++)
						{
							//jth triangle
							triangle* tj = (BBR->get_rect(i,ii)->triangles.at(j));

							//compare to other triangles
							#pragma omp parallel for
							for(int k=0; k<num_tri;k++)
							{
								//out current kth triangle
								triangle* tk = (BBR->get_rect(i,ii)->triangles.at(k));
								arma::uvec zj = tri->get_index(BBR->get_rect(i,ii)->m_globalID[j]);
								double z1 = rot_domain(zj(0)-1,2);

								arma::uvec zk = tri->get_index(BBR->get_rect(i,ii)->m_globalID[k]);
								double z2 = rot_domain(zk(0)-1,2);


								if(j != k 
									&& shadows(BBR->get_rect(i,ii)->m_globalID[k]) == 0.0 

									/*&& (rot_domain(zj(0)-1,2)+rot_domain(zj(1)-1,2)+rot_domain(zj(2)-1,2))/3 >
										(rot_domain(zk(0)-1,2)+ rot_domain(zk(1)-1,2)+rot_domain(zk(2)-1,2))/3)*/


									&&  (rot_domain(zj(0)-1,2) >  rot_domain(zk(0)-1,2) &&
										 rot_domain(zj(1)-1,2) >  rot_domain(zk(1)-1,2) &&
										 rot_domain(zj(2)-1,2) >  rot_domain(zk(2)-1,2)))
								{

									int lfactor=tk->intersects(tj);
									shadows(BBR->get_rect(i,ii)->m_globalID[k]) = (4-lfactor)/4.0;
									radiation(BBR->get_rect(i,ii)->m_globalID[k]) *=  1 - (4-lfactor)/4.0;
								

								}

							}

					}
					}
				}

				delete BBR;

				engine->put_double_vector("shadows",&shadows);
				engine->put_double_vector("radiation",&radiation);
				engine->put_double_matrix("mxRot",&rot_domain);
				
// 				//plot BBR
// 				gfx->hold_on();
// 				gfx->plot_line(BBR->bbx,BBR->bby,"'color','red'");
// 				for (int i = 0; i<BBR->n_rows;i++)
// 				{
// 					for(int j = 0; j < BBR->n_cols; j++)
// 					{
// 						gfx->plot_line(&(BBR->get_rect(i,j)->coord->unsafe_col(0)),&(BBR->get_rect(i,j)->coord->unsafe_col(1)));
// 					}
// 					
// 				}
// 				gfx->hold_off();

				//for basin view

				
				//plot
				if(viewpoint=="basin")
				{
					
					if(handle == -1)
						handle = gfx->plot_patch("[mxDomain(:,1) mxDomain(:,2) mxDomain(:,3)]","tri","radiation(:)");
					else
						handle = gfx->update_patch(handle,"[mxDomain(:,1) mxDomain(:,2) mxDomain(:,3)]","radiation(:)");
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
						handle = gfx->plot_patch("[mxRot(:,1) mxRot(:,2)]","tri","radiation(:)");
					else
						handle = gfx->update_patch(handle,"[mxRot(:,1) mxRot(:,2)]","radiation(:)");
					engine->evaluate("axis tight");

				}
					
//				engine->evaluate("colormap(flipud(gray))");

				//update time w/o UTC offset.
				ss.str("");
				ss << time;			
		 			
				ht = gfx->add_title(ss.str());
			
// 				std::cout << "paused" << std::endl;
// 				_getch();
			}
			time = time + dt;
		}
	
		delete xyz;

		std::cout << "Finished" << std::endl;

	//	_getch();
		engine->stop();

		delete engine;
		delete gfx;
		delete tri;

		return 0;

	}
	catch(std::exception e)
		{
			std::cout << e.what() << std::endl;
			_getch();
		}
		
	
	return 0;
}








