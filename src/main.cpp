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

		matlab_engine* engine = new matlab_engine();
		graphics* gfx = new graphics(engine);


		engine->start();
		engine->set_working_dir();

		std::cout << "Matlab engine created" << std::endl;

		//loads the data via matlab
		std::cout << "Loading data" << std::endl;
		engine->evaluate("load square_nodes_2m.csv");
		arma::mat* xyz = NULL;
		xyz = engine->get_double_matrix("square_nodes_2m");

		if(!xyz)
			throw std::exception(engine->get_last_error().c_str());

		engine->evaluate("clear square_nodes_2m");

		//perform the triangulation
		std::cout << "Creating triangulation..." <<std::endl;
		triangulation* tri = new triangulation(engine);
		tri->create_delaunay(xyz->unsafe_col(0),xyz->unsafe_col(1),xyz->unsafe_col(2));

// 		std::cout << "Finding obs triangle..." << std::endl;
// 		triangle* obs_tri = tri->find_containing_triangle(626345.8844,5646903.1124);
// 		if(!obs_tri)
// 			throw std::exception("Can't find containing triangle");



		//send 
		std::cout << "Sending triangulation to matlab..." <<std::endl;
		engine->put_double_matrix("tri",&(tri->matlab_tri_matrix()));

		std::cout << "Sending domain data to matlab..." <<std::endl;
		engine->put_double_matrix("mxDomain",xyz);
		
		std::cout << "Creating face normals..." <<std::endl;
		tri->compute_face_normals();
 
		//start up time
		posix_time::ptime time (gregorian::date(2010,gregorian::Nov,9), 
							posix_time::hours(13)+posix_time::minutes(00)); //start at 6am
		
		posix_time::ptime end_time (gregorian::date(2010,gregorian::Nov,9), 
			posix_time::hours(13)+posix_time::minutes(30)); 

		//time step
		posix_time::time_duration dt = posix_time::minutes(30);

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

			arma::vec obs_pt;

			//check negative solar elevation 
			if(E > 0)
			{	

				//eqns (6) & (7) in Montero
				double	z0 = M_PI-A;
				double q0 = M_PI/2.0 - E;

				K   << cos(z0)          << sin(z0)          << 0       <<arma::endr
					<< -cos(q0)*sin(z0) << cos(q0)*cos(z0)  << sin(q0) << arma::endr
					<< sin(q0)*sin(z0)  << -cos(z0)*sin(q0) << cos(q0) << arma::endr;

				//perform the euler rotation
				for(size_t i = 0; i < tri->size(); i++)
				{
					//-------------
					//Rotate the data
					//------------------

					//get the ith triangle coordinate

					point rot_v0,rot_v1,rot_v2;
					arma::vec coord(3);

					size_t index = tri->operator()(i).global_id[0];
					//vertex 0
					coord(0) = (*xyz)(tri->operator()(i).global_id[0]-1,0);
					coord(1) = (*xyz)(tri->operator()(i).global_id[0]-1,1);
					coord(2) = (*xyz)(tri->operator()(i).global_id[0]-1,2);

					coord = K*coord;

					rot_v0.x=coord(0);
					rot_v0.y=coord(1);
					rot_v0.z=coord(2);

					//vertex 1
					index = tri->operator()(i).global_id[1];
					coord(0) = (*xyz)(tri->operator()(i).global_id[1]-1,0);
					coord(1) = (*xyz)(tri->operator()(i).global_id[1]-1,1);
					coord(2) = (*xyz)(tri->operator()(i).global_id[1]-1,2);

					coord = K*coord;

					rot_v1.x=coord(0);
					rot_v1.y=coord(1);
					rot_v1.z=coord(2);

					//vertex 2
					index = tri->operator()(i).global_id[2];
					coord(0) = (*xyz)(tri->operator()(i).global_id[2]-1,0);
					coord(1) = (*xyz)(tri->operator()(i).global_id[2]-1,1);
					coord(2) = (*xyz)(tri->operator()(i).global_id[2]-1,2);

					coord = K*coord;

					rot_v2.x=coord(0);
					rot_v2.y=coord(1);
					rot_v2.z=coord(2);

					(*tri)(i).set_vertex_values(rot_v0, rot_v1, rot_v2);
					(*tri)(i).shadow = 0.0;
					(*tri)(i).radiation = 0.0;
					(*tri)(i).z_prime = (rot_v0.z + rot_v1.z + rot_v2.z)/3.0;

					(*tri)(i).compute_azimuth();
					(*tri)(i).compute_slope();

					//---------------
					//set the initial radiation calculation
					//---------------

					//radiation data
					//solar vector
					//xyz cartesian
					arma::vec S;
					  S << cos(E) * sin(A) << arma::endr
						<< cos(E) * cos(A) << arma::endr
						<< sin(E) << arma::endr;

					//angle b/w sun and facenormal
					double angle = acos(arma::dot(S,(*tri)(i).get_facenormal()));

					//stop it wrapping around
  					if(angle > 3.14159/2.0)
 						angle = 3.14159/2.0 ;
		
					double rad =  1370.0/1.0344 *  cos(angle) *0.75; //use a "default" transmittance

					(*tri)(i).radiation = cos(angle);
				}

				//put the rotated domain to matlab
				std::cout <<"Sending rotated matrix to matlab..." <<std::endl;

				int num_nodes = xyz->n_rows;
				arma::mat rot_domain(num_nodes,3);
				rot_domain.zeros();
				for(size_t i=0;i<tri->size();i++)
				{
					rot_domain((*tri)(i).global_id[0]-1,0) = (*tri)(i).get_vertex(0).x;
					rot_domain((*tri)(i).global_id[0]-1,1) = (*tri)(i).get_vertex(0).y;
					rot_domain((*tri)(i).global_id[0]-1,2) = (*tri)(i).get_vertex(0).z;

					rot_domain((*tri)(i).global_id[1]-1,0) = (*tri)(i).get_vertex(1).x;
					rot_domain((*tri)(i).global_id[1]-1,1) = (*tri)(i).get_vertex(1).y;
					rot_domain((*tri)(i).global_id[1]-1,2) = (*tri)(i).get_vertex(1).z;

					rot_domain((*tri)(i).global_id[2]-1,0) = (*tri)(i).get_vertex(2).x;
					rot_domain((*tri)(i).global_id[2]-1,1) = (*tri)(i).get_vertex(2).y;
					rot_domain((*tri)(i).global_id[2]-1,2) = (*tri)(i).get_vertex(2).z;
				}

				engine->put_double_matrix("mxRot",&rot_domain);

/* 				std::cout <<"Building BBR..." <<std::endl;
 				//build bounding rect
 				bounding_rect* BBR = new bounding_rect(engine);
 				BBR->make(&(rot_domain.unsafe_col(0)),&(rot_domain.unsafe_col(1)),20,50);
			
				//for each triangle				for(int i = 0; i< tri->size();i++)
				{
					//for each bounding box segment
					#pragma omp parallel for
					for(int j=0; j < BBR->n_rows;j++)
					{
						for(int k = 0; k < BBR->n_cols; k++)
						{
							triangle& t = (*tri)(i);

							if (  BBR->pt_in_rect(t.get_vertex(0).x, t.get_vertex(0).y, BBR->get_rect(j,k)) || //pt1
								  BBR->pt_in_rect(t.get_vertex(1).x, t.get_vertex(1).y, BBR->get_rect(j,k)) || //pt2
								  BBR->pt_in_rect(t.get_vertex(2).x, t.get_vertex(2).y, BBR->get_rect(j,k)) )  //pt3
							{
								//t.shadow = j*k;
								BBR->get_rect(j,k)->triangles.push_back(&t);
							}
						}

					}

				}

				
				std::cout << "Generating shadows" << std::endl;

				//for each rect
				for(int i = 0; i<BBR->n_rows; i++)
				{
					for(int ii = 0; ii<BBR->n_cols;ii++)
					{
						//sort descending
						std::sort(BBR->get_rect(i,ii)->triangles.begin(),BBR->get_rect(i,ii)->triangles.end(),
							[](triangle* a, triangle* b)->bool
						{
							return a->z_prime > b->z_prime;

// 							double a_avg = (a->get_vertex(0).z + a->get_vertex(1).z + a->get_vertex(2).z)/3.0;
// 							double b_avg = (b->get_vertex(0).z + b->get_vertex(1).z + b->get_vertex(2).z)/3.0;
// 							return a_avg > b_avg;

// 							if(a->get_vertex(0).z > b->get_vertex(0).z || 
// 								a->get_vertex(0).z > b->get_vertex(1).z || 
// 								a->get_vertex(0).z >  b->get_vertex(2).z ||
// 
// 								a->get_vertex(1).z > b->get_vertex(0).z || 
// 								a->get_vertex(1).z > b->get_vertex(1).z || 
// 								a->get_vertex(1).z >  b->get_vertex(2).z ||
// 
// 								a->get_vertex(2).z > b->get_vertex(0).z || 
// 								a->get_vertex(2).z > b->get_vertex(1).z || 
// 								a->get_vertex(2).z >  b->get_vertex(2).z )
// 							{
// 								return true;
// 							}
// 							else
// 							{
// 								return false;
// 							}
							
						});


						int num_tri = BBR->get_rect(i,ii)->triangles.size();
						//for each triangle in each rect
						for(int j=0; j<num_tri;j++)
						{
							//jth triangle
							triangle* tj = (BBR->get_rect(i,ii)->triangles.at(j));

							//compare to other triangles
							for(int k=j; k<num_tri;k++)
							{
									triangle* tk = (BBR->get_rect(i,ii)->triangles.at(k));
// 									double tj_avg = (tj->get_vertex(0).z + tj->get_vertex(1).z + tj->get_vertex(2).z)/3.0;
// 									double tk_avg = (tk->get_vertex(0).z + tk->get_vertex(1).z + tk->get_vertex(2).z)/3.0;

 								    if(	 tj->z_prime > tk->z_prime)
 									{
										//does tj is above tk, and tk is shadded by tj
										int lfactor=tk->intersects(tj);

										//only update it if are making it "more" shadowed.
										if (tk->shadow < lfactor)
											tk->shadow = lfactor;

										//this is getting multiplied not matter what, even if it is a poorer estimate.
										//TODO: must fix
										//tk->radiation *= ((4.0-lfactor)/4.0);
									}
							}
						}
					}
				}

				*/
				
				arma::vec shadows(tri->size());
				arma::vec radiation(tri->size());
				arma::vec zprime(tri->size());

				//this maybe wrong?????
				for(size_t i=0;i<tri->size();i++)
				{
					shadows(i) = (*tri)(i).shadow;
					radiation(i) = (*tri)(i).radiation;
					zprime(i) = (*tri)(i).z_prime;
			
				}

				engine->put_double_vector("shadows",&shadows);
				engine->put_double_vector("radiation",&radiation);
				engine->put_double_vector("zprime",&zprime);

			


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

				gfx->colorbar();

				//for basin view
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
				engine->evaluate("set(gcf,'color','black');set(gca,'visible','off');");
//				engine->evaluate("colormap(flipud(gray))");

				//update time w/o UTC offset.
				ss.str("");
				ss << time;			


		 			
				ht = gfx->add_title(ss.str(),14,"white");
				
				
				posix_time::time_facet* fname_time_facet =new posix_time::time_facet("%Y-%m-%d-%H-%M-%S");
				std::stringstream fname_time;
				fname_time.imbue(std::locale(fname_time.getloc(),fname_time_facet));
				fname_time << time;

				gfx->save_to_file(fname_time.str());			
// 				std::cout << "paused" << std::endl;
// 				_getch();
			}
			time = time + dt;
		}
	
		engine->stop();
	}
	catch(std::exception e)
	{
		std::cout << e.what() << std::endl;																						

	}
	std::cout << "Finished" << std::endl;
	_getch();
	
	return 0;
}








