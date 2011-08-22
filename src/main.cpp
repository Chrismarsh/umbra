// 	Copyright (C) 2011  Chris Marsh
// 
// 	This program is free software: you can redistribute it and/or modify
// 	it under the terms of the GNU General Public License as published by
// 	the Free Software Foundation, either version 3 of the License, or
// 	(at your option) any later version.
// 
// 	This program is distributed in the hope that it will be useful,
// 	but WITHOUT ANY WARRANTY; without even the implied warranty of
// 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// 	GNU General Public License for more details.
// 
// 	You should have received a copy of the GNU General Public License
// 	along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <string.h>
#include <vector>
#include <sstream>
#include <cmath>
#include <armadillo>
#define _USE_MATH_DEFINES
#include <math.h>
#undef min
#undef max
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>
#include <boost/smart_ptr.hpp>


#include <libmaw.h>

#include "triangulation.h"
#include "bounding_rect.h"

using namespace boost;

int main()

{
	try{

		maw::matlab_engine* engine = new maw::matlab_engine();
		maw::graphics* gfx = new maw::graphics(engine);


		engine->start();
		engine->set_working_dir();

		std::cout << "Matlab engine created" << std::endl;

		//loads the data via matlab
		std::cout << "Loading data" << std::endl;
		engine->evaluate("load tin_2mtol_1mdem_nodes.csv");
		maw::d_mat xyz = engine->get_double_matrix("tin_2mtol_1mdem_nodes");

		//load skyview data
		std::cout << "Loading skyview data" << std::endl;
		engine->evaluate("load tin_2mtol_1mdem_skyview.csv");
		
		maw::d_mat skyview = engine->get_double_matrix("tin_2mtol_1mdem_skyview");

		if(!xyz)
			throw std::runtime_error(engine->get_last_error().c_str());

		if(!skyview)
			throw std::runtime_error(engine->get_last_error().c_str());

		engine->evaluate("clear tin_2mtol_1mdem_nodes");
		engine->evaluate("clear tin_2mtol_1mdem_skyview");

		//perform the triangulation
		std::cout << "Creating triangulation..." <<std::endl;
		triangulation* tri = new triangulation(engine);
		tri->create_delaunay(&(xyz->unsafe_col(0)),&(xyz->unsafe_col(1)),&(xyz->unsafe_col(2)));

 		std::cout << "Finding obs triangle..." << std::endl;
 		triangle* obs_tri = tri->find_containing_triangle(626345.8844,5646903.1124);
		if(!obs_tri)
			throw std::runtime_error("Can't find containing triangle");

		std::vector<double> obs_shortwave_remoteshadow;
		std::vector<double> obs_shortwave_selfshadow;

		//send 
		std::cout << "Sending triangulation to matlab..." <<std::endl;
		engine->put_double_matrix("tri",tri->matlab_tri_matrix());

		std::cout << "Sending domain data to matlab..." <<std::endl;
		engine->put_double_matrix("mxDomain",xyz);
		
		std::cout << "Creating 3D bounding box..." << std::endl;
		engine->evaluate("[~,cornerpoints,~,~,~] = minboundbox(mxDomain(:,1),mxDomain(:,2),mxDomain(:,3))");
		maw::d_mat cornerpoints = engine->get_double_matrix("cornerpoints");
		engine->evaluate("clear cornerpoints");

		std::cout << "Creating face normals..." <<std::endl;
		tri->compute_face_normals();
 
		std::cout << "Loading radiation data..." << std::endl;


 		engine->evaluate("load feb_1_data.csv");
	//	engine->evaluate("load season2011met.csv");
	//	engine->evaluate("load aprilmayjune.csv");
		
    	maw::d_mat radiation_data = engine->get_double_matrix("feb_1_data");
	//	arma::mat* radiation_data = engine->get_double_matrix("season2011met");
	//	maw::d_mat radiation_data = engine->get_double_matrix("aprilmayjune");
	//	engine->evaluate("clear feb_1_data");
 		int data_counter = 0;

		

		//start up time
//  		posix_time::ptime time (gregorian::date(2010,gregorian::Oct,17), 
//  							posix_time::hours(19)+posix_time::minutes(45)); //start at 6am
// 		posix_time::ptime end_time (gregorian::date(2011,gregorian::Jun,14), 
// 			posix_time::hours(12)+posix_time::minutes(15)); 

		posix_time::ptime time (gregorian::date(2011,gregorian::Feb,1), 
			posix_time::hours(7)+posix_time::minutes(00)); 
		posix_time::ptime end_time (gregorian::date(2011,gregorian::Feb,1), 
			posix_time::hours(18)+posix_time::minutes(45)); 

// 		posix_time::ptime time (gregorian::date(2011,gregorian::Apr,1), 
// 			posix_time::hours(0)+posix_time::minutes(0)); //start at 6am
// 		posix_time::ptime end_time (gregorian::date(2011,gregorian::Jun,14), 
// 			posix_time::hours(12)+posix_time::minutes(15)); 
		
		
		posix_time::ptime chkpt = time; //start with our first time 

		//time step
		posix_time::time_duration dt = posix_time::minutes(15);

		//UTC offset. Don't know how to use datetime's UTC converter yet....
		posix_time::time_duration UTC_offset = posix_time::hours(6);
		
		//set the format
		posix_time::time_facet* facet = new posix_time::time_facet("%Y/%m/%d %H:%M:%S");
		std::stringstream ss;
		ss.imbue(std::locale(ss.getloc(),facet));
		
		std::string viewpoint = "basin";
		bool plot = false;
		//plot handle
		double handle=-1.0;
		//title handle
		double ht = -1.0;

		if(plot)
		{
			//setup the plot
 			engine->evaluate("ff=figure; set(gcf,'units','normalized','outerposition',[0 0 1 1]);");
  			engine->evaluate("set(ff,'Renderer','OpenGL')");

			//for as basin
 			if(viewpoint == "basin")
  				engine->evaluate(" campos(  1.0e+006 .*[ 0.6651    5.6380    0.0080] )");
		}



		maw::d_vec cummulative_error(new arma::vec(tri->size()));
		cummulative_error->zeros();


		while (time <= end_time)
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
			double	A  = Az * M_PI/180.0;
			double E  = El * M_PI/180.0;

			
			//check negative solar elevation 
			if(El > 5)
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

					coord(0) = (*tri)(i).center.x;
					coord(1) = (*tri)(i).center.y;
					coord(2) = (*tri)(i).center.z;
					coord=K*coord;
					(*tri)(i).rot_center.x = coord(0);
					(*tri)(i).rot_center.y = coord(1);
					(*tri)(i).rot_center.z = coord(2);

					
					(*tri)(i).set_vertex_values(rot_v0, rot_v1, rot_v2);
					(*tri)(i).shadow = 0.0;
					(*tri)(i).z_prime = (*tri)(i).rot_center.z; //(rot_v0.z + rot_v1.z + rot_v2.z)/3.0;

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


					
 					double diff_obs = (*radiation_data)(data_counter,1);//col 1 = diffuse

					int i1 = ((*tri)(i).global_id[0])-1;
					int i2 = ((*tri)(i).global_id[1])-1;
					int i3 = ((*tri)(i).global_id[2])-1;
					double avg_skyview = ( (*skyview)(i1) + (*skyview)(i2) + (*skyview)(i3) ) /3.0;
					(*tri)(i).radiation_diff = diff_obs * avg_skyview; //correct for skyview factor
					
					arma::vec flat_normal(3);//xyz
					flat_normal(0)=0;
					flat_normal(1)=0;
					flat_normal(2)=1; 

					//correct the observation for the flat plane
					//double corrected_obs_direct = (*radiation_data)(data_counter,2) * cos(acos(arma::dot(S,flat_normal)))  / cos(M_PI/2.0 - E);
					//double corrected_obs_direct = (*radiation_data)(data_counter,2) / cos(M_PI/2.0 - E);
				
					//angle b/w sun and facenormal
 					double angle = acos(arma::dot(S,(*tri)(i).get_facenormal()));
 					angle = cos(angle);

					if(angle < 0.0)
						angle = 0.0;

					double obs_dir = (*radiation_data)(data_counter,2);
					double corrected_obs_direct = obs_dir  /  sin(E);
					(*tri)(i).radiation_dir  =  corrected_obs_direct * angle; //correct for the triangle aspect & slope
																			
					(*tri)(i).cosi = angle;
				}

				//put the rotated domain to matlab
				std::cout <<"Sending rotated matrix to matlab..." <<std::endl;

				int num_nodes = xyz->n_rows;
				maw::d_mat rot_domain(new arma::mat(num_nodes,3));
				rot_domain->zeros();

				for(size_t i=0;i<tri->size();i++)
				{
					(*rot_domain)((*tri)(i).global_id[0]-1,0) = (*tri)(i).get_vertex(0).x;
					(*rot_domain)((*tri)(i).global_id[0]-1,1) = (*tri)(i).get_vertex(0).y;
					(*rot_domain)((*tri)(i).global_id[0]-1,2) = (*tri)(i).get_vertex(0).z;

					(*rot_domain)((*tri)(i).global_id[1]-1,0) = (*tri)(i).get_vertex(1).x;
					(*rot_domain)((*tri)(i).global_id[1]-1,1) = (*tri)(i).get_vertex(1).y;
					(*rot_domain)((*tri)(i).global_id[1]-1,2) = (*tri)(i).get_vertex(1).z;

					(*rot_domain)((*tri)(i).global_id[2]-1,0) = (*tri)(i).get_vertex(2).x;
					(*rot_domain)((*tri)(i).global_id[2]-1,1) = (*tri)(i).get_vertex(2).y;
					(*rot_domain)((*tri)(i).global_id[2]-1,2) = (*tri)(i).get_vertex(2).z;
				}

				engine->put_double_matrix("mxRot",rot_domain);

 				std::cout <<"Building BBR..." <<std::endl;
 				//build bounding rect
 				bounding_rect BBR(engine);// = new bounding_rect(engine);
/*				maw::d_mat rot_cornerpoints(new arma::mat(cornerpoints->n_rows,2));*/
				maw::d_vec rot_cp_x(new arma::vec(cornerpoints->n_rows));
				maw::d_vec rot_cp_y(new arma::vec(cornerpoints->n_rows));


				for(int i =0; i<cornerpoints->n_rows; i++)
				{
					arma::vec coord(3);
					coord(0) = (*cornerpoints)(i,0);
					coord(1) = (*cornerpoints)(i,1);
					coord(2) = (*cornerpoints)(i,2);
					
					coord = K * coord;

					(*rot_cp_x)(i) = coord(0);
					(*rot_cp_y)(i) = coord(1);
// 					(*rot_cornerpoints)(i,0) = coord(0);
// 					(*rot_cornerpoints)(i,1) = coord(1);
				}

 				//BBR.make(&(rot_domain.unsafe_col(0)),&(rot_domain.unsafe_col(1)),50,50);
// 				maw::d_vec cp_x = arma::
// 					(&(rot_cornerpoints->col(0)));
// 				maw::d_vec cp_y(&(rot_cornerpoints->col(1)));

				BBR.make(rot_cp_x,rot_cp_y,50,50);
			
				std::cout << "Generating shadows" << std::endl;

				//for each triangle				
				for(int i = 0; i< tri->size();i++)
				{
					//for each bounding box segment
					#pragma omp parallel for
					for(int j=0; j < BBR.n_rows;j++)
					{
						for(int k = 0; k < BBR.n_cols; k++)
						{
							triangle& t = (*tri)(i);

							if (  BBR.pt_in_rect(t.get_vertex(0).x, t.get_vertex(0).y, BBR.get_rect(j,k)) || //pt1
								  BBR.pt_in_rect(t.get_vertex(1).x, t.get_vertex(1).y, BBR.get_rect(j,k)) || //pt2
								  BBR.pt_in_rect(t.get_vertex(2).x, t.get_vertex(2).y, BBR.get_rect(j,k)) )  //pt3
							{
								//t.shadow = j*k;
								BBR.get_rect(j,k)->triangles.push_back(&t);
							}
						}

					}

				}

				


				//for each rect
				#pragma omp parallel for
				for(int i = 0; i<BBR.n_rows; i++)
				{
					for(int ii = 0; ii<BBR.n_cols;ii++)
					{
						//sort descending
						std::sort(BBR.get_rect(i,ii)->triangles.begin(),BBR.get_rect(i,ii)->triangles.end(),
							[](triangle* a, triangle* b)->bool
						{
							return a->z_prime > b->z_prime;					
						});


						int num_tri = BBR.get_rect(i,ii)->triangles.size();
						//for each triangle in each rect
						for(int j=0; j<num_tri;j++)
						{
							//jth triangle
							triangle* tj = (BBR.get_rect(i,ii)->triangles.at(j));

							//compare to other triangles
							for(int k=j+1; k<num_tri;k++)
							{
									triangle* tk = (BBR.get_rect(i,ii)->triangles.at(k));
 
 								    if(	 tj->z_prime > tk->z_prime)
 									{
 										//tj is above tk, and tk is shadded by tj?
										int lfactor=tk->intersects(tj);

										//only update it if are making it "more" shadowed.
										if (tk->shadow < lfactor)
											tk->shadow = lfactor;

										//this is getting multiplied not matter what, even if it is a poorer estimate.
										//TODO: must fix
										//tk->radiation = ((4.0-tk->shadow)/4.0);
									}
							}
						}
					}
				}

			

 				maw::d_vec shadows(new arma::vec(tri->size()));
 				maw::d_vec radiation_cast(new arma::vec(tri->size()));
 				maw::d_vec radiation_self(new arma::vec(tri->size()));
				maw::d_vec cosi(new arma::vec(tri->size()));
			
				#pragma omp parallel for
				for(int i=0;i<tri->size();i++)
				{
					
					double radiation = (*tri)(i).radiation_dir * (1.0-(*tri)(i).shadow) + (*tri)(i).radiation_diff;
					double rad_self  = (*tri)(i).radiation_dir + (*tri)(i).radiation_diff;

 					(*shadows)(i) = (*tri)(i).shadow;
 					(*radiation_cast)(i) = radiation;
 					(*radiation_self)(i) = rad_self;
					(*cosi)(i) = (*tri)(i).cosi;
					(*cummulative_error)(i) += (rad_self - radiation)*900.0;
				}

 				engine->put_double_vector("shadows",shadows);
 				engine->put_double_vector("radiation",radiation_cast);
 				engine->put_double_vector("self",radiation_self);
 				engine->put_double_vector("cummError",cummulative_error);
				engine->put_double_vector("cosi",cosi);
 				engine->evaluate("shadows=1-shadows");

				
			//save the obs_triangle values
			//need to expand for more than 1 obs pt
				double radiation = (*radiation_cast)(obs_tri->triangle_id);
				double rad_self  = (*radiation_self)(obs_tri->triangle_id);

				std::cout << "Rad: " << radiation << "\t Self: " << rad_self << std::endl; 
				obs_shortwave_remoteshadow.push_back(radiation);
				obs_shortwave_selfshadow.push_back(rad_self);


// 				//plot BBR
// 				gfx->hold_on();
// 				gfx->plot_line(BBR.bbx,BBR.bby,"'color','red'");
// 				for (int i = 0; i<BBR.n_rows;i++)
// 				{
// 					for(int j = 0; j < BBR.n_cols; j++)
// 					{
// 						gfx->plot_line(&(BBR.get_rect(i,j)->coord->unsafe_col(0)),&(BBR.get_rect(i,j)->coord->unsafe_col(1)));
// 					}
// 					
// 				}
// 				gfx->hold_off();

				
				if(plot)
				{

				
					//for basin view
					if(viewpoint=="basin")
					{
					
						if(handle == -1)
							handle = gfx->plot_patch("[mxDomain(:,1) mxDomain(:,2) mxDomain(:,3)]","tri","shadows"); //self-radiation cummError
						else
							handle = gfx->update_patch(handle,"[mxDomain(:,1) mxDomain(:,2) mxDomain(:,3)]","shadows");
					}
					else
					{
						//as sun

						if(handle == -1)
							handle = gfx->plot_patch("[mxRot(:,1) mxRot(:,2)]","tri","shadows");
						else
							handle = gfx->update_patch(handle,"[mxRot(:,1) mxRot(:,2)]","shadows");
						engine->evaluate("axis tight");

					}

				//engine->evaluate("hold on;plot3(626345.8844,5646903.1124,2234.66666,'o','MarkerFaceColor','white','MarkerSize',10)");
// 				engine->evaluate("set(gcf,'color','black');set(gca,'visible','off');");
// 
// 				//engine->evaluate("colormap(flipud(jet))");
// 				gfx->colorbar();
// 
// 				//update time w/o UTC offset.
//  				ss.str("");
//  				ss << time;			
// // 		 			
//  				ht = gfx->add_title(ss.str(),14,"white");
// 				
				
// 				posix_time::time_facet* fname_time_facet = new posix_time::time_facet("%Y-%m-%d-%H-%M-%S");
// 				std::stringstream fname_time;
// 				fname_time.imbue(std::locale(fname_time.getloc(),fname_time_facet));
// 				fname_time << time;


// 
// 				gfx->colorbar("off");
// 
// 				
/*				gfx->save_to_file(fname_time.str());		*/

				}
  				
				posix_time::time_facet* fname_time_facet = new posix_time::time_facet("%Y-%m-%d-%H-%M-%S");
				std::stringstream fname_time;
				fname_time.imbue(std::locale(fname_time.getloc(),fname_time_facet));
				fname_time << time;

				engine->evaluate(std::string("save ") + fname_time.str());
			}
			else
			{
				//sun below horizon?
				obs_shortwave_remoteshadow.push_back(0);
				obs_shortwave_selfshadow.push_back(0);
			}


//  			std::cout << "Paused..." << std::endl;
//  			std::cin.get();
			time = time + dt;
			data_counter++;


		}
		arma::vec remoteshadow;
		arma::vec selfshadow;
		remoteshadow.reshape(obs_shortwave_remoteshadow.size(),1); //this is stupid and hacky. fix it later
		selfshadow.reshape(obs_shortwave_selfshadow.size(),1); //this is stupid and hacky. fix it later
		for(int i = 0; i< obs_shortwave_remoteshadow.size();i++)
		{
			remoteshadow(i) = obs_shortwave_remoteshadow[i];
			selfshadow(i) = obs_shortwave_selfshadow[i];
		}

		remoteshadow.save("Feb1_remoteshadow_values.csv", arma::raw_ascii);
		selfshadow.save("Feb1_selfshadow_values.csv", arma::raw_ascii);

		*cummulative_error = *cummulative_error / 1000000.0;

		engine->put_double_vector("error",cummulative_error);
		engine->evaluate("save Feb1_2m_error_mat");
		engine->stop();
	}
	catch(std::runtime_error e)
	{
		std::cout << e.what() << std::endl;																						

	}
 	std::cout << "Finished" << std::endl;

	return 0;
}








