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


#pragma once

#include <vector>
#include <armadillo>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>


#include "point.h"

class triangle
{
public:
	//use xyz triples in the vector
	//store the index like matlab
	triangle(point vertex1, point vertex2, point vertex3, size_t cur_rec_depth=1);
	triangle(size_t cur_rec_depth);


	/*Setters*/

	//recomputes the subtriangles
	void update_subtri();
	void set_vertex_values( point vertex1, point vertex2, point vertex3);
	void set_facenormal(arma::vec& normal);

	/*Getters*/
	
	//returns the v-th vertex[0-2] of the triangle
//	point operator()(size_t v);
	point get_vertex(size_t vertex);
	point get_center();
	double azimuth();
	double slope();

	void compute_azimuth();
	void compute_slope();

	//get the t-th subtriangle
	triangle& sub_tri(size_t t);


	bool contains(double x, double y);
	bool contains(point xy);
	int intersects(triangle* t);
	arma::vec get_facenormal();

	//information for the physical model
	double radiation;
	double shadow;
	double z_prime;

	point center;
	point rot_center;
	//this is the index used by matlab's triangulation
	//it is [1 .. N] where N is number of triangles
	//it starts at 1 because Matlab's indexing starts at 1
	size_t global_id[3];
		

private:

	//list of the vertexes
	point m_vertex_list[3];

	//center of the triangle
	point m_center;

	//always 4 subtriangles at the moment
	triangle** m_sub_tri;

	//set the number of sub triangles
	size_t m_cur_rec_depth;

	//surface normal vector
	arma::vec m_surface_normal;

	//helper functions:

	//get the mid point of a line segment
	point* midpoint(point& p1, point& p2);
	point  calc_center(triangle* t);


	//in radians.
	double m_slope;

	//positive, clockwise, from north
	double m_azimuth;





};