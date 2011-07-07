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
#include <iostream>
#include <armadillo>

#include <libmaw.h>
#include "triangle.h"

class triangulation
{

public:
	triangulation(matlab_engine* engine);
	~triangulation();

	//Create a delaunay triangulation
	void create_delaunay(arma::vec& x, arma::vec& y, arma::vec& z);

	//return the size of the triangluation
	size_t size();

	//set the vertex data. It is assumed that t-th triangle's global-id is an index into data
	void set_vertex_data(arma::mat& data);

	//returns the t-th triangle
	triangle& operator()(size_t t);
	
	//computer the face normals for all the triangles.
	void compute_face_normals();

	triangle* find_containing_triangle(double x, double y);

	//creates the triangulation matrix for matlab
	arma::mat matlab_tri_matrix();



private:
	//the triangles
	std::vector<triangle*> m_triangles;

	size_t m_size;  //number of triangulations

	//ptr to the matlab engine
	matlab_engine* m_engine;

};




