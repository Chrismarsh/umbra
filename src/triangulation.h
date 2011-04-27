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
//matlab engine
#include <engine.h>
#include <armadillo>

#include "matlab_engine.h"
#include "triangle.h"

class triangulation
{

public:
	triangulation(matlab* engine);
	~triangulation();
	void create_delaunay(arma::vec& x, arma::vec& y, arma::vec& z);
	size_t get_num_tri();
	void set_vertex_data(arma::mat& data);
	arma::uvec get_index(size_t t);
	//returns the t-th triangle
	triangle& operator()(size_t t);
	
	//get the matlab compliant triangulation index that contains indexes into the data array
	arma::umat& get_tri_index();
	void compute_face_normals();

	triangle* find_containing_triangle(double x, double y);

	



private:



	//these two structures need to stay in sync so that matlab can plot them
	//the triangles
	std::vector<triangle*> m_triangles;
	//triangulation, stored as a list of indexes into the x,y,z data. This is like Matlab
	arma::umat* m_tri;  //size * 3
	size_t m_size;  //number of triangulations

	//ptr to the matlab engine
	matlab* m_engine;

};




