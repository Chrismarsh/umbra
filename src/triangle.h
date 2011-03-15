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

#include "point.h"
#include "matlab_engine.h"

class triangle
{
public:
	//use xyz triples in the vector
	//store the index like matlab
	triangle( ptr_point vertex1, ptr_point vertex2, ptr_point vertex3, size_t cur_rec_depth=0);
	triangle(size_t cur_rec_depth);


	//returns the v-th vertex[0-2] of the triangle
	point operator()(size_t v);

	void update_subtri();
	triangle& sub_tri(size_t t);

	void set_vertex_values( ptr_point vertex1, ptr_point vertex2, ptr_point vertex3);
	bool contains(double x, double y);
	bool contains(point xy);
	bool intersects(triangle* t);
	point get_center();
private:
	ptr_point m_vertex_list[3];
	
	point m_center;
	//always 4 subriangles
	triangle** m_sub_tri;

	//need this so as to be able to call matlab helper functions
	//matlab* m_engine;

	//set the number of sub triangles
	size_t m_cur_rec_depth;
};