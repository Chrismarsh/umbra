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

#include <armadillo>
#include <vector>

#include "triangle.h"
#include <libmaw.h>

// 
// %OBB computation 
// 	     -----------------------
// 
// 	         4       top      3
// 	         +--------+--------+
// 	         |        |        |
// 	         |        |        |
// 	   left  +--------+--------+    right
// 	         |        |        |
// 	         |        |        |
// 	         +-----------------+
// 	         1     bottom      2
// 	 +y
// 	 ^
// 	 |
// 	 +-> +x

class rect
{
public:
	rect( arma::mat* coord )
	{
		this->coord = coord;
	}
	~rect()
	{
		delete coord;
	}
	arma::mat* coord;
	std::vector<triangle*> triangles;
	std::vector<int> m_globalID;
};
class bounding_rect
{
public:
	bounding_rect(matlab_engine* m_engine);
	~bounding_rect();
	void make(const arma::vec* x, const arma::vec* y, int n_rows, int n_cols);
	rect* bounding_rect::get_rect( int i, int j );
	int n_rows;
	int n_cols;
	bool pt_in_rect(double x, double y, rect* r);

private:
	//std::vector<rect*> m_rectangles;
	std::vector<std::vector<rect*> > m_grid;
	matlab_engine* m_engine;
};


