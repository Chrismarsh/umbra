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


#include "triangulation.h"

void triangulation::create_delaunay(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z )
{
	if(m_engine)
	{

		int num_nodes = x.size();

		mxArray* xyz = mxCreateDoubleMatrix(num_nodes,3, mxREAL);
		double* ptr = mxGetPr(xyz);

		//create temp arrays to send to matlab
		for (int i=0;i<num_nodes;i++)
		{
			ptr[i*3+0] = x[i];
			ptr[i*3+1] = y[i];
			ptr[i*3+2] = z[i];	
		}
		m_engine->put("xyz",xyz);
		m_engine->evaluate("tri=DelaunayTri(xyz(:,1),xyz(:,2))");
		
		//clean up our temp array
		mxDestroyArray(xyz);
		xyz=NULL;
		
		//for some reason grabbing triangulation directly doesn't work so throw it in a temp var
		m_engine->evaluate("t=tri.Triangulation");
		//get our triangulation structure from matlab
		mxArray* tri = m_engine->get("t");
		if(!tri)
			throw std::exception(m_engine->get_last_error().c_str());

		//will be n * 3
		const mwSize* size = mxGetDimensions(tri);
		m_size = size[0]; //first element is the # of rows
		
		double* t = mxGetPr(tri);

		for (int i = 0;i<m_size; i++)
		{
				int v1=int(t[i*3+0]);
				int v2 = int(t[i*3+1]);
				int v3 = int(t[i*3+2]);

				m_triangles.push_back(new triangle(v1,v2,v3));
				//triangle centres?S
		}
	

	}
}

int triangulation::get_size()
{
	return m_size;
}

triangulation::triangulation( matlab* engine )
{
	m_engine = engine;
	m_size = 0;
}

triangulation::~triangulation()
{
	std::vector<triangle*>::iterator it = m_triangles.begin();

	for(; it != m_triangles.end(); it++)
	{
		delete *it;
	}

}

std::vector<int> triangulation::get_tri( int t )
{
	std::vector<int> v;
	
	v.push_back(m_triangles[t]->get_vertex(0));
	v.push_back(m_triangles[t]->get_vertex(1));
	v.push_back(m_triangles[t]->get_vertex(2));

	return v;
}
