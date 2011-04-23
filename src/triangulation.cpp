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

void triangulation::create_delaunay( arma::vec& x, arma::vec& y, arma::vec& z)
{
	if(m_engine)
	{
	
		int num_nodes = x.n_rows;

		mxArray* xy = mxCreateDoubleMatrix(num_nodes,3, mxREAL);
		double* ptr = mxGetPr(xy);

		//create temp arrays to send to matlab
		for (int row=0;row<num_nodes;row++)
		{
			//matlab is col major storage	
			ptr[row+0*num_nodes] = x(row);
			ptr[row+1*num_nodes] = y(row);
		}
		m_engine->put("xy",xy);
		m_engine->evaluate("tri=DelaunayTri(xy(:,1),xy(:,2))");

		//clean up our temp array
		mxDestroyArray(xy);
		xy=NULL;
		ptr=NULL;
		 
		//change this later to the struct lookup
		m_engine->evaluate("t=tri.Triangulation");
		//get our triangulation structure from matlab
		mxArray* tri = NULL;
		tri = m_engine->get("t");
	

		//will be n * 3
		const mwSize* size = mxGetDimensions(tri);
 		m_size = size[0]; //first element is the # of rows
		
		m_triangles.resize(m_size);

		double* t = mxGetPr(tri);
		m_tri = new arma::umat(m_size,3);
		for (size_t i = 0;i<m_size; i++)
		{
			//col major lookup!!!
			//get the row index data from the matlab triangulation structure
			//note: all these indexes will be off by 1 because matlab indexing starts at 1 and not 0
				int v1 = int(t[i+0*m_size]); //v1 -= 1;
				int v2 = int(t[i+1*m_size]); //v2 -= 1;
				int v3 = int(t[i+2*m_size]); //v3 -= 1;

				m_tri->operator()(i,0) = v1;
				m_tri->operator()(i,1) = v2;
				m_tri->operator()(i,2) = v3;

				point vertex1,vertex2,vertex3;
				vertex1.x = x(v1-1);
				vertex1.y = y(v1-1);
				vertex1.z = z(v1-1);

				vertex2.x = x(v2-1);
				vertex2.y = y(v2-1);
				vertex2.z = z(v2-1);

				vertex3.x = x(v3-1);
				vertex3.y = y(v3-1);
				vertex3.z = z(v2-1);

				m_triangles[i] = new triangle(vertex1,vertex2,vertex3);
				//m_triangles.push_back(new triangle(vertex1,vertex2,vertex3));
		//std::cout << i << std::endl;
		}
	
		//clean up in matlab
		m_engine->evaluate("clear t xy tri;");

	}
}

size_t triangulation::get_num_tri()
{
	return m_size;
}

triangulation::triangulation( matlab* engine )
{
	m_engine = engine;
	m_size = 0;
	m_tri = NULL ;
}

triangulation::~triangulation()
{
	for(auto it = m_triangles.begin(); it != m_triangles.end(); it++)
	{
		delete *it;
	}

	delete m_tri;
	m_engine = NULL;

	
}

arma::uvec triangulation::get_index( size_t t )
{
	arma::uvec v(3);
	
	v(0) = m_tri->operator()(t,0);
	v(1) = m_tri->operator()(t,1);
	v(2) = m_tri->operator()(t,2);

	return v;
}

triangle& triangulation::operator()(size_t t)
{
	return *(m_triangles[t]);
}

void triangulation::set_vertex_data( arma::mat& data )
{
	if(data.n_cols != 3)
		throw std::exception("Wrong number of columns");

	for(size_t i = 0;i<m_size;i++)
	{
		size_t v1 = m_tri->operator()(i,0);
		size_t v2 = m_tri->operator()(i,1);
		size_t v3 = m_tri->operator()(i,2);

		point vertex1,vertex2,vertex3;
		vertex1.x = data(v1-1,0);
		vertex1.y = data(v1-1,1);
		vertex1.z = data(v1-1,2);

		vertex2.x = data(v2-1,0);
		vertex2.y = data(v2-1,1);
		vertex2.z = data(v2-1,2);

		vertex3.x = data(v3-1,0);
		vertex3.y = data(v3-1,1);
		vertex3.z = data(v3-1,2);


		m_triangles[i]->set_vertex_values(vertex1, vertex2, vertex3);


	}
}

arma::umat& triangulation::get_tri_index()
{
	return *m_tri;
}

void triangulation::compute_face_normals()
{

	for(auto it = m_triangles.begin(); it != m_triangles.end(); it++)
	{
		triangle* t = *it;
		point vertex1 = t->get_vertex_value(0);
		point vertex2 = t->get_vertex_value(1);
		point vertex3 = t->get_vertex_value(2);
		
		arma::vec AB;
		AB << vertex2.x - vertex1.x << vertex2.y - vertex1.y << vertex2.z - vertex1.z << arma::endr;
		arma::vec AC;
		AC << vertex3.x - vertex1.x << vertex3.y - vertex1.y << vertex3.z - vertex1.z << arma::endr;

		arma::vec normal = cross(AB,AC);
		t->set_facenormal(normal);

	}



}

