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

void triangulation::create_delaunay( arma::vec* x, arma::vec* y, arma::vec* z)
{
	if(m_engine)
	{

		int num_nodes = x->n_rows;

		mxArray* xy = mxCreateDoubleMatrix(num_nodes,3, mxREAL);
		double* ptr = mxGetPr(xy);

		//create temp arrays to send to matlab
		for (int row=0;row<num_nodes;row++)
		{
			//matlab is col major storage	
			ptr[row+0*num_nodes] = (*x)(row);
			ptr[row+1*num_nodes] = (*y)(row);
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
		mxArray* tri = m_engine->get("t");
	
		//will be n * 3
		const mwSize* size = mxGetDimensions(tri);
 		m_size = size[0]; //first element is the # of rows
		
		double* t = mxGetPr(tri);

		for (size_t i = 0;i<m_size; i++)
		{
			//col major lookup!!!
			//get the row index data from the matlab triangulation structure
			//note: all these indexes will be off by 1 because matlab indexing starts at 1 and not 0
				int v1 = int(t[i+0*m_size]); //v1 -= 1;
				int v2 = int(t[i+1*m_size]); //v2 -= 1;
				int v3 = int(t[i+2*m_size]); //v3 -= 1;

				point vertex1,vertex2,vertex3;
				vertex1.x = (*x)(v1-1);
				vertex1.y = (*y)(v1-1);
				vertex1.z = (*z)(v1-1);

				vertex2.x = (*x)(v2-1);
				vertex2.y = (*y)(v2-1);
				vertex2.z = (*z)(v2-1);

				vertex3.x = (*x)(v3-1);
				vertex3.y = (*y)(v3-1);
				vertex3.z = (*z)(v2-1);

				m_triangles.push_back(new triangle(vertex1,vertex2,vertex3,0));
				m_triangles[i]->global_id[0] = v1;
				m_triangles[i]->global_id[1] = v2;
				m_triangles[i]->global_id[2] = v3;

				m_triangles[i]->triangle_id = i;
		}
	
		//clean up in matlab
		m_engine->evaluate("clear t xy tri;");
	}
}

size_t triangulation::size()
{
	return m_size;
}

triangulation::triangulation( maw::matlab_engine* engine )
{
	m_engine = engine;
	m_size = 0;

}

triangulation::~triangulation()
{
	for(auto it = m_triangles.begin(); it != m_triangles.end(); it++)
	{
		delete *it;
	}
	
}

triangle& triangulation::operator()(size_t t)
{
	return *(m_triangles[t]);
}

void triangulation::set_vertex_data( arma::mat& data )
{
	if(data.n_cols != 3)
		throw std::runtime_error("Wrong number of columns");

	//loop over all triangles

	for(auto it=m_triangles.begin(); it != m_triangles.end(); ++it)
	{
		size_t v1 = (*it)->global_id[0];
		size_t v2 = (*it)->global_id[1];
		size_t v3= (*it)->global_id[2];

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

		(*it)->set_vertex_values(vertex1, vertex2, vertex3);

	}

}

void triangulation::compute_face_normals()
{

	m_engine->evaluate("[NormalVx NormalVy NormalVz PosVx PosVy PosVz]=computeNormalVectorTriangulation(mxDomain,tri,'center-cells');");
	m_engine->evaluate("normals=[NormalVx NormalVy NormalVz];");
	m_engine->evaluate("center=[PosVx PosVy PosVz]");
	m_engine->evaluate("clear PosVx PosVy PosVz NormalVx NormalVy NormalVz");

	arma::mat* normals = m_engine->get_double_matrix("normals");
	arma::mat* centers = m_engine->get_double_matrix("center");

	size_t counter = 0;
	for(auto it=m_triangles.begin();it!=m_triangles.end();it++)
	{
		arma::vec n(3);
			n(0) = normals->row(counter)(0);
			n(1) = normals->row(counter)(1);
			n(2) = normals->row(counter)(2);

		arma::vec c(3); 
			c(0) = centers->row(counter)(0);
			c(1) = centers->row(counter)(1);
			c(2) = centers->row(counter)(2);

		(*it)->set_facenormal(n);
		point p(c(0),c(1),c(2));
		(*it)->center=p;
		++counter;
	}
	m_engine->evaluate("clear normals center");

// 	for(auto it = m_triangles.begin(); it != m_triangles.end(); it++)
// 	{
// 		point vertex1 = (*it)->get_vertex(0);
// 		point vertex2 = (*it)->get_vertex(1);
// 		point vertex3 = (*it)->get_vertex(2);
// 		
// 		arma::vec AB;
// 		AB << vertex2.x - vertex1.x << vertex2.y - vertex1.y << vertex2.z - vertex1.z << arma::endr;
// 		arma::vec AC;
// 		AC << vertex3.x - vertex1.x << vertex3.y - vertex1.y << vertex3.z - vertex1.z << arma::endr;
// 
// 		arma::vec normal = cross(AB,AC);
// 
// 		//make sure we have the up facing normal
//  		if (normal(2) < 0.0)
//  			normal = -normal;
// 
// 		(*it)->set_facenormal(normal);
// 
// 	}
}

triangle* triangulation::find_containing_triangle( double x,double y )
{
	for(std::vector<triangle*>::iterator it=m_triangles.begin(); it != m_triangles.end(); it++)
	{
		if( (*it)->contains(x,y) )
			return *it;
	}

	return NULL;
}

arma::mat triangulation::matlab_tri_matrix()
{
	arma::mat tri(m_size,3);
	size_t i=0;
	for(auto it=m_triangles.begin();it!=m_triangles.end();it++)
	{
		tri(i,0) = (*it)->global_id[0];
		tri(i,1) = (*it)->global_id[1];
		tri(i,2) = (*it)->global_id[2];
		i++;
	}

	return tri;
}
