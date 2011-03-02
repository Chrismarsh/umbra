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
		m_engine->evaluate("tri=DelaunayTri(xyz(:,1),xyz(:,2)");
		
		//clean up our temp array
		mxDestroyArray(xyz);
		xyz=NULL;
		
		//get our triangulation structure from matlab
		mxArray* tri = m_engine->get("tri.Triangulation");

		//will be n * 3
		const mwSize* size = mxGetDimensions(tri);
		m_size = size[0]; //first element is the # of rows
		
		double* t = mxGetPr(tri);

		for (int i = 0;i<m_size; i++)
		{
				m_triangles.push_back(new triangle(int(t[i*3+0]),int(t[i*3+1]),int(t[i*3+2])));
				//triangle centres?S
		}
	

	}
}