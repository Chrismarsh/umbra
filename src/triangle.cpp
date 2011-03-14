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


#include "triangle.h"

triangle::triangle( int vertex1, int vertex2, int vertex3)
{
	vertex_index[0] = vertex1;
	vertex_index[1] = vertex2;
	vertex_index[2] = vertex3;

	m_shadow = false;
	sub_tri = NULL;
}

triangle::triangle()
{
	m_shadow = false;
	sub_tri = NULL;
}
int triangle::get_vertex( int v )
{
	//add boundary check here
	return (vertex_index[v]);
}

void triangle::set_shadow( bool s )
{
	m_shadow = s;
}

bool triangle::contains( double x, double y )
{
	double x1 = *(vertex_value[0].x);
	double y1 = *(vertex_value[0].y);

	double x2 = *(vertex_value[1].y);
	double y2 = *(vertex_value[1].x);

	double x3 = *(vertex_value[2].x);
	double y3 = *(vertex_value[2].y);

	double Area_PP1P2  = 0.5 *abs(x*y1-x*y2+x1*y2-x1*y+x2*y-x2*y1);
	double Area_PP2P3  = 0.5 *abs(x*y2-x*y3+x2*y3-x2*y+x3*y-x3*y2);
	double Area_PP3P1  = 0.5 *abs(x*y3-x*y1+x3*y1-x3*y+x1*y-x1*y3);
	double Area_P1P2P3 = 0.5 *abs(x1*y2-x1*y3+x2*y3-x2*y1+x3*y1-x3*y2);

	double areasum = Area_PP1P2 + Area_PP2P3 + Area_PP3P1;
	if( abs(areasum - Area_P1P2P3) < 0.0001)
	{
		return true;
	}
	else//check subtriangle1
	{
		if(sub_tri) //need to stop the recursion
		{
			if( sub_tri[0].contains(x,y) || 
				sub_tri[1].contains(x,y) ||
				sub_tri[2].contains(x,y) ||
				sub_tri[3].contains(x,y))
				return true;
		}
		
		return false;

	}
}

void triangle::set_vertex_values( ptr_point v1, ptr_point v2, ptr_point v3 )
{
	vertex_value[0].x = v1.x;
	vertex_value[0].y = v1.y;
	vertex_value[0].z = v1.z;

	vertex_value[1].x = v2.x;
	vertex_value[1].y = v2.y;
	vertex_value[1].z = v2.z;

	vertex_value[2].x = v3.x;
	vertex_value[2].y = v3.y;
	vertex_value[2].z = v3.z;

}

point triangle::get_vertex_value( int v )
{
	point p;
	p.x = *(vertex_value[v].x);
	p.y = *(vertex_value[v].y);
	p.z = *(vertex_value[v].z);
	 
	return p;
}

void triangle::update_subtri()
{
	int longest;
 	if (sub_tri)
 		delete[] sub_tri;	
// 	

	sub_tri = new triangle[4]();

	//only one level of recursion at the moment
	for(int i = 0; i<3;i++)
		sub_tri->sub_tri=NULL;


	arma::vec l_sides(3);
	//0-1
	l_sides[0] = sqrt( (*(vertex_value[0].x)-*(vertex_value[1].x))*(*(vertex_value[0].x)-*(vertex_value[1].x)) + ((*(vertex_value[0].y)-*(vertex_value[1].y))*(*(vertex_value[0].y)-*(vertex_value[1].y))) * 1.0);
	//1-2
	l_sides[1] = sqrt( (*(vertex_value[1].x)-*(vertex_value[2].x))*(*(vertex_value[1].x)-*(vertex_value[2].x)) + ((*(vertex_value[1].y)-*(vertex_value[2].y))*(*(vertex_value[1].y)-*(vertex_value[2].y)))* 1.0);
	//0-2
	l_sides[2] = sqrt( (*(vertex_value[0].x)-*(vertex_value[2].x))*(*(vertex_value[0].x)-*(vertex_value[2].x)) + ((*(vertex_value[0].y)-*(vertex_value[2].y))*(*(vertex_value[0].y)-*(vertex_value[2].y)))* 1.0);

	if(l_sides[0] > l_sides[1] || l_sides[0] > l_sides[2])
	{
		longest = 0;
		//midpoint of the longest edge (which is 0-1)
		point *midpt_01 = new point;
		midpt_01->x = sqrt( (*(vertex_value[0].x)-*(vertex_value[1].x))*(*(vertex_value[0].x)-*(vertex_value[1].x)) + ((*(vertex_value[0].y)-*(vertex_value[1].y))*(*(vertex_value[0].y)-*(vertex_value[1].y))) * 1.0) /2.0;
		double m = (*(vertex_value[1].y)-*(vertex_value[0].y))/(*(vertex_value[1].x)-*(vertex_value[0].x));
		double b = *(vertex_value[0].y)-m**(vertex_value[0].x);
		midpt_01->y = m*midpt_01->x+b;


		//midpoint of the opposite edge (which is 0-2)
		point* midpt_02 = new point;
		midpt_02->x = sqrt( (*(vertex_value[0].x)-*(vertex_value[2].x))*(*(vertex_value[0].x)-*(vertex_value[2].x)) + ((*(vertex_value[0].y)-*(vertex_value[2].y))*(*(vertex_value[0].y)-*(vertex_value[2].y))) * 1.0) /2.0;
		m = (*(vertex_value[2].y)-*(vertex_value[0].y))/(*(vertex_value[2].x)-*(vertex_value[0].x));
		b = *(vertex_value[0].y)-m**(vertex_value[0].x);
		midpt_01->y = m*midpt_01->x+b;

		//midpoint of the opposite edge2 (which is 1-2)
		point* midpt_12 = new point;
		midpt_12->x = sqrt( (*(vertex_value[1].x)-*(vertex_value[2].x))*(*(vertex_value[1].x)-*(vertex_value[2].x)) + ((*(vertex_value[1].y)-*(vertex_value[2].y))*(*(vertex_value[1].y)-*(vertex_value[2].y))) * 1.0) /2.0;
		m = (*(vertex_value[2].y)-*(vertex_value[1].y))/(*(vertex_value[2].x)-*(vertex_value[0].x));
		b = *(vertex_value[1].y)-m**(vertex_value[1].x);
		midpt_12->y = m*midpt_12->x+b;

		//have 4 sub triangles now:
		//define CCW, left->right
		//t1 = 0-01-02-01
		//t2 = 01-2-02-01
		//t3 = 01-12-2-01
		//t4 = 01-1-12-01

		//t1 = 0-01-02-01
		ptr_point v1;
		v1.x = vertex_value[0].x;
		v1.y = vertex_value[0].y;
		ptr_point v2;
		v2.x = &(midpt_01->x);
		v2.y = &(midpt_01->y);	
		ptr_point v3;
		v3.x = &(midpt_02->x);
		v3.y = &(midpt_02->y);
		sub_tri[0].set_vertex_values(v1,v2,v3);

		//t2 = 01-2-02-01
		//ptr_point v1;
		v1.x = &(midpt_01->x);
		v1.y = &(midpt_01->y);
		//ptr_point v2;
		v2.x = vertex_value[2].x;
		v2.y = vertex_value[2].y;
		//ptr_point v3;
		v3.x = &(midpt_02->x);
		v3.y = &(midpt_02->y);
		sub_tri[1].set_vertex_values(v1,v2,v3);

		//t3 = 01-12-2-01
		//ptr_point v1;
		v1.x = &(midpt_01->x);
		v1.y = &(midpt_01->y);
		//ptr_point v2;
		v2.x = &(midpt_12->x);
		v2.y = &(midpt_12->y);
		//ptr_point v3;
		v3.x = vertex_value[2].x;
		v3.y = vertex_value[2].y;
		sub_tri[2].set_vertex_values(v1,v2,v3);

		//t4 = 01-1-12-01
		//v1;
		v1.x = &(midpt_01->x);
		v1.y = &(midpt_01->y);
		//ptr_point v2;
		v2.x = vertex_value[1].x;
		v2.y = vertex_value[1].y;
		//ptr_point v3;
		v3.x = &(midpt_12->x);
		v3.y = &(midpt_12->y);
		sub_tri[3].set_vertex_values(v1,v2,v3);

		

	}
	else if(l_sides[1] > l_sides[2] || l_sides[1] > l_sides[0])
	{
		longest = 1;

		//midpoint of the longest edge (which is 1-2)
		point* midpt_12 = new point;
		midpt_12->x = sqrt( (*(vertex_value[1].x)-*(vertex_value[2].x))*(*(vertex_value[1].x)-*(vertex_value[2].x)) + ((*(vertex_value[1].y)-*(vertex_value[2].y))*(*(vertex_value[1].y)-*(vertex_value[2].y))) * 1.0) /2.0;
		double m = (*(vertex_value[2].y)-*(vertex_value[1].y))/(*(vertex_value[2].x)-*(vertex_value[0].x));
		double b = *(vertex_value[1].y)-m**(vertex_value[1].x);
		midpt_12->y = m*midpt_12->x+b;


		//midpoint of the opposite edge (which is 0-1)
		point* midpt_01 = new point;
		midpt_01->x = sqrt( (*(vertex_value[0].x)-*(vertex_value[1].x))*(*(vertex_value[0].x)-*(vertex_value[1].x)) + ((*(vertex_value[0].y)-*(vertex_value[1].y))*(*(vertex_value[0].y)-*(vertex_value[1].y))) * 1.0) /2.0;
		m = (*(vertex_value[1].y)-*(vertex_value[0].y))/(*(vertex_value[1].x)-*(vertex_value[0].x));
		b = *(vertex_value[0].y)-m**(vertex_value[0].x);
		midpt_01->y = m*midpt_01->x+b;

		//midpoint of the opposite edge2 (which is 0-2)
		point* midpt_02 = new point;
		midpt_02->x = sqrt( (*(vertex_value[0].x)-*(vertex_value[2].x))*(*(vertex_value[0].x)-*(vertex_value[2].x)) + ((*(vertex_value[0].y)-*(vertex_value[2].y))*(*(vertex_value[0].y)-*(vertex_value[2].y))) * 1.0) /2.0;
		m = (*(vertex_value[2].y)-*(vertex_value[0].y))/(*(vertex_value[2].x)-*(vertex_value[0].x));
		b = *(vertex_value[0].y)-m**(vertex_value[0].x);
		midpt_02->y = m*midpt_02->x+b;

		//have 4 sub triangles now:
		//define CCW, left->right
		//t1 = 0-01-12-0
		//t2 = 01-1-12-01
		//t3 = 12-2-02-12
		//t4 = 0-12-02-0

		//t1 = 0-01-12-0
		ptr_point v1;
		v1.x = vertex_value[0].x;
		v1.y = vertex_value[0].y;
		ptr_point v2;
		v2.x = &(midpt_01->x);
		v2.y = &(midpt_01->y);	
		ptr_point v3;
		v3.x = &(midpt_12->x);
		v3.y = &(midpt_12->y);
		sub_tri[0].set_vertex_values(v1,v2,v3);

		//t2 = 01-1-12-01
//		ptr_point v1;
		v1.x = &(midpt_01->x);
		v1.y = &(midpt_01->y);
//		ptr_point v2;
		v2.x = vertex_value[1].x;
		v2.y = vertex_value[1].y;	
//		ptr_point v3;
		v3.x = &(midpt_12->x);
		v3.y = &(midpt_12->y);
		sub_tri[1].set_vertex_values(v1,v2,v3);

		//t3 = 12-2-02-12
//		ptr_point v1;
		v1.x = &(midpt_12->x);
		v1.y = &(midpt_12->y);
//		ptr_point v2;
		v2.x = vertex_value[2].x;
		v2.y = vertex_value[2].y;	
//		ptr_point v3;
		v3.x = &(midpt_02->x);
		v3.y = &(midpt_02->y);
		sub_tri[2].set_vertex_values(v1,v2,v3);

		//t4 = 0-12-02-0
//		ptr_point v1;
		v1.x = vertex_value[0].x;
		v1.y = vertex_value[0].y;
//		ptr_point v2;
		v2.x = &(midpt_12->x);
		v2.y = &(midpt_12->y);	
//		ptr_point v3;
		v3.x = &(midpt_02->x);
		v3.y = &(midpt_02->y);
		sub_tri[3].set_vertex_values(v1,v2,v3);


	}
	else 
	{
		longest = 2;

		//midpoint of the longest edge (which is 0-2)
		point* midpt_02 = new point;
		midpt_02->x = sqrt( (*(vertex_value[0].x)-*(vertex_value[2].x))*(*(vertex_value[0].x)-*(vertex_value[1].x)) + ((*(vertex_value[0].y)-*(vertex_value[2].y))*(*(vertex_value[0].y)-*(vertex_value[2].y))) * 1.0) /2.0;
		double m = (*(vertex_value[2].y)-*(vertex_value[0].y))/(*(vertex_value[2].x)-*(vertex_value[0].x));
		double b = *(vertex_value[0].y)-m**(vertex_value[0].x);
		midpt_02->y = m*midpt_02->x+b;


		//midpoint of the opposite edge (which is 0-1)
		point* midpt_01 = new point;
		midpt_01->x = sqrt( (*(vertex_value[0].x)-*(vertex_value[2].x))*(*(vertex_value[0].x)-*(vertex_value[1].x)) + ((*(vertex_value[0].y)-*(vertex_value[1].y))*(*(vertex_value[0].y)-*(vertex_value[1].y))) * 1.0) /2.0;
		m = (*(vertex_value[1].y)-*(vertex_value[0].y))/(*(vertex_value[1].x)-*(vertex_value[0].x));
		b = *(vertex_value[0].y)-m**(vertex_value[0].x);
		midpt_01->y = m*midpt_01->x+b;

		//midpoint of the opposite edge2 (which is 1-2)
		point* midpt_12 = new point;
		midpt_12->x = sqrt( (*(vertex_value[1].x)-*(vertex_value[2].x))*(*(vertex_value[1].x)-*(vertex_value[2].x)) + ((*(vertex_value[1].y)-*(vertex_value[2].y))*(*(vertex_value[1].y)-*(vertex_value[2].y))) * 1.0) /2.0;
		m = (*(vertex_value[2].y)-*(vertex_value[1].y))/(*(vertex_value[2].x)-*(vertex_value[1].x));
		b = *(vertex_value[1].y)-m**(vertex_value[1].x);
		midpt_12->y = m*midpt_12->x+b;

		//have 4 sub triangles now:
		//define CCW, left->right
		//t1 = 0-01-02
		//t2 = 01 -1 -02
		//t3 = 1 - 12-02
		//t4 = 12 - 2 -02

		//t1 =0-01-02
		ptr_point v1;
		v1.x = vertex_value[0].x;
		v1.y = vertex_value[0].y;
		ptr_point v2;
		v2.x = &(midpt_01->x);
		v2.y = &(midpt_01->y);	
		ptr_point v3;
		v3.x = &(midpt_02->x);
		v3.y = &(midpt_02->y);
		sub_tri[0].set_vertex_values(v1,v2,v3);

		//t2 = 01 -1 -02
//		ptr_point v1;
		v1.x = &(midpt_01->x);
		v1.y = &(midpt_01->y);
//		ptr_point v2;
		v2.x = vertex_value[1].x;
		v2.y = vertex_value[1].y;
//		ptr_point v3;
		v3.x = &(midpt_02->x);
		v3.y = &(midpt_02->y);
		sub_tri[1].set_vertex_values(v1,v2,v3);

		//t3 = 1 - 12-02
//		ptr_point v1;
		v1.x = vertex_value[1].x;
		v1.y = vertex_value[1].y;
//		ptr_point v2;
		v2.x = &(midpt_12->x);
		v2.y = &(midpt_12->y);
//		ptr_point v3;
		v3.x = &(midpt_02->x);
		v3.y = &(midpt_02->y);
		sub_tri[2].set_vertex_values(v1,v2,v3);

		//t4 = 12 - 2 -02
//		ptr_point v1;
		v1.x = &(midpt_12->x);
		v1.y = &(midpt_12->y);
//		ptr_point v2;
		v2.x = vertex_value[2].x;
		v2.y = vertex_value[2].y;
//		ptr_point v3;
		v3.x = &(midpt_02->x);
		v3.y = &(midpt_02->y);
		sub_tri[3].set_vertex_values(v1,v2,v3);

	}	
		



}

bool triangle::intersects( triangle* t )
{
	//for each sub-triangle of t
// 	for(int i = 0; i<4; i++)
// 	{
// 		if (this->contains(t->get_vertex_value(0).x,t->get_vertex_value(0).y) ||
// 			this->contains(t->get_vertex_value(1).x,t->get_vertex_value(1).y) ||
// 			this->contains(t->get_vertex_value(2).x,t->get_vertex_value(2).y) ||
// 			this->contains(t->sub_tri[0].get_vertex_value(0).x,t->sub_tri[0].get_vertex_value(0).y) ||
// 
// 	}
return false;
}