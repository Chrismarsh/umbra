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

triangle::triangle( ptr_point vertex1, ptr_point vertex2, ptr_point vertex3,size_t cur_rec_depth/*=0*/)
{

	m_cur_rec_depth = cur_rec_depth+1;
// 	m_sub_tri[0] = NULL;
// 	m_sub_tri[1] = NULL;
// 	m_sub_tri[2] = NULL;
// 	m_sub_tri[3] = NULL;

	set_vertex_values(vertex1,vertex2, vertex3);

	if(cur_rec_depth == 0)
		update_subtri();
	else
	{
		for(int i = 0;i<4;i++)
			m_sub_tri[i] = NULL;
	}
}

triangle::triangle(size_t cur_rec_depth)
{
	m_sub_tri = NULL;

	m_cur_rec_depth = cur_rec_depth+1;


}



bool triangle::contains( double x, double y )
{
	double x1 = *(m_vertex_list[0].x);
	double y1 = *(m_vertex_list[0].y);

	double x2 = *(m_vertex_list[1].y);
	double y2 = *(m_vertex_list[1].x);

	double x3 = *(m_vertex_list[2].x);
	double y3 = *(m_vertex_list[2].y);

	double Area_PP1P2  = 0.5 *abs(x*y1-x*y2+x1*y2-x1*y+x2*y-x2*y1);
	double Area_PP2P3  = 0.5 *abs(x*y2-x*y3+x2*y3-x2*y+x3*y-x3*y2);
	double Area_PP3P1  = 0.5 *abs(x*y3-x*y1+x3*y1-x3*y+x1*y-x1*y3);
	double Area_P1P2P3 = 0.5 *abs(x1*y2-x1*y3+x2*y3-x2*y1+x3*y1-x3*y2);

	double areasum = Area_PP1P2 + Area_PP2P3 + Area_PP3P1;
	if( abs(areasum - Area_P1P2P3) < 0.0001)
		return true;
	else
		return false;

}

bool triangle::contains( point xy )
{
	return contains(xy.x,xy.y);
}

void triangle::set_vertex_values( ptr_point vertex1, ptr_point vertex2,ptr_point vertex3)
{
	m_vertex_list[0].x = vertex1.x;
	m_vertex_list[0].y = vertex1.y;
//	m_vertex_list[0].z = vertex1.z;

	m_vertex_list[1].x = vertex2.x;
	m_vertex_list[1].y = vertex2.y;
//	m_vertex_list[1].z = vertex2.z;

	m_vertex_list[2].x = vertex3.x;
	m_vertex_list[2].y = vertex3.y;
//	m_vertex_list[2].z = vertex3.z;

// 	arma::mat c(3,3);
// 	c << *(m_vertex_list[0].x) << *(m_vertex_list[0].y) << 0 << arma::endr
// 	  << *(m_vertex_list[1].x) << *(m_vertex_list[1].y) << 0 << arma::endr
// 	  << *(m_vertex_list[2].x) << *(m_vertex_list[2].y) << 0 << arma::endr;
	

	//ok lets try a faster circumcenter algorithm *not* implimented in matlab :/

	arma::vec Pa(3);
	Pa(0) = *(m_vertex_list[0].x);
	Pa(1) = *(m_vertex_list[0].y);
	Pa(2) = 0;

	arma::vec Pb(3);
	Pb(0) = *(m_vertex_list[1].x);
	Pb(1) = *(m_vertex_list[1].y);
	Pb(2) = 0;

	arma::vec Pc(3);
	Pc(0) = *(m_vertex_list[2].x);
	Pc(1) = *(m_vertex_list[2].y);
	Pc(2) = 0;

	arma::vec AB = Pb - Pa;
	arma::vec AC = Pc - Pa;
	arma::vec BC = Pc - Pb;

	arma::vec N = arma::cross(AC,AB);
	arma::vec L1 = arma::cross(AB,N);
	arma::vec L2 = arma::cross(BC,N);
	arma::vec P21 = (Pc - Pa)/2;
	arma::vec P1 = (Pa + Pb)/2;

	//temp
	arma::mat ML1(L1);
	arma::mat ML2(L2);

	arma::mat ML = arma::join_rows(ML1,-ML2);

	arma::vec lambda = arma::solve(ML,P21);
		
	arma::vec pos = P1+lambda(0)*L1;
		

// 
// 	m_engine->put_double_matrix("t_set",&c);
// 	m_engine->evaluate("c=tri_center( [t_set(1,1) t_set(1,2) t_set(1,3)],[t_set(2,1) t_set(2,2) t_set(2,3)],[t_set(3,1) t_set(3,2) t_set(3,3)],'circumcenter')");
//	arma::vec* pos = m_engine->get_double_vector("c");

	m_center.x = (pos)(0);
	m_center.y = (pos)(1);
//	m_center.z = (*pos)(2);
}

void triangle::update_subtri()
{
	int longest;
//  	if (m_sub_tri)
//  		delete[] m_sub_tri;	
	

	//only one level of recursion at the moment
	// set each sub
	m_sub_tri = new triangle*[4];
	for(int i = 0; i<4;i++)
	{
		m_sub_tri[i] = new triangle(m_cur_rec_depth);
		

		//this->m_sub_tri[i]->m_cur_rec_depth = this->m_cur_rec_depth+1;
	}


	arma::vec l_sides(3);
	//0-1
	l_sides[0] = sqrt( (*(m_vertex_list[0].x)-*(m_vertex_list[1].x))*(*(m_vertex_list[0].x)-*(m_vertex_list[1].x)) + ((*(m_vertex_list[0].y)-*(m_vertex_list[1].y))*(*(m_vertex_list[0].y)-*(m_vertex_list[1].y))) * 1.0);
	//1-2
	l_sides[1] = sqrt( (*(m_vertex_list[1].x)-*(m_vertex_list[2].x))*(*(m_vertex_list[1].x)-*(m_vertex_list[2].x)) + ((*(m_vertex_list[1].y)-*(m_vertex_list[2].y))*(*(m_vertex_list[1].y)-*(m_vertex_list[2].y)))* 1.0);
	//0-2
	l_sides[2] = sqrt( (*(m_vertex_list[0].x)-*(m_vertex_list[2].x))*(*(m_vertex_list[0].x)-*(m_vertex_list[2].x)) + ((*(m_vertex_list[0].y)-*(m_vertex_list[2].y))*(*(m_vertex_list[0].y)-*(m_vertex_list[2].y)))* 1.0);

	if(l_sides[0] > l_sides[1] || l_sides[0] > l_sides[2])
	{
		longest = 0;
		//midpoint of the longest edge (which is 0-1)
		point *midpt_01 = new point;
		midpt_01->x = sqrt( (*(m_vertex_list[0].x)-*(m_vertex_list[1].x))*(*(m_vertex_list[0].x)-*(m_vertex_list[1].x)) + ((*(m_vertex_list[0].y)-*(m_vertex_list[1].y))*(*(m_vertex_list[0].y)-*(m_vertex_list[1].y))) * 1.0) /2.0;
		double m = (*(m_vertex_list[1].y)-*(m_vertex_list[0].y))/(*(m_vertex_list[1].x)-*(m_vertex_list[0].x));
		double b = *(m_vertex_list[0].y)-m**(m_vertex_list[0].x);
		midpt_01->y = m*midpt_01->x+b;


		//midpoint of the opposite edge (which is 0-2)
		point* midpt_02 = new point;
		midpt_02->x = sqrt( (*(m_vertex_list[0].x)-*(m_vertex_list[2].x))*(*(m_vertex_list[0].x)-*(m_vertex_list[2].x)) + ((*(m_vertex_list[0].y)-*(m_vertex_list[2].y))*(*(m_vertex_list[0].y)-*(m_vertex_list[2].y))) * 1.0) /2.0;
		m = (*(m_vertex_list[2].y)-*(m_vertex_list[0].y))/(*(m_vertex_list[2].x)-*(m_vertex_list[0].x));
		b = *(m_vertex_list[0].y)-m**(m_vertex_list[0].x);
		midpt_01->y = m*midpt_01->x+b;

		//midpoint of the opposite edge2 (which is 1-2)
		point* midpt_12 = new point;
		midpt_12->x = sqrt( (*(m_vertex_list[1].x)-*(m_vertex_list[2].x))*(*(m_vertex_list[1].x)-*(m_vertex_list[2].x)) + ((*(m_vertex_list[1].y)-*(m_vertex_list[2].y))*(*(m_vertex_list[1].y)-*(m_vertex_list[2].y))) * 1.0) /2.0;
		m = (*(m_vertex_list[2].y)-*(m_vertex_list[1].y))/(*(m_vertex_list[2].x)-*(m_vertex_list[0].x));
		b = *(m_vertex_list[1].y)-m**(m_vertex_list[1].x);
		midpt_12->y = m*midpt_12->x+b;

		//have 4 sub triangles now:
		//define CCW, left->right
		//t1 = 0-01-02-01
		//t2 = 01-2-02-01
		//t3 = 01-12-2-01
		//t4 = 01-1-12-01

		//t1 = 0-01-02-01
		ptr_point v1;
		v1.x = m_vertex_list[0].x;
		v1.y = m_vertex_list[0].y;
		ptr_point v2;
		v2.x = &(midpt_01->x);
		v2.y = &(midpt_01->y);	
		ptr_point v3;
		v3.x = &(midpt_02->x);
		v3.y = &(midpt_02->y);
		m_sub_tri[0]->set_vertex_values(v1,v2,v3);

		//t2 = 01-2-02-01
		//ptr_point v1;
		v1.x = &(midpt_01->x);
		v1.y = &(midpt_01->y);
		//ptr_point v2;
		v2.x = m_vertex_list[2].x;
		v2.y = m_vertex_list[2].y;
		//ptr_point v3;
		v3.x = &(midpt_02->x);
		v3.y = &(midpt_02->y);
		m_sub_tri[1]->set_vertex_values(v1,v2,v3);

		//t3 = 01-12-2-01
		//ptr_point v1;
		v1.x = &(midpt_01->x);
		v1.y = &(midpt_01->y);
		//ptr_point v2;
		v2.x = &(midpt_12->x);
		v2.y = &(midpt_12->y);
		//ptr_point v3;
		v3.x = m_vertex_list[2].x;
		v3.y = m_vertex_list[2].y;
		m_sub_tri[2]->set_vertex_values(v1,v2,v3);

		//t4 = 01-1-12-01
		//v1;
		v1.x = &(midpt_01->x);
		v1.y = &(midpt_01->y);
		//ptr_point v2;
		v2.x = m_vertex_list[1].x;
		v2.y = m_vertex_list[1].y;
		//ptr_point v3;
		v3.x = &(midpt_12->x);
		v3.y = &(midpt_12->y);
		m_sub_tri[3]->set_vertex_values(v1,v2,v3);

		

	}
	else if(l_sides[1] > l_sides[2] || l_sides[1] > l_sides[0])
	{
		longest = 1;

		//midpoint of the longest edge (which is 1-2)
		point* midpt_12 = new point;
		midpt_12->x = sqrt( (*(m_vertex_list[1].x)-*(m_vertex_list[2].x))*(*(m_vertex_list[1].x)-*(m_vertex_list[2].x)) + ((*(m_vertex_list[1].y)-*(m_vertex_list[2].y))*(*(m_vertex_list[1].y)-*(m_vertex_list[2].y))) * 1.0) /2.0;
		double m = (*(m_vertex_list[2].y)-*(m_vertex_list[1].y))/(*(m_vertex_list[2].x)-*(m_vertex_list[0].x));
		double b = *(m_vertex_list[1].y)-m**(m_vertex_list[1].x);
		midpt_12->y = m*midpt_12->x+b;


		//midpoint of the opposite edge (which is 0-1)
		point* midpt_01 = new point;
		midpt_01->x = sqrt( (*(m_vertex_list[0].x)-*(m_vertex_list[1].x))*(*(m_vertex_list[0].x)-*(m_vertex_list[1].x)) + ((*(m_vertex_list[0].y)-*(m_vertex_list[1].y))*(*(m_vertex_list[0].y)-*(m_vertex_list[1].y))) * 1.0) /2.0;
		m = (*(m_vertex_list[1].y)-*(m_vertex_list[0].y))/(*(m_vertex_list[1].x)-*(m_vertex_list[0].x));
		b = *(m_vertex_list[0].y)-m**(m_vertex_list[0].x);
		midpt_01->y = m*midpt_01->x+b;

		//midpoint of the opposite edge2 (which is 0-2)
		point* midpt_02 = new point;
		midpt_02->x = sqrt( (*(m_vertex_list[0].x)-*(m_vertex_list[2].x))*(*(m_vertex_list[0].x)-*(m_vertex_list[2].x)) + ((*(m_vertex_list[0].y)-*(m_vertex_list[2].y))*(*(m_vertex_list[0].y)-*(m_vertex_list[2].y))) * 1.0) /2.0;
		m = (*(m_vertex_list[2].y)-*(m_vertex_list[0].y))/(*(m_vertex_list[2].x)-*(m_vertex_list[0].x));
		b = *(m_vertex_list[0].y)-m**(m_vertex_list[0].x);
		midpt_02->y = m*midpt_02->x+b;

		//have 4 sub triangles now:
		//define CCW, left->right
		//t1 = 0-01-12-0
		//t2 = 01-1-12-01
		//t3 = 12-2-02-12
		//t4 = 0-12-02-0

		//t1 = 0-01-12-0
		ptr_point v1;
		v1.x = m_vertex_list[0].x;
		v1.y = m_vertex_list[0].y;
		ptr_point v2;
		v2.x = &(midpt_01->x);
		v2.y = &(midpt_01->y);	
		ptr_point v3;
		v3.x = &(midpt_12->x);
		v3.y = &(midpt_12->y);
		m_sub_tri[0]->set_vertex_values(v1,v2,v3);

		//t2 = 01-1-12-01

		v1.x = &(midpt_01->x);
		v1.y = &(midpt_01->y);

		v2.x = m_vertex_list[1].x;
		v2.y = m_vertex_list[1].y;	

		v3.x = &(midpt_12->x);
		v3.y = &(midpt_12->y);
		m_sub_tri[1]->set_vertex_values(v1,v2,v3);

		//t3 = 12-2-02-12

		v1.x = &(midpt_12->x);
		v1.y = &(midpt_12->y);

		v2.x = m_vertex_list[2].x;
		v2.y = m_vertex_list[2].y;	

		v3.x = &(midpt_02->x);
		v3.y = &(midpt_02->y);
		m_sub_tri[2]->set_vertex_values(v1,v2,v3);

		//t4 = 0-12-02-0

		v1.x = m_vertex_list[0].x;
		v1.y = m_vertex_list[0].y;

		v2.x = &(midpt_12->x);
		v2.y = &(midpt_12->y);	

		v3.x = &(midpt_02->x);
		v3.y = &(midpt_02->y);
		m_sub_tri[3]->set_vertex_values(v1,v2,v3);


	}
	else 
	{
		longest = 2;

		//midpoint of the longest edge (which is 0-2)
		point* midpt_02 = new point;
		midpt_02->x = sqrt( (*(m_vertex_list[0].x)-*(m_vertex_list[2].x))*(*(m_vertex_list[0].x)-*(m_vertex_list[1].x)) + ((*(m_vertex_list[0].y)-*(m_vertex_list[2].y))*(*(m_vertex_list[0].y)-*(m_vertex_list[2].y))) * 1.0) /2.0;
		double m = (*(m_vertex_list[2].y)-*(m_vertex_list[0].y))/(*(m_vertex_list[2].x)-*(m_vertex_list[0].x));
		double b = *(m_vertex_list[0].y)-m**(m_vertex_list[0].x);
		midpt_02->y = m*midpt_02->x+b;


		//midpoint of the opposite edge (which is 0-1)
		point* midpt_01 = new point;
		midpt_01->x = sqrt( (*(m_vertex_list[0].x)-*(m_vertex_list[2].x))*(*(m_vertex_list[0].x)-*(m_vertex_list[1].x)) + ((*(m_vertex_list[0].y)-*(m_vertex_list[1].y))*(*(m_vertex_list[0].y)-*(m_vertex_list[1].y))) * 1.0) /2.0;
		m = (*(m_vertex_list[1].y)-*(m_vertex_list[0].y))/(*(m_vertex_list[1].x)-*(m_vertex_list[0].x));
		b = *(m_vertex_list[0].y)-m**(m_vertex_list[0].x);
		midpt_01->y = m*midpt_01->x+b;

		//midpoint of the opposite edge2 (which is 1-2)
		point* midpt_12 = new point;
		midpt_12->x = sqrt( (*(m_vertex_list[1].x)-*(m_vertex_list[2].x))*(*(m_vertex_list[1].x)-*(m_vertex_list[2].x)) + ((*(m_vertex_list[1].y)-*(m_vertex_list[2].y))*(*(m_vertex_list[1].y)-*(m_vertex_list[2].y))) * 1.0) /2.0;
		m = (*(m_vertex_list[2].y)-*(m_vertex_list[1].y))/(*(m_vertex_list[2].x)-*(m_vertex_list[1].x));
		b = *(m_vertex_list[1].y)-m**(m_vertex_list[1].x);
		midpt_12->y = m*midpt_12->x+b;

		//have 4 sub triangles now:
		//define CCW, left->right
		//t1 = 0-01-02
		//t2 = 01 -1 -02
		//t3 = 1 - 12-02
		//t4 = 12 - 2 -02

		//t1 =0-01-02
		ptr_point v1;
		v1.x = m_vertex_list[0].x;
		v1.y = m_vertex_list[0].y;
		ptr_point v2;
		v2.x = &(midpt_01->x);
		v2.y = &(midpt_01->y);	
		ptr_point v3;
		v3.x = &(midpt_02->x);
		v3.y = &(midpt_02->y);
		m_sub_tri[0]->set_vertex_values(v1,v2,v3);

		//t2 = 01 -1 -02

		v1.x = &(midpt_01->x);
		v1.y = &(midpt_01->y);

		v2.x = m_vertex_list[1].x;
		v2.y = m_vertex_list[1].y;

		v3.x = &(midpt_02->x);
		v3.y = &(midpt_02->y);
		m_sub_tri[1]->set_vertex_values(v1,v2,v3);

		//t3 = 1 - 12-02

		v1.x = m_vertex_list[1].x;
		v1.y = m_vertex_list[1].y;

		v2.x = &(midpt_12->x);
		v2.y = &(midpt_12->y);

		v3.x = &(midpt_02->x);
		v3.y = &(midpt_02->y);
		m_sub_tri[2]->set_vertex_values(v1,v2,v3);

		//t4 = 12 - 2 -02
		v1.x = &(midpt_12->x);
		v1.y = &(midpt_12->y);

		v2.x = m_vertex_list[2].x;
		v2.y = m_vertex_list[2].y;

		v3.x = &(midpt_02->x);
		v3.y = &(midpt_02->y);
		m_sub_tri[3]->set_vertex_values(v1,v2,v3);

	}	
		



}
triangle& triangle::sub_tri(size_t t)
{
	return *(m_sub_tri[t]);
}

bool triangle::intersects( triangle* t )
{
	bool intersec=false;
	if (m_sub_tri)
	{
		//for each sub-triangle of t
		for(int i =0 ; i<4;i++)
		{
			if ( intersec == false && m_sub_tri[i]->intersects(t) == true)
				intersec = true;
		}
	}
	else
	{
		if(this->contains(t->get_center()))
			return true;
		else
			return false;
	}

	return intersec;
}

point triangle::operator()(size_t v )
{
	point p;
	p.x = *(m_vertex_list[v].x);
	p.y = *(m_vertex_list[v].y);
//	p.z = *(m_vertex_list[v].z);
	return p;
}

point triangle::get_center()
{
	return m_center;
}
