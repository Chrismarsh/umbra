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

triangle::triangle( point vertex1, point vertex2, point vertex3, size_t cur_rec_depth/*=1*/)
{
	m_cur_rec_depth = cur_rec_depth;
 	m_sub_tri = NULL;
	set_vertex_values(vertex1,vertex2, vertex3);

	radiation = 0.0;
	shadow = 0.0;
	z_prime = 0.0;

}

triangle::triangle(size_t cur_rec_depth)
{
	m_sub_tri = NULL;

	m_cur_rec_depth = cur_rec_depth;

	radiation = 0.0;
	shadow = 0.0;
	z_prime = 0.0;
}


bool triangle::contains( point xy )
{
	return contains(xy.x,xy.y);
}

void triangle::set_vertex_values( point vertex1, point vertex2, point vertex3)
{
	m_vertex_list[0].x = vertex1.x;
	m_vertex_list[0].y = vertex1.y;
	m_vertex_list[0].z = vertex1.z;

	m_vertex_list[1].x = vertex2.x;
	m_vertex_list[1].y = vertex2.y;
	m_vertex_list[1].z = vertex2.z;

	m_vertex_list[2].x = vertex3.x;
	m_vertex_list[2].y = vertex3.y;
	m_vertex_list[2].z = vertex3.z;

	//center of the triangle
	point pos = calc_center(this);		

	m_center.x = pos.x;
	m_center.y = pos.y;
// 	std::cout.precision(16);
// 	std::cout << "vertex" <<std::endl;
// 	std::cout << m_vertex_list[0].x <<" " << m_vertex_list[0].y << std::endl
// 			  << m_vertex_list[1].x <<" " << m_vertex_list[1].y << std::endl
// 			  << m_vertex_list[2].x <<" " << m_vertex_list[2].y << std::endl;
// 	std::cout << "center" <<std::endl;
// 	std::cout << pos.x << " " << pos.y << std::endl;

	if (m_cur_rec_depth != 0)
	{
		update_subtri();
	}
	
}

void triangle::update_subtri()
{
	int longest;
 	if (m_sub_tri)
	{
		for(int i = 0;i<4;i++)
		{
			delete m_sub_tri[i];
		}
 		delete[] m_sub_tri;	
	}

	// set each sub
	m_sub_tri = new triangle*[4];
	for(int i = 0; i<4;i++)
	{
		m_sub_tri[i] = new triangle(m_cur_rec_depth-1);
	}

// 	std::cout.precision(16);
// 	std::cout << "Main triangle" << std::endl;
// 	std::cout << "[" << m_vertex_list[0].x << " " << m_vertex_list[0].y << "];" << std::endl
// 				<< "[" << m_vertex_list[1].x << " " << m_vertex_list[1].y << "];" << std::endl
// 				<< "[" << m_vertex_list[2].x << " " << m_vertex_list[2].y << "];" << std::endl;

	arma::vec l_sides(3);
	//0-1
	l_sides[0] = sqrt( ((m_vertex_list[0].x)-(m_vertex_list[1].x))*((m_vertex_list[0].x)-(m_vertex_list[1].x)) + (((m_vertex_list[0].y)-(m_vertex_list[1].y))*((m_vertex_list[0].y)-(m_vertex_list[1].y))) * 1.0);
	//1-2
	l_sides[1] = sqrt( ((m_vertex_list[1].x)-(m_vertex_list[2].x))*((m_vertex_list[1].x)-(m_vertex_list[2].x)) + (((m_vertex_list[1].y)-(m_vertex_list[2].y))*((m_vertex_list[1].y)-(m_vertex_list[2].y)))* 1.0);
	//0-2
	l_sides[2] = sqrt( ((m_vertex_list[0].x)-(m_vertex_list[2].x))*((m_vertex_list[0].x)-(m_vertex_list[2].x)) + (((m_vertex_list[0].y)-(m_vertex_list[2].y))*((m_vertex_list[0].y)-(m_vertex_list[2].y)))* 1.0);

	if(l_sides[0] > l_sides[1] || l_sides[0] > l_sides[2])
	{
		longest = 0;
		//midpoint of the longest edge (which is 0-1)
		point *midpt_01 = midpoint(m_vertex_list[0],m_vertex_list[1]);
		
		//midpoint of the opposite edge (which is 0-2)
		point* midpt_02 =midpoint(m_vertex_list[0],m_vertex_list[2]);


		//midpoint of the opposite edge2 (which is 1-2)
		point* midpt_12 = midpoint(m_vertex_list[1],m_vertex_list[2]);
	
		//have 4 sub triangles now:
		//define CCW, left->right
		//t1 = 0-01-02-01
		//t2 = 01-2-02-01
		//t3 = 01-12-2-01
		//t4 = 01-1-12-01

		//t1 = 0-01-02-01
		point v1;
		v1.x = m_vertex_list[0].x;
		v1.y = m_vertex_list[0].y;
		point v2;
		v2.x = (midpt_01->x);
		v2.y = (midpt_01->y);	
		point v3;
		v3.x = (midpt_02->x);
		v3.y = (midpt_02->y);
		m_sub_tri[0]->set_vertex_values(v1,v2,v3);

		//t2 = 01-2-02-01
		//ptr_point v1;
		v1.x = (midpt_01->x);
		v1.y = (midpt_01->y);
		//ptr_point v2;
		v2.x = m_vertex_list[2].x;
		v2.y = m_vertex_list[2].y;
		//ptr_point v3;
		v3.x = (midpt_02->x);
		v3.y = (midpt_02->y);
		m_sub_tri[1]->set_vertex_values(v1,v2,v3);

		//t3 = 01-12-2-01
		//ptr_point v1;
		v1.x = (midpt_01->x);
		v1.y = (midpt_01->y);
		//ptr_point v2;
		v2.x = (midpt_12->x);
		v2.y = (midpt_12->y);
		//ptr_point v3;
		v3.x = m_vertex_list[2].x;
		v3.y = m_vertex_list[2].y;
		m_sub_tri[2]->set_vertex_values(v1,v2,v3);

		//t4 = 01-1-12-01
		//v1;
		v1.x = (midpt_01->x);
		v1.y = (midpt_01->y);
		//ptr_point v2;
		v2.x = m_vertex_list[1].x;
		v2.y = m_vertex_list[1].y;
		//ptr_point v3;
		v3.x = (midpt_12->x);
		v3.y = (midpt_12->y);
		m_sub_tri[3]->set_vertex_values(v1,v2,v3);

		

	}
	else if(l_sides[1] > l_sides[2] || l_sides[1] > l_sides[0])
	{
		longest = 1;

		//midpoint of the longest edge (which is 1-2)
		point* midpt_12 = midpoint(m_vertex_list[1],m_vertex_list[2]);


		//midpoint of the opposite edge (which is 0-1)
		point* midpt_01 = midpoint(m_vertex_list[0],m_vertex_list[1]);

		//midpoint of the opposite edge2 (which is 0-2)
		point* midpt_02 = midpoint(m_vertex_list[0],m_vertex_list[2]);

		//have 4 sub triangles now:
		//define CCW, left->right
		//t1 = 0-01-12-0
		//t2 = 01-1-12-01
		//t3 = 12-2-02-12
		//t4 = 0-12-02-0

		//t1 = 0-01-12-0
		point v1;
		v1.x = m_vertex_list[0].x;
		v1.y = m_vertex_list[0].y;
		point v2;
		v2.x = (midpt_01->x);
		v2.y = (midpt_01->y);	
		point v3;
		v3.x = (midpt_12->x);
		v3.y = (midpt_12->y);
		m_sub_tri[0]->set_vertex_values(v1,v2,v3);

		//t2 = 01-1-12-01

		v1.x = (midpt_01->x);
		v1.y = (midpt_01->y);

		v2.x = m_vertex_list[1].x;
		v2.y = m_vertex_list[1].y;	

		v3.x = (midpt_12->x);
		v3.y = (midpt_12->y);
		m_sub_tri[1]->set_vertex_values(v1,v2,v3);

		//t3 = 12-2-02-12

		v1.x = (midpt_12->x);
		v1.y = (midpt_12->y);

		v2.x = m_vertex_list[2].x;
		v2.y = m_vertex_list[2].y;	

		v3.x = (midpt_02->x);
		v3.y = (midpt_02->y);
		m_sub_tri[2]->set_vertex_values(v1,v2,v3);

		//t4 = 0-12-02-0

		v1.x = m_vertex_list[0].x;
		v1.y = m_vertex_list[0].y;

		v2.x = (midpt_12->x);
		v2.y = (midpt_12->y);	

		v3.x = (midpt_02->x);
		v3.y = (midpt_02->y);
		m_sub_tri[3]->set_vertex_values(v1,v2,v3);


	}
	else 
	{
		longest = 2;

		//midpoint of the longest edge (which is 0-2)
		point* midpt_02 = midpoint(m_vertex_list[0],m_vertex_list[2]);

		//midpoint of the opposite edge (which is 0-1)
		point* midpt_01 = midpoint(m_vertex_list[0],m_vertex_list[1]);

		//midpoint of the opposite edge2 (which is 1-2)
		point* midpt_12 = midpoint(m_vertex_list[1],m_vertex_list[2]);

		//have 4 sub triangles now:
		//define CCW, left->right
		//t1 = 0-01-02
		//t2 = 01 -1 -02
		//t3 = 1 - 12-02
		//t4 = 12 - 2 -02

		//t1 =0-01-02
		point v1;
		v1.x = m_vertex_list[0].x;
		v1.y = m_vertex_list[0].y;
		point v2;
		v2.x = (midpt_01->x);
		v2.y = (midpt_01->y);	
		point v3;
		v3.x = (midpt_02->x);
		v3.y = (midpt_02->y);
		m_sub_tri[0]->set_vertex_values(v1,v2,v3);

		//t2 = 01 -1 -02

		v1.x = (midpt_01->x);
		v1.y = (midpt_01->y);

		v2.x = m_vertex_list[1].x;
		v2.y = m_vertex_list[1].y;

		v3.x = (midpt_02->x);
		v3.y = (midpt_02->y);
		m_sub_tri[1]->set_vertex_values(v1,v2,v3);

		//t3 = 1 - 12-02

		v1.x = m_vertex_list[1].x;
		v1.y = m_vertex_list[1].y;

		v2.x = (midpt_12->x);
		v2.y = (midpt_12->y);

		v3.x = (midpt_02->x);
		v3.y = (midpt_02->y);
		m_sub_tri[2]->set_vertex_values(v1,v2,v3);

		//t4 = 12 - 2 -02
		v1.x = (midpt_12->x);
		v1.y = (midpt_12->y);

		v2.x = m_vertex_list[2].x;
		v2.y = m_vertex_list[2].y;

		v3.x = (midpt_02->x);
		v3.y = (midpt_02->y);
		m_sub_tri[3]->set_vertex_values(v1,v2,v3);

	}	
}
triangle& triangle::sub_tri(size_t t)
{
	return *(m_sub_tri[t]);
}


int triangle::intersects( triangle* t )
{
	//I have no children
	if(m_sub_tri == NULL)
	{
		//does t contain my points?
		bool intersect = t->contains(this->get_center());

		if(!intersect)
		{
			if( t->contains(this->get_vertex(0)) ||
				t->contains(this->get_vertex(1)) ||
				t->contains(this->get_vertex(2)))
			{
				intersect = true;
			}
		}

		return intersect == true ? 1:0;
	}
	else
		//i have children
	{
		int sum=0;
		for(size_t i = 0; i<4; i++)
		{
			sum += m_sub_tri[i]->intersects(t);
		}
		return sum;
	}
}

bool triangle::contains(double x, double y)
{
	double x1=m_vertex_list[0].x;
	double y1=m_vertex_list[1].y;

	double x2=m_vertex_list[1].x;
	double y2=m_vertex_list[1].y;

	double x3=m_vertex_list[2].x;
	double y3=m_vertex_list[2].y;
	double lambda1= ((y2-y3)*(x-x3)+(x3-x2)*(y-y3))/((y2-y3)*(x1-x3)+(x3-x2)*(y1-y3));
	double lambda2= ((y3-y1)*(x-x3)+(x1-x3)*(y-y3))/((y3-y1)*(x2-x3)+(x1-x3)*(y2-y3));
	double lambda3=1-lambda1-lambda2;

	return lambda1 > 0 && lambda1 <1 
		&& lambda2 > 0 && lambda2 < 1 
		&& lambda3 > 0 && lambda3 <1;


}

point triangle::get_vertex( size_t vertex )
{
	return m_vertex_list[vertex];
}
// 
// point triangle::operator()(size_t v )
// {
// 	return m_vertex_list[vertex];
// }


point triangle::get_center()
{
	return m_center;
}

point* triangle::midpoint( point& p1, point& p2 )
{
	point* mp = new point;
	//new x value 1/2 way there
	mp->x = p1.x + (p2.x - p1.x)/2.0;
	//slope
	double m = (p2.y-p1.y)/(p2.x-p1.x);
	double b = p1.y-m*p1.x;
	mp->y = m*mp->x+b;

	return mp;
}

point triangle::calc_center( triangle* t )
{

	arma::vec Pa(3);
	Pa(0) = (t->m_vertex_list[0].x);
	Pa(1) = (t->m_vertex_list[0].y);
	Pa(2) = 0;

	arma::vec Pb(3);
	Pb(0) = (t->m_vertex_list[1].x);
	Pb(1) = (t->m_vertex_list[1].y);
	Pb(2) = 0;

	arma::vec Pc(3);
	Pc(0) = (t->m_vertex_list[2].x);
	Pc(1) = (t->m_vertex_list[2].y);
	Pc(2) = 0;

	arma::vec AB = Pb - Pa;
	arma::vec AC = Pc - Pa;
	arma::vec BC = Pc - Pb;

//circumcenter
// 	arma::vec N = arma::cross(AC,AB);
// 	arma::vec L1 = arma::cross(AB,N);
// 	arma::vec L2 = arma::cross(BC,N);
// 	arma::vec P21 = (Pc - Pa)/2;
// 	arma::vec P1 = (Pa + Pb)/2;

	
//incenter

	arma::vec uab  = AB / arma::norm(AB,1);
	arma::vec uac  = AC / arma::norm(AC,1);
	arma::vec ubc  = BC / arma::norm(BC,1);
	arma::vec uba = -uab;

	arma::vec L1 = uab + uac;
	arma::vec L2 = uba + ubc;
	arma::vec P21 = Pb-Pa;
	arma::vec P1 = Pa;
     
	arma::mat ML1(L1);
	arma::mat ML2(L2);

	arma::mat ML = arma::join_rows(ML1,-ML2);

	arma::vec lambda = arma::solve(ML,P21);

	arma::vec pos = P1+lambda(0)*L1;

	point p;
	p.x = pos(0);
	p.y = pos(1);
	return p;
}

void triangle::set_facenormal( arma::vec& normal )
{
	m_surface_normal = normal;
}

arma::vec triangle::get_facenormal()
{
	return m_surface_normal;
}

void triangle::compute_azimuth()
{
	//convert normal to spherical
	double r = sqrt(m_surface_normal(0)*m_surface_normal(0) + m_surface_normal(1)*m_surface_normal(1) + m_surface_normal(2)*m_surface_normal(2));
	double theta = acos(m_surface_normal(2)/r); 

	//y=north
	double phi = atan2(m_surface_normal(1),m_surface_normal(0)); // + M_PI /*-3.14159/2*/; //south == 0
	m_azimuth = phi - M_PI/2; //set north = 0

	if(m_azimuth < 0.0)
		m_azimuth += 2*M_PI;
}

void triangle::compute_slope()
{
	//z surface normal
	arma::vec n(3);
	n(0) = 0.0; //x
	n(1) = 0.0; //y
	n(2) = 1.0;

	m_slope = acos( arma::norm_dot(m_surface_normal,n));
}

double triangle::azimuth()
{
	return m_azimuth;
}

double triangle::slope()
{
	return m_slope;
}