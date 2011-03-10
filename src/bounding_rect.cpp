#include "bounding_rect.h"



void bounding_rect::make(const arma::vec* x, const arma::vec* y, int n_segments)
{
	arma::vec midpoint_bottom_x(n_segments);
	arma::vec midpoint_top_x(n_segments);

	arma::vec midpoint_bottom_y(n_segments);
	arma::vec midpoint_top_y(n_segments);

	
	m_engine->evaluate("[bbx,bby,~,~]=minboundrect(mxRot(:,1),mxRot(:,2));");
	arma::vec& bbx = *(m_engine->get_double_vector("bbx"));
	arma::vec& bby = *(m_engine->get_double_vector("bby"));
	m_engine->evaluate("clear bbx bby");

	this->n_segments = n_segments;

	this->bbx = new arma::vec(bbx);
	this->bby = new arma::vec(bby);


	double m = (bby(1)-bby(0))/(bbx(1)-bbx(0));

	double step=(bbx(2)-bbx(3))/n_segments;
	for (int i = 0; i<n_segments;i++)
	{
		double xpos_bottom = bbx(0)+step*(i+1);
		double xpos_top     =bbx(3)+step*(i+1);
		         
		midpoint_bottom_x[i]=xpos_bottom;
		midpoint_bottom_y[i]=m*(xpos_bottom-bbx(0))+bby(0); //2pt line eqn
			       
		midpoint_top_x[i]=xpos_top;
		midpoint_top_y[i]=m*(xpos_top-bbx(3))+bby(3); //2pt line eqn

	}

	m_rectangles.reserve(n_segments);

	for (int i = 0; i<n_segments;i++)
	{
		if (i == 0)
		{
			arma::mat* c = new arma::mat(5,2);
			arma::mat& coord = *c;
			coord	<< bbx(0)				 << bby(0)					<< arma::endr  //bottom left
					<< midpoint_bottom_x[0]	 << midpoint_bottom_y[0]	<< arma::endr  //bottom right mid point
					<< midpoint_top_x[0]	 << midpoint_top_y[0]		<< arma::endr  //top right mid point
					<< bbx(3)				 << bby(3)					<< arma::endr  //top right
					<< bbx(0)				 << bby(0)					<< arma::endr;
			
			m_rectangles.push_back(new rect(c));
		}
		//       last rect
		else if( i==n_segments-1)
		{
			arma::mat* c = new arma::mat(5,2);
			arma::mat& coord = *c;

			coord	<< 	midpoint_bottom_x[i-1]	<< 	midpoint_bottom_y[i-1]		<< arma::endr  //bottom left
					<< 	bbx(1)					<<  bby(1)						<< arma::endr  //bottom right mid point
					<<  bbx(2)					<< 	bby(2)						<< arma::endr  //top right mid point
					<< 	midpoint_top_x[i-1]		<< 	midpoint_top_y[i-1]			<< arma::endr  //top right
					<<  midpoint_bottom_x[i-1]	<< 	midpoint_bottom_y[i-1]		<< arma::endr;

			m_rectangles.push_back(new rect(c));

		}
		else
		{

			arma::mat* c = new arma::mat(5,2);
			arma::mat& coord = *c;

			coord	<< 	midpoint_bottom_x[i-1]	<< 	midpoint_bottom_y[i-1]		<< arma::endr  //bottom left
					<< 	midpoint_bottom_x[i]	<<  midpoint_bottom_y[i]		<< arma::endr  //bottom right mid point
					<<  midpoint_top_x[i]		<< 	midpoint_top_y[i]			<< arma::endr  //top right mid point
					<< 	midpoint_top_x[i-1]		<< 	midpoint_top_y[i-1]			<< arma::endr  //top right
					<<  midpoint_bottom_x[i-1]	<< 	midpoint_bottom_y[i-1]		<< arma::endr;

			m_rectangles.push_back(new rect(c));
		}
	
	}
}

rect* bounding_rect::get_rect( int i )
{
	 return m_rectangles.at(i);
}

bounding_rect::bounding_rect( matlab* m_engine )
{
	this->m_engine = m_engine;
}

bool bounding_rect::pt_in_rect( double x, double y, rect* r )
{
	for (int i = 0; i<4; i++)
	{
		double x0 = r->coord->operator()(i,0);
		double y0 = r->coord->operator()(i,1);
		double x1 = r->coord->operator()(i+1,0);
		double y1 = r->coord->operator()(i+1,1);

		double pt = ((y - y0)*(x1 - x0) - (x - x0)*(y1 - y0));
		if (pt < 0)
			return false;
		else if (pt ==0)
			return false;
	}
	return true;
}