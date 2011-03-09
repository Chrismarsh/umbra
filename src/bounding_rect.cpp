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
			arma::mat rect(5,2);
			
			rect	<< bbx(0)				 << bby(0)					<< arma::endr  //bottom left
					<< midpoint_bottom_x[0]	 << midpoint_bottom_y[0]	<< arma::endr  //bottom right mid point
					<< midpoint_top_x[0]	 << midpoint_top_y[0]		<< arma::endr  //top right mid point
					<< bbx(3)				 << bby(3)					<< arma::endr  //top right
					<< bbx(0)				 << bby(0)					<< arma::endr;

			m_rectangles.push_back(rect);
		}
		//       last rect
		else if( i==n_segments-1)
		{
			arma::mat rect(5,2);

			rect	<< 	midpoint_bottom_x[i-1]	<< 	midpoint_bottom_y[i-1]		<< arma::endr  //bottom left
					<< 	bbx(1)					<<  bby(1)						<< arma::endr  //bottom right mid point
					<<  bbx(2)					<< 	bby(2)						<< arma::endr  //top right mid point
					<< 	midpoint_top_x[i-1]		<< 	midpoint_top_y[i-1]			<< arma::endr  //top right
					<<  midpoint_bottom_x[i-1]	<< 	midpoint_bottom_y[i-1]		<< arma::endr;

			m_rectangles.push_back(rect);

		}
		else
		{

			arma::mat rect(5,2);

			rect	<< 	midpoint_bottom_x[i-1]	<< 	midpoint_bottom_y[i-1]		<< arma::endr  //bottom left
					<< 	midpoint_bottom_x[i]	<<  midpoint_bottom_y[i]		<< arma::endr  //bottom right mid point
					<<  midpoint_top_x[i]		<< 	midpoint_top_y[i]			<< arma::endr  //top right mid point
					<< 	midpoint_top_x[i-1]		<< 	midpoint_top_y[i-1]			<< arma::endr  //top right
					<<  midpoint_bottom_x[i-1]	<< 	midpoint_bottom_y[i-1]		<< arma::endr;

			m_rectangles.push_back(rect);
		}
	
	}
}

arma::mat bounding_rect::get_rect( int i )
{
	return m_rectangles[i];
}
