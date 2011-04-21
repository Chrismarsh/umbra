#pragma once

#include <armadillo>
#include <vector>

#include "triangle.h"
#include "matlab_engine.h"

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
	rect::rect( arma::mat* coord )
	{
		this->coord = coord;
	}
	arma::mat* coord;
	std::vector<triangle*> triangles;
	std::vector<int> m_globalID;
};
class bounding_rect
{
public:
	bounding_rect(matlab* m_engine);
	void make(const arma::vec* x, const arma::vec* y, int n_rows, int n_cols);
	rect* bounding_rect::get_rect( int i, int j );
	int n_rows;
	int n_cols;
	arma::vec *bbx;
	arma::vec *bby;
	bool pt_in_rect(double x, double y, rect* r);

private:
	//std::vector<rect*> m_rectangles;
	std::vector<std::vector<rect*> > m_grid;
	matlab* m_engine;
};


