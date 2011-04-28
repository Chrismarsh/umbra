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
	rect( arma::mat* coord )
	{
		this->coord = coord;
	}
	~rect()
	{
		delete coord;
	}
	arma::mat* coord;
	std::vector<triangle*> triangles;
	std::vector<int> m_globalID;
};
class bounding_rect
{
public:
	bounding_rect(matlab_engine* m_engine);
	~bounding_rect();
	void make(const arma::vec* x, const arma::vec* y, int n_rows, int n_cols);
	rect* bounding_rect::get_rect( int i, int j );
	int n_rows;
	int n_cols;
	bool pt_in_rect(double x, double y, rect* r);

private:
	//std::vector<rect*> m_rectangles;
	std::vector<std::vector<rect*> > m_grid;
	matlab_engine* m_engine;
};


