#pragma once

#include <armadillo>
#include <vector>

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
	arma::mat coord;
	std::vector<std::triangle*> triangles;
};
class bounding_rect
{
public:
	bounding_rect(matlab* m_engine);
	void make(const arma::vec* x, const arma::vec* y, int n_segments);
	arma::mat get_rect(int i);
	int n_segments;

private:
	std::vector<arma::mat> m_rectangles;
	matlab* m_engine;
};


