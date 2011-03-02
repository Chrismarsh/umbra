#pragma once


//simple vertex class for a triangle. The x,y,z values are pointers into a data source
class vertex
{
public:
	vertex& operator=(vertex const& t);

	double* x;
	double* y;
	double* z;
};

vertex& vertex::operator=( vertex const& t )
{
	//copy ptrs
	x = t.x;
	y = t.y;
	z = t.z;

}
