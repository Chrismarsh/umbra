#pragma once

class point
{
public:
	point()
	{
		x=0;
		y=0;
		z=0;
	}

	double x;
	double y;
	double z;
};

class ptr_point
{
public:
	ptr_point()
	{
		x=0;
		y=0;
		z=0;
	}

	double* x;
	double* y;
	double* z;
};