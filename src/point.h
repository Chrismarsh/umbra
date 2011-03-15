#pragma once

class point
{
public:
	point()
	{
		x=0;
		y=0;
	}

	double x;
	double y;
};

class ptr_point
{
public:
	ptr_point()
	{
		x=0;
		y=0;
	}

	double* x;
	double* y;
};

