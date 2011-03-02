#include "matlab_engine.h"
void matlab::start()
{
	if (!(m_engine = engOpen("\0"))) 
	{
		throw std::exception("Can't start MATLAB engine");
	}
}

void matlab::stop()
{
	if(m_engine)
	{
		if (engEvalString(m_engine, "close;") == 1)
			throw std::exception("Command failed");
	}
}

void matlab::evaluate( std::string command )
{
	if(m_engine)
	{
		if (engEvalString(m_engine, command.c_str()) ==1 )
			throw std::exception("Command failed");
	}
	else
	{
		throw std::exception("No MATLAB engine open");
	}
}

void matlab::put( std::string name, mxArray* var )
{
	if(m_engine)
	{
		if( engPutVariable(m_engine, name.c_str(), var) == 1)
			throw std::exception("Command failed");
	}
	else
	{
		throw std::exception("No MATLAB engine open");
	}
}

mxArray* matlab::get( std::string name )
{

		if(m_engine)
		{
			return engGetVariable(m_engine, name.c_str());
		}
		else
		{
			throw std::exception("No MATLAB engine open");
		}
}

matlab::~matlab()
{

}

matlab::matlab()
{
	m_engine = NULL;
}
