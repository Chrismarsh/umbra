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

void matlab::set_working_dir( std::string dir )
{
	if(m_engine)
	{
		engEvalString(m_engine, (std::string("cd '")+dir+std::string("'")).c_str());
	}
	else
	{
		throw std::exception("No MATLAB engine open");
	}
}

void matlab::set_working_dir()
{
	char path[_MAX_PATH];
	_getcwd(path, _MAX_PATH);

	set_working_dir(std::string(path));

}

std::string matlab::get_last_error()
{
	if(m_engine)
	{
		int retval = engEvalString(m_engine,"myErr=lasterror");
		// get the struct
		mxArray *err = engGetVariable(m_engine,"myErr");
		char str[512];
		if(mxIsStruct(err))
		{
			// get the error message string field
			mxArray *errStr = mxGetField(err,0,"message");
			if( (errStr != NULL) && mxIsChar(errStr) )
			{
				
				// get the string
				retval =
					mxGetString(errStr,str,sizeof(str)/sizeof(str[0]));

			}
		}
		mxDestroyArray(err);

		return std::string(str);
	}
	else
	{
		return "";
	}
}


