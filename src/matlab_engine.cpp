// % Copyright (c) 2011, Chris Marsh
// 	% All rights reserved.
// 	% 
// 	% Redistribution and use in source and binary forms, with or without 
// 	% modification, are permitted provided that the following conditions are 
// 	% met:
// % 
// 	%     * Redistributions of source code must retain the above copyright 
// 	%       notice, this list of conditions and the following disclaimer.
// 	%     * Redistributions in binary form must reproduce the above copyright 
// 	%       notice, this list of conditions and the following disclaimer in 
// 	%       the documentation and/or other materials provided with the distribution
// 	%     * Neither the name of the University of Saskatchewan nor the names 
// 	%       of its contributors may be used to endorse or promote products derived 
// 	%       from this software without specific prior written permission.
// 	%       
// 	% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// 	% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// 	% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
// 	% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
// 	% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
// 	% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
// 	% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
// 	% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
// 	% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
// 	% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
// 	% POSSIBILITY OF SUCH DAMAGE.


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


