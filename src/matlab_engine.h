#pragma once

#include <engine.h>
#include <iostream>
#include <direct.h> // for getcwd
#include <stdlib.h>// for MAX_PATH
#include <string>

class matlab
{
public:
	matlab();
	~matlab();

	void put(std::string name, mxArray* var );
	mxArray* get(std::string name);
	void evaluate(std::string command);

	//set working directory for matlab engine to a specified directory
	void set_working_dir(std::string dir);

	//set working directory for matlab engine to the current applications working directory
	void set_working_dir(); 

	std::string get_last_error();

	void start();

	void stop();

	

private:
	Engine *m_engine;

};