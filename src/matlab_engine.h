#pragma once

#include <engine.h>
#include <iostream>
class matlab
{
public:
	matlab();
	~matlab();

	void put(std::string name, mxArray* var );
	mxArray* get(std::string name);
	void evaluate(std::string command);

	void start();

	void stop();

	

private:
	Engine *m_engine;

};