#ifndef SUPPORT_H
#define SUPPORT_H

#include "global.h"
#include <string>

class LogInfo {

public:

	LogInfo();
	~LogInfo();
	void print_features(Configuration &config);
	int print_this(const char* pstring, int level);
	int print_this(string pstring, int level);

	void set_printLevel(int level){m_printLevel = level;};

private:
	int m_printLevel; ///sets the level(inclusive) up to where print information should be displayed, the higher the more strict

};

#endif
