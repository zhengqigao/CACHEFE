#ifndef __DEMO__H__
#define __DEMO__H__

#include"APA.h"

namespace demo {
	bool rundemo(int simtype, int par_nsimiter = 1000);
	bool MCdemo(int simtype);
	bool cleardemo();
}

#endif