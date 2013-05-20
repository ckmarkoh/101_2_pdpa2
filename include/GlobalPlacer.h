#ifndef GLOBALPLACER_H
#define GLOBALPLACER_H

#include "Placement.h"

class GlobalPlacer 
{
public:
    GlobalPlacer(Placement &placement);
	void place();

	double getNetWA(Net& n );
	void testNetWA();

	//double sigMoid(Module& m);

private:
    Placement& _placement;
	
};

#endif // GLOBALPLACER_H
