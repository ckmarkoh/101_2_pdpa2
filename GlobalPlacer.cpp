#include "GlobalPlacer.h"
#include "ExampleFunction.h"
#include "NumericalOptimizer.h"
#include "Util.h"
#include <stdlib.h>
#include <math.h>
GlobalPlacer::GlobalPlacer(Placement &placement)
	:_placement(placement)
{

}


void GlobalPlacer::place()
{
	///////////////////////////////////////////////////////////////////
	// The following example is only for analytical methods.
	// if you use other methods, you can skip and delete it directly.
	//////////////////////////////////////////////////////////////////

	ExampleFunction ef; // require to define the object function and gradient function
//	double avg_x = 0;
//	double avg_y = 0;
/*	for (unsigned i=0; i<_placement.numModules();i++){
		cout<<i <<" "<<_placement.module(i).name()<<" "<<(int)_placement.module(i).centerX()<<" "<<(int)_placement.module(i).centerY()<<endl;

		
	}
	for (unsigned i=0; i<_placement.numNets();i++){
		cout<<i <<" "<<_placement.net(i).numPins()<<endl;
	}*/

	testNetWA();

    vector<double> x(2);//(_placement.numModules()*2); // solution vector, size: num_blocks*2 
                         // each 2 variables represent the X and Y dimensions of a block
    x[0] = 50; // initialize the solution vector
    x[1] = 50;
	


    NumericalOptimizer no(ef);
    no.setX(x); // set initial solution
    no.setNumIteration(5); // user-specified parameter
    no.setStepSizeBound(5); // user-specified parameter
    no.solve(); // Conjugate Gradient solver

    cout << "Current solution:" << endl;
    for (unsigned i = 0; i < no.dimension(); i++) {
        cout << "x[" << i << "] = " << no.x(i) << endl;
    }
	//cout<<"1"<<endl;
    cout << "Objective: " << no.objective() << endl;
	////////////////////////////////////////////////////////////////
}
void GlobalPlacer::testNetWA(){
/*	Pin* p1= new Pin();
	Pin* p2= new Pin();
	Pin* p3= new Pin();
	Pin* p4= new Pin();

	p1->setPosition(10,3);
	p2->setPosition(3,9);
	p3->setPosition(1,4);
	p4->setPosition(6,2);

	Net n;
	n.addPin(p1);
	n.addPin(p2);
	n.addPin(p3);
	n.addPin(p4);
	
	cout<<"WA:"<<getNetWA(n)<<endl;*/
	for( double i=0;i<100;i++){
		cout<<sigMoid(10,i,90)<<endl;
	}
}
double GlobalPlacer::getNetWA(Net& n ){
	double r=10;
	double xd0=0;
	double xd1=0;
	double yd0=0;
	double yd1=0;
	double xn0=0;
	double xn1=0;
	double yn0=0;
	double yn1=0;
	for(unsigned i=0;i<n.numPins();i++){
		Pin & p=n.pin(i);
		double x=p.x();
		double y=p.y();
		double exp_1x=pow(E,x/r);
		double exp_0x=pow(E,(-x)/r);
		double exp_1y=pow(E,y/r);
		double exp_0y=pow(E,(-y)/r);
		xd1+=(x*exp_1x);
		xn1+=exp_1x;

		xd0+=(x*exp_0x);
		xn0+=exp_0x;

		yd1+=(y*exp_1y);
		yn1+=exp_1y;

		yd0+=(y*exp_0y);
		yn0+=exp_0y;
	
	}
	return xd1/xn1-xd0/xn0+yd1/yn1-yd0/yn0;
}



