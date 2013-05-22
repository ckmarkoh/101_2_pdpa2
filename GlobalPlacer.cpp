#include "GlobalPlacer.h"
#include "ExampleFunction.h"
#include "NumericalOptimizer.h"
#include "Util.h"
#include <stdlib.h>
#include <math.h>
#include <algorithm>

#include <math.h>
#include<iomanip>

#define R 0.9
#define T 1000
#define FROZEN 0.001
#define N 500

using namespace std;
double Bin::_width=0;
double Bin::_height=0;

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
//	ExampleFunction ef; // require to define the object function and gradient function

//	cout<<"total_x:"<<total_x<<endl;
//	double row_n=total_x/( _placement.boundryRight()-_placement.boundryLeft()  );
//	cout<<"row_n:"<<row_n<<endl;
//	cout<<"require_height:"<<max_y*row_n<<endl;
/*	for (unsigned i=0; i<_placement.numNets();i++){
		cout<<i <<" "<<_placement.net(i).numPins()<<endl;
	}*/

//	testNetWA();

  vector<double> x(1);//(_placement.numModules()*2); // solution vector, size: num_blocks*2 
                         // each 2 variables represent the X and Y dimensions of a block
    x[0] = 50; // initialize the solution vector
   // x[1] = 50;
	/*	
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
	*/
//	unsigned ax=11;
//	unsigned ay=6;
//	cout<<ax/ay<<ax%ay<<endl;
	initialize();
	printf( "\nHPWL: %.0f\n",_placement.computeHpwl() );	



	simu_anneal_sub(R,T,FROZEN,N);
//	testBin();
}

void GlobalPlacer::initialize(){
	//double total_x = 0;
	double max_x = 0;
	double max_y = 0;
	for (unsigned i=0; i<_placement.numModules();i++){
		double module_x=_placement.module(i).width();
		double module_y=_placement.module(i).height();
//		cout<<i <<" "<<_placement.module(i).name()<<" "<<(int)module_x<<" "<<(int)module_y<<endl;
//		total_x+=module_x;
		max_y=max(max_y,module_y);
		max_x=max(max_x,module_x);
	}
	cout<<"max_x:"<<max_x<<endl;
	cout<<"max_y:"<<max_y<<endl;
	Bin::setWidthHeight(max_x,max_y);
	reset_bin(unsigned( (_placement.boundryRight()-_placement.boundryLeft())/max_x) ,
				unsigned((_placement.boundryTop()-_placement.boundryBottom())/max_y) );
	cout<<"binx:"<<_bin_numx<<" biny:"<<_bin_numy <<endl;
	cout<<"***********0*************"<<endl;
	unsigned i=0;
	while(i<_placement.numModules()){
		unsigned j=0;
		while(j<_bin_numx){
			unsigned k=0;
			while(k<_bin_numy){
				if(_bins_vec[j][k]->capacity()>0){
				//	cout<<"j:"<<j<<" k:"<<k<<" i:"<<i<<endl;
				//	cout<<_bins_vec[j][k]->size()<<endl;
				//	cout<<_placement.module(i).name()<<endl;
					_bins_vec[j][k]->addModule(&_placement.module(i++));
					if(i>=_placement.numModules()){break;}
				}
				else{k++;}
			}
			if(i>=_placement.numModules()){break;}
			j++;
		}
		if(i>=_placement.numModules()){break;}
	}
	cout<<"***********1*************"<<endl;

}

void GlobalPlacer::printAllBins(int zi){
	
	cout<<"***********start_"<<zi<<"*************"<<endl;
	for(unsigned j=0;j<_bin_numx;j++){
		for(unsigned k=0;k<_bin_numy;k++){
			for(unsigned i=0;i<_bins_vec[j][k]->size();i++ ){
			//	cout<<"j:"<<j<<" k:"<<k<<" i:"<<i<<endl;
			//	cout<<_bins_vec[j][k]->size()<<endl;
				cout<<_bins_vec[j][k]->getModule(i)->name()<<endl;
			}
		}
	}
	cout<<"***********end_"<<zi<<"*************"<<endl;
}

void GlobalPlacer::testBin(){
	printAllBins(0);
	random_neighbor();
	printAllBins(1);
	restore_backup();
	printAllBins(2);
}


double GlobalPlacer::getNetEU(Net& n ){
	double r=1;
	double eu=0;
	assert(n.numPins()>=2);
	for(unsigned i=0; i<n.numPins(); i++){
		Pin & p1=n.pin(i);
		for(unsigned j=i+1; j<n.numPins(); j++){
			Pin & p2=n.pin(j);
			eu+=pow(p1.x()-p2.x(),2)+pow(p1.y()-p2.y(),2);
		}	
	}
	return r*eu;
}
double GlobalPlacer::getNetHPWL(Net& n ){
	Pin & p=n.pin(0);
	double min_x=p.x();
	double min_y=p.y();
	double max_x=p.x();
	double max_y=p.y();
	for(unsigned i=1;i<n.numPins();i++){
		p=n.pin(i);
		min_x=min(p.x(),min_x);
		min_y=min(p.y(),min_y);
		max_x=max(p.x(),max_x);
		max_y=max(p.y(),max_y);
	}
	return (max_x-min_x)+(max_y-min_y);
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

void GlobalPlacer::testNetWA(){
	Pin* p1= new Pin();
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
	
	cout<<"WA:"<<getNetEU(n)<<endl;
//	for( double i=0;i<100;i++){
//		cout<<sigMoid(10,i,90)<<endl;
//	}
}
bool GlobalPlacer::simu_anneal_sub(double r,double t,double frozen,double n){
	unsigned times=0;
	bool complete=false;
	unsigned fall_time=0;
	unsigned restore_time=0;
	unsigned climb_time=0;
	double tempture=t;
	while((tempture*=r)>frozen){
		for(unsigned nz=0;nz<n;nz++){
	//		cout<<times++<<endl;
			times++;
			double origin_cost=getCost();
	//		cout<<"origin_cost:"<<origin_cost<<endl;
	//		cout<<"before random_neighbor"<<endl;
			random_neighbor();
	//		cout<<"after random_neighbor"<<endl;
			double after_cost=getCost();
	//		cout<<"after_cost:"<<after_cost<<endl;
		/*	if(after_cost==-1){
				complete=true;
		//		cout<<"complete!!!"<<endl;
				break;
			}*/

			double delta_cost=after_cost-origin_cost;
	//		cout<<"delta_cost:"<<delta_cost<<endl;
			if(delta_cost>0){
				if(_rnGen( pow(E,(delta_cost/tempture)) ) >= 1.0){
			//		cout<<"  restore"<<endl;
					restore_time++;
					restore_backup();
	//				cout<<"restore"<<endl;
				}
				else
				{
		//			confirm_backup();
					climb_time++;
	//				cout<<"climb"<<endl;
				}
			}
			else{
		//		confirm_backup();
				fall_time++;
	//			cout<<"fall"<<endl;
			}
		}
	}
//	cout<<"stage_complete"<<endl;
	cout<<"times:"<<times<<endl;
	cout<<"fall_times:"<<fall_time<<endl;
	cout<<"climb_times:"<<climb_time<<endl;
	cout<<"restore_times:"<<restore_time<<endl;
//	cout<<"========report=========="<<endl;
//	_myusage.report(1,1);
//	cout<<"**********************end sa:************************8"<<endl;
//	sa_total_time+=times;
	return complete;
}
double GlobalPlacer::getCost(){
		double total_cost=0;
		for (unsigned i=0; i<_placement.numNets();i++){
			total_cost+=(getNetEU(_placement.net(i))*0.0001);
		}
		return _placement.computeHpwl()/100; 
	}
void GlobalPlacer::random_neighbor(){
	store_backup();
	unsigned i_bin1=_rnGen(_bin_numx*_bin_numy);
	unsigned i_bin2=_rnGen(_bin_numx*_bin_numy);

	while(_bins_vec[i_bin1/_bin_numy][i_bin1%_bin_numy]->size()==0){
		i_bin1=_rnGen(_bin_numx*_bin_numy);
	}
	while(i_bin1==i_bin2){
		i_bin2=_rnGen(_bin_numx*_bin_numy);
	}
	assert(i_bin1!=i_bin2);
	unsigned idx1=i_bin1/_bin_numy;
	unsigned idy1=i_bin1%_bin_numy;
	unsigned idx2=i_bin2/_bin_numy;
	unsigned idy2=i_bin2%_bin_numy;
	double cap1=_bins_vec[idx1][idy1]->capacity();
	double cap2=_bins_vec[idx2][idy2]->capacity();
	unsigned sz1=_bins_vec[idx1][idy1]->size();
	unsigned sz2=_bins_vec[idx2][idy2]->size();
	enum Move_state{
			MOVE_M1=0,
			MOVE_M2,
			MOVE_EXCHANGE
		};
	Move_state move_state;
	if((cap1>0)&&(sz2==0)){
		move_state=MOVE_M1;
		}
	else if((cap1>0)&&(cap2>0)){
		if(_rnGen(3)>2){
			move_state=MOVE_M2;
		}
		else if(_rnGen(3)<1){
			move_state=MOVE_M1;
		}
		else{
			move_state=MOVE_EXCHANGE;
		}
	}
	else if((cap1>0)&&(cap2<=0)){
		move_state=MOVE_M2;
	}
	else if((cap1<=0)&&(cap2>0)){
		move_state=MOVE_M1;
	}
	else{ //((cap1<=0)&&(cap2<=0)){
		move_state=MOVE_EXCHANGE;
	}
	switch(move_state){
		//		//cout<< "rstate2:" <<rstate<<endl;
		case MOVE_M1:{
			Module* m1=_bins_vec[idx1][idy1]->removeModule(_rnGen(sz1));
			_bins_vec[idx2][idy2]->addModule(m1);

		}break;
		case MOVE_M2:{
			Module* m2=_bins_vec[idx2][idy2]->removeModule(_rnGen(sz2));
			_bins_vec[idx1][idy1]->addModule(m2);

		}break;
		case MOVE_EXCHANGE:{
			Module* m1=_bins_vec[idx1][idy1]->removeModule(_rnGen(sz1));
			Module* m2=_bins_vec[idx2][idy2]->removeModule(_rnGen(sz2));
			_bins_vec[idx2][idy2]->addModule(m1);
			_bins_vec[idx1][idy1]->addModule(m2);

		}break;
		default:{
			assert(0);
		}break;
	}
}
void GlobalPlacer::store_backup(){
	_backup_module.clear();
	for(unsigned j=0;j<_bin_numx;j++){
		for(unsigned k=0;k<_bin_numy;k++){
			for(unsigned i=0;i<_bins_vec[j][k]->size();i++ ){
			//	cout<<"j:"<<j<<" k:"<<k<<" i:"<<i<<endl;
			//	cout<<_bins_vec[j][k]->size()<<endl;
				_backup_module.push_back(_bins_vec[j][k]->getModule(i));
			}
		}
	}
}
void GlobalPlacer::restore_backup(){

	for(unsigned j=0;j<_bin_numx;j++){
		for(unsigned k=0;k<_bin_numy;k++){
			//	cout<<"j:"<<j<<" k:"<<k<<" i:"<<i<<endl;
			//	cout<<_bins_vec[j][k]->size()<<endl;
				_bins_vec[j][k]->clear();
		}
	}


	unsigned i=0;
	while(i<_backup_module.size()){
		unsigned j=0;
		while(j<_bin_numx){
			unsigned k=0;
			while(k<_bin_numy){
				if(_bins_vec[j][k]->capacity()>0){
				//	cout<<"j:"<<j<<" k:"<<k<<" i:"<<i<<endl;
				//	cout<<_bins_vec[j][k]->size()<<endl;
				//	cout<<_backup_module.module(i).name()<<endl;
					_bins_vec[j][k]->addModule( _backup_module[i++]);
					if(i>=_backup_module.size()){break;}
				}
				else{k++;}
			}
			if(i>=_backup_module.size()){break;}
			j++;
		}
		if(i>=_backup_module.size()){break;}
	}
}
void GlobalPlacer::confirm_backup(){

}
