#include "GlobalPlacer.h"
#include "ExampleFunction.h"
#include "NumericalOptimizer.h"
#include "Util.h"
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <math.h>
#include<iomanip>


#define R 0.8
#define T 1000
#define FROZEN 100

//#define R 0.9
//#define T 1000
//#define FROZEN 0.001
#define N 200

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
//printAllModules(0);
	
	initialize();
	printf( "\nHPWL: %.0f\n",_placement.computeHpwl() );	
//printAllBins(1);

//printAllModules(0);
	simu_anneal();
//	simu_anneal_sub(R,T,FROZEN,N);
//	testBin();
//printAllBins(1);
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
//	cout<<"max_x:"<<max_x<<endl;
//	cout<<"max_y:"<<max_y<<endl;
	Bin::setWidthHeight(max_x,max_y);
	reset_bin(unsigned( (_placement.boundryRight()-_placement.boundryLeft())/max_x) ,
				unsigned((_placement.boundryTop()-_placement.boundryBottom())/max_y) );
//	cout<<"binx:"<<_bin_numx<<" biny:"<<_bin_numy <<endl;
//	cout<<"***********0*************"<<endl;
//	prePlaceInitialize();
	randomInitialize();
		_initial_cost=_placement.computeHpwl();
}
void GlobalPlacer::prePlaceInitialize(){
	unsigned half_x=(unsigned)(_bin_numx)/2;
	unsigned half_y=(unsigned)(_bin_numy)/2;
	unsigned min_half=( int( min(half_x,half_y))  > 0 ) ? min(half_x,half_y) : 0;
	
	cout<<"half_x:"<<half_x<<endl;
	cout<<"half_y:"<<half_y<<endl;
	cout<<"min_half:"<<min_half<<endl;	


	vector<Module*> vec_module;
	vector<Net*> sorted_net;
	for(unsigned i=0; i<_placement.numModules(); i++){
		vec_module.push_back(& _placement.module(i));
	}
	for(unsigned i=0; i<_placement.numNets();i++){
		sorted_net.push_back(& _placement.net(i));
	}
	for(int i=0;i<(int)sorted_net.size()-1;i++){
		for(int j=(int)sorted_net.size()-1;j>=i+1;j--){
			if(sorted_net[j]->numPins() > sorted_net[j-1]->numPins()){
				swap(sorted_net[j],sorted_net[j-1]);
			}
		}
	}

	for(int i=0;i<(int)vec_module.size()-1;i++){
		for(int j=(int)vec_module.size()-1;j>=i+1;j--){
			if(vec_module[j]->numPins() > vec_module[j-1]->numPins()){
				swap(vec_module[j],vec_module[j-1]);
			}
		}
	}

	unsigned ma=0;
	//unsigned na=0;
	//unsigned pa=0;


	unsigned i=0;
	while(i<min_half ){
		unsigned xi=half_x-i;
		unsigned yi=half_y-i;
		unsigned j=0;
		while(j<4){
			unsigned k=0;
			while(k < (2*i+1)){
	//			cout<<xi<<" "<<yi<<endl;
	//			cout<<"i"<<i<< " 2i+1"<< (2*i+1) <<endl;
				if(_bins_vec[xi][yi]->capacity()>0){
						if(ma<vec_module.size()){
							_bins_vec[j][k]->addModule(vec_module[ma++]);
						}
						else{
							break;
						}
				}
				else{
					if(j==0){
						xi+=1;
					}
					else if(j==1){
						yi+=1;
					}
					else if(j==2){
						xi-=1;
					}
					else if(j==3){
						yi-=1;
					}
					else{
						assert(0);
					}
					k++;
				}
			}
			j++;
		}
		i++;
	}


//	while(i<_placement.numModules()){
	while(ma<vec_module.size()){
		unsigned j=0;
		while(j<_bin_numx){
			unsigned k=0;
			while(k<_bin_numy){
				if(_bins_vec[j][k]->capacity()>0){
				//	cout<<"j:"<<j<<" k:"<<k<<" i:"<<i<<endl;
				//	cout<<_bins_vec[j][k]->size()<<endl;
				//	cout<<_placement.module(i).name()<<endl;
//					_bins_vec[j][k]->addModule(&_placement.module(i++));
//					if(i>=_placement.numModules()){break;}
						if(ma<vec_module.size()){
							_bins_vec[j][k]->addModule(vec_module[ma++]);
						}
						else{
							break;
						}

				}
				else{k++;}
			}

//			if(ma>=vec_module.size()){break;}
			//if(i>=_placement.numModules()){break;}
			j++;
		}
//		if(ma>=vec_module.size()){break;}
		//if(i>=_placement.numModules()){break;}
	}
//	cout<<"***********1*************"<<endl;


}
void GlobalPlacer::randomInitialize(){
	unsigned ma=0;
	while(ma<_placement.numModules()){
		unsigned j=0;
		while(j<_bin_numx){
			unsigned k=0;
			while(k<_bin_numy){
				if(_bins_vec[j][k]->capacity()>0){
//					_bins_vec[j][k]->addModule(&_placement.module(ma++));
//					if(i>=_placement.numModules()){break;}
						if(ma<_placement.numModules()){
							_bins_vec[j][k]->addModule(& _placement.module(ma++));
						}
						else{
							break;
						}

				}
				else{k++;}
			}

//			if(ma>=vec_module.size()){break;}
			//if(i>=_placement.numModules()){break;}
			j++;
		}
//		if(ma>=vec_module.size()){break;}
		//if(i>=_placement.numModules()){break;}
	}
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

void GlobalPlacer::printAllModules(int zi){

	cout<<"***********start_"<<zi<<"*************"<<endl;

	for(unsigned i=0;i<_placement.numModules();i++){
		cout<<_placement.module(i).name()<<endl;
	}

	cout<<"***********end_"<<zi<<"*************"<<endl;
}

void GlobalPlacer::testBin(){

  //  cout << format("HPWL: %.10f", getCost() ) << endl;
//	printAllBins(1);

//	printAllBins(0);
	for(size_t i=1;i<100;i++){
		random_neighbor();
	  //  cout << format("HPWL: %.10f", getCost() ) << endl;
		restore_backup();
//		if(i==2){
//			printAllBins(i*2-1);
//		}
	}

//	printAllBins(1);

//    cout << format("HPWL: %.10f", getCost() ) << endl;
	
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

void GlobalPlacer::simu_anneal(){
	float r=R;
	float t=T;
	float frozen=FROZEN;
	float i=0;
	float n=N;
//	_myusage.reset();

//	_treemgr.printCost();
	while(i<50){
		simu_anneal_sub(r,t,frozen,n);
		t*=0.1;
		frozen*=0.1;
		n*=3;
		i++;
		float val =i;
//		cout << format("HPWL: %.0f", _placement.computeHpwl() ) << endl;
//		cout<<"val: "<<val<<endl;
		if(unsigned(i)%5==0){
			r=R;
			t=T;
			frozen=FROZEN;
			n=N;
			RandomNumGen::change_seed();
		}
	}
//	cout<<"overall times:"<<sa_total_time<<endl;
/*	
	float r=R1;
	float t=T1;
	float frozen=FROZEN1;
	float i=0;
	while(!simu_anneal_sub(r,t,frozen)){
		t*=0.03;
		frozen*=0.03;
		r=sqrt(r);
		float val =i++;
		stringstream ss (stringstream::in | stringstream::out);
		ss << val;
		string test = ss.str();
		string imgstate=_imgname+"_"+string(test)+".raw";
		draw_block(imgstate.c_str());
	}*/
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
		//	double origin_cost=getCost();
	//		cout<<"origin_cost:"<<origin_cost<<endl;
	//		cout<<"before random_neighbor"<<endl;
			double delta_cost=random_neighbor();
	//		cout<<"after random_neighbor"<<endl;
		//	double after_cost=getCost();
	//		cout<<"after_cost:"<<after_cost<<endl;
		/*	if(after_cost==-1){
				complete=true;
		//		cout<<"complete!!!"<<endl;
				break;
			}*/

		//	double delta_cost=after_cost-origin_cost;
		//	cout<<delta_cost_2-delta_cost<<endl;
		//	cout<<"delta_cost:"<<delta_cost<<endl;
		//	cout<<"delta_cost2:"<<delta_cost_2<<endl;
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
//	cout<<"times:"<<times<<endl;
//	cout<<"fall_times:"<<fall_time<<endl;
//	cout<<"climb_times:"<<climb_time<<endl;
//	cout<<"restore_times:"<<restore_time<<endl;
//	cout<<"========report=========="<<endl;
//	_myusage.report(1,1);
//	cout<<"**********************end sa:************************8"<<endl;
//	sa_total_time+=times;
	return complete;
}
double GlobalPlacer::getCost(){
		double total_cost=0;
		for (unsigned i=0; i<_placement.numNets();i++){
			total_cost+=getNetEU(_placement.net(i));
		}
		return total_cost;//_placement.computeHpwl()/100; 
	}
double GlobalPlacer::getNetCost(Module* m1, Module* m2){
	vector<Net*> net_vec;
	if(m1!=0){
		for(unsigned i=0;i<m1->numPins();i++){
			vector<Net*>::iterator it;
			Net* n1=&( _placement.net(m1->pin(i).netId()) );
			it = find (net_vec.begin(), net_vec.end(), n1 );
			if(it==net_vec.end()){
				net_vec.push_back(n1);
			}
		}
	}
	if(m2!=0){
		for(unsigned i=0;i<m2->numPins();i++){
			vector<Net*>::iterator it;
			Net* n2=&( _placement.net(m2->pin(i).netId()) );
			it = find (net_vec.begin(), net_vec.end(), n2 );
			if(it==net_vec.end()){
				net_vec.push_back(n2);
			}
		}
	}
	double net_cost=0;
//	cout<<net_vec.size();
	for(unsigned i=0;i<net_vec.size();i++){
		
		net_cost+=getNetEU(*net_vec[i]);
	}
	return net_cost/_initial_cost;//sqrt(_initial_cost);
}
double GlobalPlacer::random_neighbor(){
//	store_backup();
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
	assert(idx1<_bin_numx);
	assert(idx2<_bin_numx);
	assert(idy1<_bin_numy);
	assert(idy2<_bin_numy);
	double cap1=_bins_vec[idx1][idy1]->capacity();
	double cap2=_bins_vec[idx2][idy2]->capacity();
	unsigned sz1=_bins_vec[idx1][idy1]->size();
	unsigned sz2=_bins_vec[idx2][idy2]->size();
/*	enum Move_state{
			MOVE_M1=0,
			MOVE_M2,
			MOVE_EXCHANGE
		};
	Move_state move_state;*/
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
	//cout<<"move_state:"<<move_state<<endl;
	//cout<<"idx1,idy1:"<<idx1<<" "<<idy1<<endl;
	//cout<<"idx2,idy2:"<<idx2<<" "<<idy2<<endl;
	backup_bin_1=_bins_vec[idx1][idy1];
	backup_bin_2=_bins_vec[idx2][idy2];
	Module* m1=0;
	Module* m2=0;
//	double pre_cost_cap= (-1)*(_bins_vec[idx1][idy1]->capacity())+(-1)*(_bins_vec[idx2][idy2]->capacity());

	switch(move_state){
		case MOVE_M1:{
			m1=_bins_vec[idx1][idy1]->removeModule(_rnGen(sz1));
		//	_bins_vec[idx2][idy2]->addModule(m1);
		//	backup_m_1=m1;

		}break;
		case MOVE_M2:{
			m2=_bins_vec[idx2][idy2]->removeModule(_rnGen(sz2));
		//	_bins_vec[idx1][idy1]->addModule(m2);
		//	backup_m_2=m2;

		}break;
		case MOVE_EXCHANGE:{
			m1=_bins_vec[idx1][idy1]->removeModule(_rnGen(sz1));
			m2=_bins_vec[idx2][idy2]->removeModule(_rnGen(sz2));
		//	_bins_vec[idx2][idy2]->addModule(m1);
		//	_bins_vec[idx1][idy1]->addModule(m2);
		//	backup_m_1=m1;
		//	backup_m_2=m2;
		}break;
		default:{
			assert(0);
		}break;
	}

	backup_m_1=m1;
	backup_m_2=m2;
	
	//cout<<"m1:"<<m1->name()<<endl;
	//cout<<"m2:"<<m2->name()<<endl;

	double pre_cost_net=getNetCost(m1,m2);
	_bins_vec[idx2][idy2]->addModule(m1);
	_bins_vec[idx1][idy1]->addModule(m2);
	double post_cost_net=getNetCost(m1,m2);
//	double post_cost_cap= (-1)*(_bins_vec[idx1][idy1]->capacity())+(-1)*(_bins_vec[idx2][idy2]->capacity());
	double net_cost=post_cost_net-pre_cost_net;
//	double cap_cost=(post_cost_cap-pre_cost_cap);
	//cout<<"net_cost:"<<net_cost<<endl;
	//cout<<"cap_cost:"<<cap_cost<<endl;
	return net_cost;//+cap_cost;
}
void GlobalPlacer::store_backup(){
	
/*	_backup_module.clear();
	for(unsigned j=0;j<_bin_numx;j++){
		for(unsigned k=0;k<_bin_numy;k++){
			for(unsigned i=0;i<_bins_vec[j][k]->size();i++ ){
			//	cout<<"j:"<<j<<" k:"<<k<<" i:"<<i<<endl;
			//	cout<<_bins_vec[j][k]->size()<<endl;
				_backup_module.push_back(_bins_vec[j][k]->getModule(i));
			}
		}
	}*/
}
void GlobalPlacer::restore_backup(){
	switch(move_state){
		//		//cout<< "rstate2:" <<rstate<<endl;
		case MOVE_M1:{
			backup_bin_2->removeModule(backup_m_1);
			backup_bin_1->addModule(backup_m_1);

		}break;
		case MOVE_M2:{
			backup_bin_1->removeModule(backup_m_2);
			backup_bin_2->addModule(backup_m_2);

		}break;
		case MOVE_EXCHANGE:{
			backup_bin_1->removeModule(backup_m_2);
			backup_bin_2->removeModule(backup_m_1);
			backup_bin_1->addModule(backup_m_1);
			backup_bin_2->addModule(backup_m_2);
		}break;
		default:{
			assert(0);
		}break;
	}
	backup_bin_1=0;
	backup_bin_2=0;
	backup_m_1=0;
	backup_m_2=0;

/*
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
*/
}
void GlobalPlacer::confirm_backup(){

}
