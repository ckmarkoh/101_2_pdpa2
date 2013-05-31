#ifndef GLOBALPLACER_H
#define GLOBALPLACER_H

#include "Placement.h"
#include <cassert>
#include <map>
#include <vector>
#include "rnGen.h"

class Bin
{
public:
	Bin(){}
	Bin(double x,double y,unsigned i, unsigned j):_x(x),_y(y),_i(i),_j(j){}
	Module* getModule(unsigned i){
		return _modules_vec[i];
	}
	void print_all(){
		cout<<"start_bin:"<<_i<<","<<_j<<endl;
		for(size_t i=0;i<_modules_vec.size();i++){
			cout<<_modules_vec[i]->name()<<endl;
		}

		cout<<"end_bin:"<<_i<<","<<_j<<endl;
	}
	Module* removeModule(Module* m){
	//	cout<<"before_remove0:"<<endl;
	//	print_all();
		Module* m2=0;
		if(m==0){
			return 0;
		}
		//i = (i < _modules_vec.size()) ? i : _modules_vec.size();
		vector<Module* >::iterator it =_modules_vec.begin();
		while(it!=_modules_vec.end()){
			if((*it)==m){
	//			cout<<"remove1:"<<(*it)->name()<<endl;
				m2=(*it);
				_modules_vec.erase(it);
				break;
				//return (*it);
			}
			else{
				it++;
			}
		}
	//	cout<<"remove2:"<<(*it)->name()<<endl;
	//	assert(it!=_modules_vec.end());
	//	print_all();
		return m2;
	}
	Module* removeModule(unsigned i){
		//cout<<"before_remove0:"<<endl;
		//print_all();
		if(_modules_vec.size()==0){
			return 0;
		}
		i = (i <= _modules_vec.size()) ? i : _modules_vec.size();
		vector<Module* >::iterator it =_modules_vec.begin();
		//cout<<"i:"<<i<<endl;
		for(unsigned j =0 ;j<i;j++){
			//cout<<"r_count:"<<(*it)->name()<<endl;
			it++;
			//cout<<"j"<<j<<endl;
		}
		Module* m2=(*it);
		//cout<<"remove00:"<<m2->name()<<endl;
		_modules_vec.erase(it);
		//cout<<"remove01:"<<m2->name()<<endl;
		//print_all();
		return m2;
	}
	void addModule(Module *m){
		if(m!=0){
			m->setCenterPosition(centerX(),centerY());	
			_modules_vec.push_back(m);
		//	cout<<"add:"<<m->name()<<endl;
		}
	}
	double capacity(double offset=1.0){
		double cap=area();
		for(unsigned i=0;i<_modules_vec.size();i++){
			cap-=(_modules_vec[i]->area()*offset);	
		}
		return cap;
	}
	void clear(){_modules_vec.clear();}

	unsigned size(){return _modules_vec.size();}
	void setXY(double x,double y){_x=x;_y=y;}
    double centerX() {return _x + _width/2;}
	double centerY() {return _y + _height/2;}


	static double width(){return _width;}
	static double height(){return _height;}
	static void setWidthHeight(double w,double h){_width=w;_height=h;}
	static double area(){return _width*_height;}
	
private:
	vector<Module*> _modules_vec;
	double _x;
	double _y;
	static double _width;
	static double _height;
	unsigned _i;
	unsigned _j;
};

class GlobalPlacer 
{
public:
    GlobalPlacer(Placement &placement);
	void place();
	void initialize();
	void prePlaceInitialize();
	void randomInitialize();
	void testBin();
	void printAllBins(int zi=0);
	void printAllModules(int zi=0);
	double getNetWA(Net& n );
	double getNetEU(Net& n );
	double getNetHPWL(Net& n );
	void testNetWA();
	//double sigMoid(Module& m);
	void reset_bin(unsigned numx,unsigned numy){
		_bin_numx=numx;
		_bin_numy=numy;
		_bins_vec= new Bin**[numx];
		for(size_t i=0;i<numx;i++){
			_bins_vec[i]=new Bin*[numy];
			for(size_t j=0;j<numy;j++){
				_bins_vec[i][j]=new Bin( _placement.boundryLeft()+i*Bin::width(), _placement.boundryBottom()+j*Bin::height(),i,j);
			}
		}
	}
	double getCost();
	double getNetCost(Module* m1, Module* m2=0);

	void simu_anneal();
	bool simu_anneal_sub(double r,double t,double frozen,double n);
	double random_neighbor();
	void restore_backup();

	void store_backup();
	void confirm_backup();
private:
	enum Move_state{
			MOVE_M1=0,
			MOVE_M2,
			MOVE_EXCHANGE
		};
	Move_state move_state;

	RandomNumGen _rnGen;
    Placement& _placement;
	Bin*** _bins_vec;
	vector<Bin> _full_bin;
	vector<Bin> _unfull_bin;
	unsigned _bin_numx;
	unsigned _bin_numy;
	vector<Module*> _backup_module;
	double _initial_cost;
//	unsigned backup_id1;
//	unsigned backup_id2;
	Bin*	backup_bin_1;
	Bin*	backup_bin_2;	
	Module* backup_m_1;
	Module* backup_m_2;
//	Bin*	backup_bin_t_1;	
//	Bin*	backup_bin_t_2;	
};



#endif // GLOBALPLACER_H
