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
	Bin(double x,double y):_x(x),_y(y){}
	Module* getModule(unsigned i){
		return _modules_vec[i];
	}

	Module* removeModule(unsigned i){
		i = (i < _modules_vec.size()) ? i : _modules_vec.size();
		vector<Module* >::iterator it =_modules_vec.begin();
		for(unsigned j =0 ;j<i;j++){
			it++;
		}
		_modules_vec.erase(it);
		return (*it);
	}
	void addModule(Module *m){
		m->setCenterPosition(centerX(),centerY());	
		_modules_vec.push_back(m);
	
	}
	double capacity(){
		double cap=area();
		for(unsigned i=0;i<_modules_vec.size();i++){
			cap-=_modules_vec[i]->area();	
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
		
};

class GlobalPlacer 
{
public:
    GlobalPlacer(Placement &placement);
	void place();
	void initialize();
	void testBin();
	void printAllBins(int zi=0);
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
				_bins_vec[i][j]=new Bin( _placement.boundryLeft()+i*Bin::width(), _placement.boundryBottom()+j*Bin::height());
			}
		}
	}
	double getCost();
	bool simu_anneal_sub(double r,double t,double frozen,double n);
	void random_neighbor();
	void restore_backup();

	void store_backup();
	void confirm_backup();
private:

	RandomNumGen _rnGen;
    Placement& _placement;
	Bin*** _bins_vec;
	vector<Bin> _full_bin;
	vector<Bin> _unfull_bin;
	unsigned _bin_numx;
	unsigned _bin_numy;
	vector<Module*> _backup_module;

//	Bin*	backup_bin_f_1;
//	Bin*	backup_bin_f_2;	
//	Bin*	backup_bin_t_1;	
//	Bin*	backup_bin_t_2;	
};



#endif // GLOBALPLACER_H
