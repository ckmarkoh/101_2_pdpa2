/******************************************************************************
  FileName     [ rnGen.h ]
  PackageName  [ util ]
  Synopsis     [ Random number generator ]
  Author       [ Chung-Yang (Ric) Huang ]
  Copyright    [ Copyleft(c) 2007-2010 LaDs(III), GIEE, NTU, Taiwan ]
******************************************************************************/
#ifndef RN_GEN_H
#define RN_GEN_H
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <stdlib.h>  
#include <limits.h>
#include <ctime>
//#define SEED 1365591672
#define my_srandom  srandom
#define my_random   random
//#define SEED 1365691387 //e2
//#ifdef RSEED
//#define SEED time(0)
//#define SEED 12 //seg fault
//#endif
//#define SEED 1365691203
//#define SEED 1365631505 //quick success

#define SEED 1369840810 //USE THIS
//#define SEED 1365439588 //bug2 rotate+rotate
//#define SEED 1365438783 //bug1 rotate+exchange
//#define SEED 1365439384 //good1 rotate+exchange
class RandomNumGen
{
   public:
   //   RandomNumGen() { my_srandom(getpid()); }
      RandomNumGen() { 
		  cout<<"SEED:"<<SEED<<endl; 
		  my_srandom(SEED); 
	  }
      RandomNumGen(unsigned seed) { my_srandom(seed); }
      const double operator() (const double range) const {
         return double(range * (double(my_random()) / INT_MAX));
      }
      const unsigned operator() (const int range) const {
         return unsigned(range * (double(my_random()) / INT_MAX));
      }
     /* const unsigned operator() (const unsigned range ) const {
			 return unsigned(range * (double(my_random()) / INT_MAX));
	  }*/
      const unsigned operator() (const unsigned range,bool unbalanced=false) const {

	  	if(unbalanced){
			float half_range=range/2;
			//cout<<"my_random "<<my_random()<<endl;

			float x=(float(my_random())/float(RAND_MAX))*half_range;
			//cout<<"x "<<x<<endl;
			if(my_random()%2==1){
				return unsigned(half_range-sqrt(pow(half_range,2)-pow(x,2)));
			}
			else{
				return unsigned(half_range+sqrt(pow(half_range,2)-pow(x,2)));
			}
		}
		else{
			 return unsigned(range * (double(my_random()) / INT_MAX));
      	}
	  }
	  static void change_seed(){my_srandom(my_random());}

};

#endif // RN_GEN_H

