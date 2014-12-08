/**
 * @file histogrammer.h
 *
 *  Created on: 30.06.2014
 *      Author: Georg Bergner <G.Bergner@uni-meunster.de>
 * Copyright: See COPYING file that comes with this distribution
 *
 *
 */

#ifndef HISTOGRAMMER_H_
#define HISTOGRAMMER_H_
#include <vector>
#include <iostream>
#include "basicdef.h"

class Histogrammer {
private:
	typedef unsigned long int Counts;
	std::vector<Counts> data_;
	std::vector<Real> limits_;
	Counts total_;
public:
	Histogrammer(const Real min,const Real max,const size_t numsteps): data_(numsteps+1,0),limits_(numsteps+1,0),total_(0){
      initialize(min,max,numsteps);
	}

	void reset(){
		for(size_t i(data_.size());i--;){
			data_[i]=0;
		}
		total_=0;
	}

	void initialize(const Real min,const Real max,const size_t numsteps){
		data_.resize(numsteps+1);
		reset();
		limits_.resize(numsteps+1);
		if(limits_.size()<=2) return;
		if(limits_.size()==2){
			limits_[0]=min;
			limits_[1]=max;
			return;
		}
		Real dif=fabs(max-min)/Real(numsteps);
		Real cm(max);
		for(size_t i(limits_.size());i--;){
	      limits_[i]=cm;
	      cm-=dif;
		}
	}
	virtual ~Histogrammer(){

	}

	void addData(const Real val){
		total_++;
		if(data_.size()==1) return;
		size_t i(data_.size()-1);
		bool notfound(true);
		while(i-- && notfound){
			if(val>=limits_[i]&&val<limits_[i+1]){
				data_[i]++;
				notfound=false;
			}
		}
		if(notfound){
			data_[data_.size()-1]++;
		}
	}

	void print(std::ostream & os){
		if(data_.size()==0) return;
		Counts out(data_[data_.size()-1]);
		Real weight(1.0/Real(total_-out));
		if(total_<=out){
			weight=1;
		}
		for(size_t i(0);i+1<data_.size();i++){
		 os<<limits_[i]<<" "<<weight*data_[i]<<" "<<data_[i]<<"\n";
		}
		os<<"#outside: "<<out<<"\n";
	}
};

class Histogrammer2D {
private:
	typedef unsigned long int Counts;
	std::vector<Counts> data_;
	std::vector<Real> limitsx_;
	std::vector<Real> limitsy_;
	Counts total_;
public:
	Histogrammer2D(const Real minx,const Real maxx,const size_t numstepsx,const Real miny,const Real maxy,const size_t numstepsy): data_((numstepsx*numstepsy)+1,0),limitsx_(numstepsx+1,0),limitsy_(numstepsx+1,0),total_(0){
      initialize(minx,maxx,numstepsx,miny,maxy,numstepsy);
	}

	void reset(){
		for(size_t i(data_.size());i--;){
			data_[i]=0;
		}
		total_=0;
	}

	void initialize(const Real minx,const Real maxx,const size_t numstepsx,const Real miny,const Real maxy,const size_t numstepsy){
		data_.resize(numstepsx*numstepsy+1);
		reset();
		setLimits(minx,maxx,numstepsx,limitsx_);
		setLimits(miny,maxy,numstepsy,limitsy_);
	}

	void setLimits(const Real min, const Real max, const size_t numsteps, std::vector<Real>& limits ) const{
		limits.resize(numsteps+1);
		if(limits.size()<=2) return;
		if(limits.size()==2){
			limits[0]=min;
			limits[1]=max;
			return;
		}
		Real dif=fabs(max-min)/Real(numsteps);
		Real cm(max);
		for(size_t i(limits.size());i--;){
	      limits[i]=cm;
	      cm-=dif;
		}
	}

	virtual ~Histogrammer2D(){

	}

	void addData(const Real valx,const Real valy){
		total_++;
		if(data_.size()==1 || limitsx_.size()==0 || limitsy_.size()==0) return;
		const size_t numx(limitsx_.size()-1);
		size_t i(data_.size()-1);
		bool notfound(true);
		while(i-- && notfound){
			const size_t yi(i/numx);
			const size_t xi(i-yi*numx);
			if(valx>=limitsx_[xi] && valx<limitsx_[xi+1] && valy>=limitsy_[yi] && valy<limitsy_[yi+1]){
				data_[i]++;
				notfound=false;
			}
		}
		if(notfound){
			data_[data_.size()-1]++;
		}
	}

	void print(std::ostream & os){
		if(data_.size()<=1 || limitsx_.size()==0 || limitsy_.size()==0) return;
		const size_t numx(limitsx_.size()-1);
		Counts out(data_[data_.size()-1]);
		Real weight(1.0/Real(total_-out));
		if(total_<=out){
			weight=1;
		}
		for(size_t i(0);i+1<data_.size();i++){
	     const size_t yi(i/numx);
		 const size_t xi(i-yi*numx);
		 os<<limitsx_[xi]<<" "<<limitsy_[yi]<<" "<<weight*data_[i]<<" "<<data_[i]<<"\n";
		}
		os<<"#outside: "<<out<<"\n";
	}
};

#endif /* HISTOGRAMMER_H_ */
