/* 
 * File:   intersectInfo.h
 * Author: radovic1
 *
 * Created on April 1, 2015, 8:51 PM
 */
#include <vector>

#ifndef INTERSECTINFO_H
#define	INTERSECTINFO_H

using namespace std;

struct intersectInfo {
    
    vector<unsigned> path;
    double time;
    
    intersectInfo(vector<unsigned> path_, double time_){
        path= path_;
        time = time_;
    }
    
    intersectInfo(){
        time = 0;
    }
    
      const bool operator>(const intersectInfo & rhs) const {
        return (time > rhs.time);
    }

    const bool operator<(const intersectInfo & rhs) const {
        return (time < rhs.time);
    }

};


#endif	/* INTERSECTINFO_H */

