/* 
 * File:   GridParam.h
 * Author: radovic1
 *
 * Created on March 20, 2015, 8:15 PM
 */

#ifndef GRIDPARAM_H
#define	GRIDPARAM_H


class gridParam {
    double firstLonBarrier;
    double firstLatBarrier;
    double secondLonBarrier;
    double secondLatBarrier;
    
public:
    
    double getFirstLon();
    double getFirstLat();
    double getSecondLon();
    double getSecondLat();
    
    void setParameters(double firstLon, double firstLat,double secondLon, double secondLat);
    
};



#endif	/* GRIDPARAM_H */

