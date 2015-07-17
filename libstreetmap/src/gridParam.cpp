#include "gridParam.h"

double gridParam::getFirstLon(){
    return firstLonBarrier;
}

double gridParam::getFirstLat(){
    return firstLatBarrier;
}

double gridParam::getSecondLon(){
    return secondLonBarrier;
}

double gridParam::getSecondLat(){
    return secondLatBarrier;
}
    
void gridParam::setParameters(double firstLon, double firstLat, double secondLon, double secondLat){
    firstLonBarrier = firstLon;
    firstLatBarrier = firstLat;
    secondLonBarrier = secondLon;
    secondLatBarrier = secondLat;
}
