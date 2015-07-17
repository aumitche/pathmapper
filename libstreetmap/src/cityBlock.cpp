#include <iostream>
#include <utility>
#include <cstdlib>
#include <vector>
#include <string>
#include <time.h>
#include <cmath>
#include <readline/readline.h>
#include <readline/history.h>
#include <boost/unordered_map.hpp>

#include <queue>
#include <set>
#include <algorithm>
#include <functional>

#include "easygl_constants.h"
#include "graphics.h"
#include "IntersectionGraph.h"
#include "StreetsDatabaseAPI.h"
#include "LatLon.h"

#include "Feature.h"
#include "Highway.h"
#include "POI.h"
#include "Street.h"
#include "cityBlock.h"

vector<unsigned> cityBlock::getHere(){
    return here;
}

vector<unsigned> cityBlock::getNear(){
    return near;
}

vector<unsigned> cityBlock::getFar(){
    return far;
}
    
void cityBlock::addHere (unsigned int id){
    here.push_back(id);
}

void cityBlock::addNear (unsigned int idN){
    near.push_back(idN);
}
void cityBlock::addFar (unsigned int idF){
    far.push_back(idF);
}

cityBlock::cityBlock(int position){
    switch(position){
            case TOP_LEFT:
                near = { TOP_MID, CENTER, MID_LEFT};
                far = { TOP_RIGHT, MID_RIGHT, BOTTOM_RIGHT, BOTTOM_MID, BOTTOM_LEFT };
            
            case TOP_MID:
                near = { TOP_LEFT, MID_LEFT, CENTER, MID_RIGHT, TOP_RIGHT };
                far = { BOTTOM_LEFT, BOTTOM_MID, BOTTOM_RIGHT};
                
            case TOP_RIGHT:
                near = { TOP_MID, CENTER, MID_RIGHT };
                far = {TOP_LEFT, MID_LEFT, BOTTOM_LEFT, BOTTOM_MID, BOTTOM_RIGHT} ;
                
            case MID_LEFT:
                near = { TOP_LEFT, TOP_MID, CENTER, BOTTOM_MID, BOTTOM_LEFT }; 
                far = {TOP_RIGHT, MID_RIGHT, BOTTOM_RIGHT };
                
            case CENTER:
                near = { TOP_LEFT, TOP_MID, TOP_RIGHT, MID_RIGHT, BOTTOM_RIGHT, BOTTOM_MID, BOTTOM_LEFT, MID_LEFT };
                
            case MID_RIGHT:
                near = { TOP_RIGHT, TOP_MID, CENTER, BOTTOM_MID, BOTTOM_RIGHT };
                far = { TOP_LEFT, MID_LEFT, BOTTOM_LEFT };
                
            case BOTTOM_LEFT:
                near = { MID_LEFT, CENTER, BOTTOM_MID };
                far = { TOP_LEFT, TOP_MID, TOP_RIGHT, MID_RIGHT, BOTTOM_RIGHT };
                
            case BOTTOM_MID:
                near = { BOTTOM_LEFT, MID_LEFT, CENTER, MID_RIGHT, BOTTOM_RIGHT };
                far = { TOP_LEFT, TOP_MID, TOP_RIGHT };
                
            case BOTTOM_RIGHT:
                near = { MID_RIGHT, CENTER, BOTTOM_MID };
                far = { TOP_RIGHT, TOP_MID, TOP_LEFT, MID_LEFT, BOTTOM_LEFT };    
               
            default:                        // This will never happen, this exists solely for the purpose of the compiler not crying
                near = {0};
                far = {0};
    }      
}

cityBlock::cityBlock(){};           // This should never be called but it's here for safety purposes

//void cityBlock::operator=  ( const cityBlock & rhs ){
//        here = rhs.here;
//        near = rhs.near;
//        far = rhs.far;
//    }
