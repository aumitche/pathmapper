/* 
 * File:   CityBlock.h
 * Author: radovic1
 *
 * Created on March 20, 2015, 6:40 PM
 */
#include <vector>

#ifndef CITYBLOCK_H
#define	CITYBLOCK_H

#define TOP_LEFT 0
#define TOP_MID 1
#define TOP_RIGHT 2
#define MID_LEFT 3
#define CENTER 4
#define MID_RIGHT 5
#define BOTTOM_LEFT 6
#define BOTTOM_MID 7
#define BOTTOM_RIGHT 8

using namespace std;


class cityBlock {
    
    vector<unsigned> here;          // The POIs that exist in these bounds of the city
    vector<unsigned> near;          // The numerical indices of directly connected blocks
    vector<unsigned> far;           // The numerical indices of existing but not directly connected blocks
    
public:
    
    vector<unsigned> getHere();
    vector<unsigned> getNear();
    vector<unsigned> getFar();
    
    void addHere (unsigned int id);
    void addNear (unsigned int idN);
    void addFar (unsigned int idF);
    
    cityBlock(int position);
    cityBlock();
    
    void operator=  ( const cityBlock & rhs ){
        here = rhs.here;
        near = rhs.near;
        far = rhs.far;
    }
    
};


#endif	 /*CITYBLOCK_H*/ 

