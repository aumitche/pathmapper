/*THINGS NEEDED IN POI.cpp*/

#include "POI.h" 
#include "graphics.h"
#include <string>

using namespace std;

//Constructor and Destructor
POI::POI(string _namePOI, t_point _coordsPOI){
    namePOI = _namePOI; // type of feature (feature attribute)
    coordsPOI = _coordsPOI; // vector of coordinates, as t_point structs
}

POI::~POI(){
    // does nothing ayy
}

//Only need functions to retrieve data, the constructor will set everything.

string POI::getNamePOI(){
	return namePOI;
}

t_point POI::getCoordsPOI(){
	return coordsPOI;
}