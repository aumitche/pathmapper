#include "Feature.h"

using namespace std;

//Constructor and Destructor
Feature::Feature(unsigned _featNum, bool _isWater, bool _isClosed, string _featAttr, unsigned _numPoints, vector<t_point> _featCoords){
    featNum = _featNum; // the index of the feature
    isWater = _isWater; // for determining colour on map
    isClosed = _isClosed; // closed shape or line segments
    featAttr = _featAttr; // type of feature (feature attribute)
    numPoints = _numPoints; // number of coordinates
    featCoords = _featCoords; // vector of coordinates, as t_point structs
}
Feature::~Feature(){
}

//Only need functions to retrieve data, the constructor will set everything.
unsigned Feature::getNum(){
	return featNum;
}

bool Feature::getType(){
	return isWater;
}

bool Feature::getClosed(){
	return isClosed;
}

string Feature::getAttr(){
	return featAttr;
}

unsigned Feature::getPoints(){
	return numPoints;
}

vector<t_point> Feature::getCoords(){
	return featCoords;
}
