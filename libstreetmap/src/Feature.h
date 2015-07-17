/* 
 * File:   Feature.h
 * Author: guofan1
 *
 * Created on February 20, 2015, 12:49 PM
 */
#ifndef FEATURE_H
#define	FEATURE_H
#include <iostream>
#include <string>
#include <vector>
#include "StreetsDatabaseAPI.h"
#include "graphics.h"

using namespace std;
class Feature {
private:
	unsigned featNum; // the index of the feature
	bool isWater; // for determining colour on map
	bool isClosed; // closed shape or line segments
	string featAttr; // type of feature (feature attribute)
	unsigned numPoints; // number of coordinates
	vector<t_point> featCoords; // vector of coordinates, as t_point structs
public:
	//Constructor and Destructor
	Feature(unsigned _featNum, bool _isWater, bool _isClosed, string _featAttr, unsigned _numPoints, vector<t_point> _featCoords);
	~Feature();
	
	//Only need functions to retrieve data, the constructor will set everything.
	unsigned getNum(); // gets the index of the feature
	bool getType(); // gets type (water/ land) of feature
	bool getClosed(); // gets type of shape
	string getAttr(); // gets the actual attribute name
	unsigned getPoints(); // gets number of points
	vector<t_point> getCoords(); // gets the vector of coordinates
};
#endif	/* FEATURE_H */