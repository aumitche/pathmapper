/* 
 * File:   POI.h
 * Author: guofan1
 *
 * Created on February 22, 2015, 16:49 PM
 */

#ifndef POI_H
#define	POI_H
#include <iostream>
#include <string>
#include "StreetsDatabaseAPI.h"
#include "graphics.h"

using namespace std;
class POI {
private:
	string namePOI; // name of point of interest
	t_point coordsPOI; // coordinates of POI, as a t_point
public:
	//Constructor and Destructor
	POI(string _namePOI, t_point _coordsPOI);
	~POI();
	
	//Only need functions to retrieve data, the constructor will set everything.
	string getNamePOI(); // gets the actual attribute name
	t_point getCoordsPOI(); // gets the vector of coordinates
};

#endif	/* POI_H */
