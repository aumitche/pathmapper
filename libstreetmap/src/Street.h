/*
 * File:   Street.h
 * Author: radovic1
 *
 * Created on February 18, 2015, 3:16 PM
 */

#ifndef STREET_H
#define	STREET_H

#include <iostream>
#include <string>
using namespace std;

class Street {
private:
	string name;
	unsigned int id;
	double length;
	unsigned int numSegs;
	double speed;
	unsigned int type; //0 = sideroad, 1 = mainroad, 2 = highway
public:
	//Constructor and Destructor
	Street(string name_, unsigned int id_, double length_, unsigned int numSegs_, double speed_, unsigned int type_);
	~Street();

	//Only need functions to retrieve data, the constructor will set everything.
	string getName();
	unsigned int getType();
	unsigned int getID();
	double getLength();
	unsigned int getNumSegs();
	double getSpeed();
};



#endif	/* STREET_H */
