#define SIDEROADLENGTH 2000 // length is in metres
#define MAINROADLENGTH 5000

#include <iostream>
#include <utility>
#include <vector>
#include <string>
#include <time.h>
#include <cmath>
#include <boost/unordered_map.hpp>

#include "easygl_constants.h"
#include "graphics.h"
#include "IntersectionGraph.h"
#include "StreetsDatabaseAPI.h"
#include "LatLon.h"
#include "Street.h"

using namespace std;

// Creates a new class called "street", includes name, street id, length, number of segments, speed limit

// Constructor
Street::Street(string name_, unsigned int id_, double length_, unsigned int numSegs_, double speed_, unsigned int type_) {
	name = name_;
	id = id_;
	length = length_;
	numSegs = numSegs_;
	speed = speed_;
	type = type_;
}

Street::~Street() {
	//No dynamic stuff. Do nothing.
}

// accessor functions
string Street::getName() {
	return name;
}

unsigned int Street::getType() {
	return type;
}

unsigned int Street::getID() {
	return id;
}

double Street::getLength() {
	return length;
}

unsigned int Street::getNumSegs() {
	return numSegs;
}

double Street::getSpeed() {
  return speed;
}
