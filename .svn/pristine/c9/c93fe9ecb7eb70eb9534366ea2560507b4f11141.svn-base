#include "Highway.h"

using namespace std;

Highway::Highway(){
    name = "";
    vector<unsigned> streetIDs_, segIDs_;
    id = streetIDs_;
    streetSegs = segIDs_;
}

Highway::~Highway(){}

string Highway::getName(){
    return name;
}

vector<unsigned> Highway::getStreetID(){
    return id;
}

vector<unsigned> Highway::getSegID(){
    return streetSegs;
}

void Highway::setName(string highwayName){
    name = highwayName;
}

void Highway::setSegID(vector<unsigned> & segIDs){
    streetSegs.insert(streetSegs.end(), segIDs.begin(), segIDs.end());
}

void Highway::setStreetID(unsigned streetID){
    id.push_back(streetID);
}
