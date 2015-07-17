#pragma once //protects against multiple inclusions of this header file

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

#include "Feature.h"
#include "Highway.h"
#include "POI.h"
#include "Street.h"

//use these defines whenever you need earth radius
//or conversion from degrees to radians
#define EARTH_RADIUS_IN_METERS 6372797.560856
#define DEG_TO_RAD 0.017453292519943295769236907684886

using namespace std;

extern vector<Street> streets;
//function to load bin or osm map
bool load_map(string map_name);

//close the loaded map
void close_map();

//function to return intersection id for an intersection name
unsigned find_intersection_id_from_name(string intersection_name);

//function to return street id for an intersection name
unsigned find_street_id_from_name(string street_name);

//function to return the street segments for a given intersection
vector<unsigned> find_intersection_street_segments(string intersection_name);
vector<unsigned> find_intersection_street_segments(unsigned intersection_id);

//function to return street names/IDs at an intersection
vector<string> find_intersection_street_names(string intersection_name);
vector<string> find_intersection_street_names(unsigned intersection_id);
vector<unsigned> find_intersection_street_IDs(unsigned intersection_id);

//can you get from intersection1 to intersection2 using a single street segment (hint: check for 1-way streets too)
bool are_directly_connected(string intersection_name1, string intersection_name2);

//find all intersections connected by one street segment from given intersection (hint: check for 1-way streets too)
vector<unsigned> find_adjacent_intersections(string intersection_name);

//for a given street, return all the street segments
vector<unsigned> find_street_street_segments(string street_name);

//for a given street ID, return all the street segments
vector<unsigned> find_street_street_segments_from_id(unsigned int id);

//for a given street, find all the intersections
vector<unsigned> find_all_street_intersections(string street_name);

//find distance between two coordinates
double find_distance_between_two_points(LatLon point1, LatLon point2);

//find the length of a given street segments
double find_street_segment_length(unsigned street_segment_id);

//find the length of a whole street
double find_street_length(string street_name);

//find the length of a whole street, given its ID
double find_street_length_from_id(unsigned int id);

//find the travel time to drive a street segment (time(minutes) = distance(Km)/speed_limit(Km/h)*60)
double find_segment_travel_time(unsigned street_segment_id);

//find the nearest point of interest (by name) to a given position
unsigned int find_closest_point_of_interest(LatLon my_position);
