#define SIDEROADLENGTH 2000 // length is in metres
#define MAINROADLENGTH 5000

#include "m1.h"
#include "m2.h"
#include "m3.h"
#include "m4.h"

using namespace std;

// We want to use the unordered map to reduce computation time
unordered_map <string, unsigned> intersectNameToID; // unordered map of intersection names and IDs
unordered_map <string, unsigned> streetToID; // unordered map of street names and IDs
unordered_map <string, vector<unsigned>> intersectionsForStreet;
gridParam parameters;
vector<cityBlock> cityGrid;


// These global variables will be needed in multiple functions.
vector<Street> streets;
vector<vector<unsigned>> streetIDs;

double averageSpeed(unsigned streetID);
unsigned find_intersection_id_from_name(string intersection_name);
unsigned find_street_id_from_name(string street_name);
vector<unsigned> find_intersection_street_segments(string intersection_name);
vector<unsigned> find_intersection_street_segments(unsigned intersection_id);
vector<string> find_intersection_street_names(string intersection_name);
vector<string> find_intersection_street_names(unsigned intersection_id);
vector<unsigned> find_intersection_street_IDs(unsigned intersection_id);
bool are_directly_connected(string intersection_name1, string intersection_name2);
vector<unsigned> find_adjacent_intersections(string intersection_name);
vector<unsigned> find_street_street_segments(string street_name);
vector<unsigned> find_all_street_intersections(string street_name);
double find_distance_between_two_points(LatLon point1, LatLon point2);
double find_street_segment_length(unsigned street_segment_id);
double find_street_length(string street_name);
double find_segment_travel_time(unsigned street_segment_id);
unsigned int find_closest_point_of_interest(LatLon my_position);


// load the map
bool load_map(string map_name) /*done*/ {

    bool load_success = loadStreetDatabaseBIN(map_name);

    // Our data structures
    // Apparently all these can be treated as arrays - look into this more later (consider using iterators)


    //We will prepare a vector, each element of which contains an integer vector (street id).
	//This for loop populates it with empty vectors.
	//This will be extremely useful in find_street_street_segments.
	for (unsigned int i = 0; i < getNumberOfStreets(); i++) {
		vector<unsigned> segs;	//All of these vectors are empty.
		streetIDs.push_back(segs);
	}

    //Now we are going to actually sort these ids into the appropriate vector element.
    unsigned int streetNum;
    for (unsigned int i = 0; i < getNumberOfStreetSegments(); i++) {
        streetNum = getStreetSegmentStreetID(i);
        streetIDs[streetNum].push_back(i);	//Remember, each element of streetIDs is itself a vector!
    }

    // Makes an unordered map for the intersections (names + ID)
    for (unsigned int i = 0; i < getNumberOfIntersections(); i++) {
        string hashedIntersection;
        string partOfIntersection;
        int count = 0;
        string intersectionName = getIntersectionName(i);

        int ampersandPos = intersectionName.find_first_of("&");    //Finds the location of the first ampersand in the intersection name (if it exists)

        do {
            partOfIntersection = intersectionName.substr(0, ampersandPos - 1);
            if (count == 0) { //For the first instance of a name, simply assign the hashed intersection name to be the first name
                hashedIntersection = partOfIntersection;
                count++;
            } else { // Compare the string values and append them so that they are in descending order (from start to end)
                if (hashedIntersection > partOfIntersection) {
                    hashedIntersection = hashedIntersection.append((" &" + partOfIntersection));
                } else {
                    hashedIntersection = partOfIntersection.append((" &" + hashedIntersection));
                }
            }

            intersectionName.erase(0, ampersandPos + 1); // Erase the part of the intersection name that has been appended
            ampersandPos = intersectionName.find_first_of("&");
            if (ampersandPos == -1) { //If there are no more ampersands in the name, compare and append the name as necessary
                if (hashedIntersection > intersectionName) {
                    hashedIntersection = hashedIntersection.append((" &" + intersectionName));
                } else {
                    hashedIntersection = intersectionName.append((" &" + hashedIntersection));
                }
            }

        } while (ampersandPos != -1); // while there are ampersands in the name, continue to acquire the position of the next ampersand in the name

        intersectNameToID[hashedIntersection] = i; // Assign the map at the address of the hashed intersection to the appropriate intersection ID
    }



    // Make an unordered map for the streets (name + ID)
    // Iterate through the entirety of the street names, and create an unordered map that associates the name (as an index) to the ID
    for (unsigned int i = 0; i < getNumberOfStreets(); i++) {
        string streetName = getStreetName(i);
        streetToID[streetName] = i;
    }

    //We will prepare a vector, each element of which contains an integer vector (street id).
    //This for loop populates it with empty vectors.
    //This will be extremely useful in find_street_street_segments.
    for (unsigned int i = 0; i < getNumberOfStreets(); i++) {
        vector<unsigned> segs; //All of these vectors are empty.
        streetIDs.push_back(segs);
    }

//    //Now we are going to actually sort these ids into the appropriate vector element.
//    unsigned int streetNum;
//    for (unsigned int i = 0; i < getNumberOfStreetSegments(); i++) {
//        streetNum = getStreetSegmentStreetID(i);
//        streetIDs[streetNum].push_back(i); //Remember, each element of streetIDs is itself a vector!
//    }


    for (unsigned int i = 0; i < getNumberOfStreets(); i++) {
        string streetName = getStreetName(i);
        vector<unsigned> intersection_vect; // Create a vector to hold the ids for all the intersections connected to a given street
        vector<unsigned> street_segments; // Find all the street segments that are associated with a street of a given name
        street_segments = find_street_street_segments(streetName);
        StreetSegmentEnds checking_ends;
        for (unsigned int j = 0; j < street_segments.size(); j++) {
            checking_ends = getStreetSegmentEnds(street_segments[j]); // Iterate through each of the street segments in the segments vector and push them back onto the intersection vector
            intersection_vect.push_back(checking_ends.to);
            intersection_vect.push_back(checking_ends.from);
        }

        std::sort(intersection_vect.begin(), intersection_vect.end()); // After all of the segments have been placed in the vector, we wish to sort the vector then delete any duplicates
        intersection_vect.erase(std::unique(intersection_vect.begin(), intersection_vect.end()), intersection_vect.end()); //Removes duplicate IDs
        intersectionsForStreet[streetName] = intersection_vect; // Use the street name to index the map and store the associated vector
    }

    vector<Street> sideroad;
    vector<Street> mainroad;
    vector<Street> highway;

    // Make a vector of Streets where the index is the street id
    unsigned int i = 0;
    for (i = 0; i < getNumberOfStreets(); i++) {
        string name = getStreetName(i);
        double length;
        unsigned int numSegs = find_num_segs(i);
        double speed = averageSpeed(i);
//        if (i % 100 == 0) cout << name << " " << speed << endl;
        unsigned int type;
        if (!i) {
            length = 1;
            type = 0;
        } else {
            length = find_street_length_from_id(i);
            if (speed >= 85) {
                type = 2;
            } else if (length > MAINROADLENGTH) {
                type = 1;
            } else {
                type = 0;
            }
        }
        Street street(name, i, length, numSegs, speed, type);
        streets.push_back(street);

        if (type == 2) {
            highway.push_back(street);
        } else if (type == 1) {
            mainroad.push_back(street);
        } else {
            sideroad.push_back(street);
        }
    }

    cout << "INFO: Sorting features...";
    sortFeatures();
    cout << "DONE" << endl;
    cout << "INFO: Sorting POIs...";
    sortPOI();
    cout << "DONE" << endl;
    cout << "INFO: Filtering highways...";
    filterHighways(highway);
    cout << "DONE" << endl;

	//Draw the map in an EasyGL window
    cout << "INFO: Drawing map...";
	draw_map("Mapper", sideroad, mainroad, highway);
    cout << "DONE" << endl;

    cout << "INFO: Load successful" << endl;
    return load_success;
}

//close the map

void close_map() {
    closeStreetDatabase();
    cout << "INFO: Window closed" << endl;
    // destroy any data structures you created in load_map
    // ...
}

// Finds the average speed limit for a given street
double averageSpeed(unsigned streetID) {
    vector<unsigned> segs = streetIDs[streetID];
    double averageStreetSpeed = 0.0;
//    vector<unsigned> streetSegments = find_street_street_segments_from_id(streetID);
    for (unsigned int i = 0; i < segs.size(); i++) {
//        if (streetID % 100 == 0) cout << getStreetSegmentSpeedLimit(i) << ", ";
        averageStreetSpeed += getStreetSegmentSpeedLimit(segs[i]);
    }
//    if (streetID % 100 == 0) cout << endl;
    averageStreetSpeed /= segs.size();
    return (averageStreetSpeed);
}

//function to return intersection id for an intersection name

unsigned find_intersection_id_from_name(string intersection_name) /*TESTED*/ {

    // This function implements hashing of the input name identical to that done in load map
    // As a result we can simply index into the unordered map to receive the street ID
    // Cascading if statements are used to sort the street name in descending order if there are three involved
    string hashedIntersection;
    string partOfIntersection;
    int count = 0;
    int ampersandPos = intersection_name.find_first_of("&");

    do {
        partOfIntersection = intersection_name.substr(0, ampersandPos - 1);
        if (count == 0) {
            hashedIntersection = partOfIntersection;
            count++;
        } else {
            if (hashedIntersection > partOfIntersection) {
                hashedIntersection = hashedIntersection.append((" &" + partOfIntersection));
            } else {
                hashedIntersection = partOfIntersection.append((" &" + hashedIntersection));
            }
        }

        intersection_name.erase(0, ampersandPos + 1);
        ampersandPos = intersection_name.find_first_of("&");
        if (ampersandPos == -1) {
            if (hashedIntersection > intersection_name) {
                hashedIntersection = hashedIntersection.append((" &" + intersection_name));
            } else {
                hashedIntersection = intersection_name.append((" &" + hashedIntersection));
            }
        }

        } while (ampersandPos != -1);

    return (intersectNameToID[hashedIntersection]);
}

// find street ID from name

unsigned find_street_id_from_name(string street_name) /*TESTED*/ {
    //Remember that hack we did in load_map? Well here it is!
    return (streetToID[street_name]);
} /*done*/

vector<unsigned> find_intersection_street_segments(string intersection_name) /*TESTED*/ {
    //convert the intersection name to an id and use the other find_intersection_street_segments function

    unsigned int isID; //intersection ID
    isID = find_intersection_id_from_name(intersection_name);

    vector<unsigned> intersections;
    intersections = find_intersection_street_segments(isID);
    return intersections;
}

//get the IDs of all street segments meeting at a given intersection

vector<unsigned> find_intersection_street_segments(unsigned intersection_id) /*TESTED*/ {
    // MOST LIKELY THE PREVIOUS FUNCTION CAN CALL THIS ONE and that saves some space in coding...

    // 1. call getIntersectionStreetSegmentCount(unsigned intersectionID) from the API
    // to set size of the vector
    // 2. create vector of this size
    // 3. using some sort of loop, call getIntersectionStreetSegment(unsigned intersectionID,unsigned idx)
    // from the API to get the street segment ID and add to the vector
    // 4. return the vector

    vector<unsigned> intersections;
    for (unsigned int i = 0; i < getIntersectionStreetSegmentCount(intersection_id); i++) {
        intersections.push_back(getIntersectionStreetSegment(intersection_id, i));
    }

    return intersections;
}

//find distance between two coordinates
double find_distance_between_two_points(LatLon point1, LatLon point2) /*TESTED*/ {
    //First convert all latitudes and longitudes to radians (only radians work in math functions).
    double lat1 = point1.lat * DEG_TO_RAD;
    double lon1 = point1.lon * DEG_TO_RAD;
    double lat2 = point2.lat * DEG_TO_RAD;
    double lon2 = point2.lon * DEG_TO_RAD;

	//The following procedure was provided in the documentation for Milestone 1.
    double avgLat = (lat1 + lat2) / 2; //Average latitude

	double cosine = cos(avgLat);

	double x1 = lon1 * cosine;
	double y1 = lat1;
	double x2 = lon2 * cosine;
	double y2 = lat2;

    double distance = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);
    distance = sqrt(distance);
    distance *= EARTH_RADIUS_IN_METERS;

    return distance;
}

//function to return street names at an intersection, given its name

vector<string> find_intersection_street_names(string intersection_name) /*TESTED*/ {
    unsigned int id = find_intersection_id_from_name(intersection_name);
    vector<string> streetNames = find_intersection_street_names(id); //Use the id version to do all the work
    return streetNames;
}

//If we want the street names, simply get a vector of IDs and convert them.
vector<string> find_intersection_street_names(unsigned intersection_id) {
    vector<unsigned> streetIDs = find_intersection_street_IDs(intersection_id);
    vector<string> streetNames;
    for (unsigned int i = 0; i < streetIDs.size(); i++) {
        streetNames.push_back(getStreetName(streetIDs[i]));
    }
    return streetNames;
}

vector<unsigned> find_intersection_street_IDs(unsigned intersection_id) /*TESTED*/ {
    vector<unsigned> segments = find_intersection_street_segments(intersection_id);
    vector<unsigned> streetIDs;
    unsigned int streetID;

	//The vector is currently empty. Get the name of the first segment and throw it in.
	streetID = getStreetSegmentStreetID(segments[0]);
	streetIDs.push_back(streetID);
	bool found = false; //this is a flag which is raised if the segment currently being examined belongs to a street already added to the vector.
	for (unsigned int i = 0; i < segments.size(); i++) {
		found = false;
		for (unsigned int j = 0; j < streetIDs.size(); j++) {
			streetID = getStreetSegmentStreetID(segments[i]);
			if (streetID == streetIDs[j]) { //This street had another segment already recorded.
				found = true;
			}
		}
		if (!found) { //If this street has not been seen before, add it to the list.
			streetID = getStreetSegmentStreetID(segments[i]);
			streetIDs.push_back(streetID);
		}
	}
    return streetIDs;
}

//can you get from intersection1 to intersection2 using a single street segment (hint: check for 1-way streets too)

bool are_directly_connected(string intersection_name1, string intersection_name2) /*FINISHED, NOT TESTED*/ {
	vector<unsigned> adjacentInts = find_adjacent_intersections(intersection_name1); // find the intersections directly connected to intersection 1

	//First check if the user decided to input the same intersection twice. If so, return true.
	if (intersection_name1 == intersection_name2) {
		return true;
	}

	unsigned Int2ID = find_intersection_id_from_name(intersection_name2); // gets the ID of intersection 2 for comparison later
	for (unsigned i = 0; i<adjacentInts.size(); i++) {
		if (adjacentInts[i] == Int2ID) {
			return true; // if ID of intersection 2 matches that of one of the IDs in the vector adjacentInts, then connected
		}
	}
    //If we survive to this point, then no matches have been found. Return false.
    return false;
} //NOT TESTED

//find all intersections connected by one street segment from given intersection (hint: check for 1-way streets too)

vector<unsigned> find_adjacent_intersections(string intersection_name) {
    unsigned isID = find_intersection_id_from_name(intersection_name); // get intersection ID from the name
    vector<unsigned> segments = find_intersection_street_segments(isID); //create+initialize vector holding street segments at the intersection
    vector<unsigned> intersections;
    // find the ends of each street and add it if it's not the source intersection
    // THE ONLY THINGS THAT DON'T GET ADDED ARE STREET SEGMENTS THAT ARE ONEWAY AND LEAD TO THE QUERIED INTERSECTION
    for (unsigned i = 0; i < segments.size(); i++) { // loop through the
        StreetSegmentEnds segEnds = getStreetSegmentEnds(segments[i]);
        if (segEnds.from == isID) { // checks if FROM is the intersection
            intersections.push_back(segEnds.to); // adds TO
        } else {
            if (!getStreetSegmentOneWay(segments[i])) {
                intersections.push_back(segEnds.from); // adds FROM if not one way
            }
        }
    }
    return intersections;
} /*NOT DONE*/

//for a given street, return all the street segments
vector<unsigned> find_street_street_segments(string street_name) /*TESTED*/ {
    unsigned int streetID;
    //Use the unordered_map to return the street id.
    streetID = streetToID[street_name];

    return streetIDs[streetID];
}

vector<unsigned> find_street_street_segments_from_id(unsigned int id) {
    return streetIDs[id];
}

//for a given street, find all the intersections

vector<unsigned> find_all_street_intersections(string street_name) {
    // 0. create the intersection vector
    // 1. call the above thing to get all the street segments (creates a segment vector)
    // 2. go through the vector - call getStreetSegmentEnds(unsigned streetSegmentID) on each segment
    // 3. check if duplicate: if not duplicated, add to intersection vector
    // 4. return intersection vector
    return (intersectionsForStreet[street_name]);
} //NOT FINISHED

//find the length of a given street segments

double find_street_segment_length(unsigned street_segment_id) /*TESTED*/ {
	double distance = 0; // creates and initializes distance variable to 0
	StreetSegmentEnds segment = getStreetSegmentEnds(street_segment_id); // gets the ends of the street segments
	int curveCount = getStreetSegmentCurvePointCount(street_segment_id); // gets total number of curve points
	LatLon intStart = getIntersectionPosition(segment.from); // gets the coordinates of the first intersection
	LatLon intEnd = getIntersectionPosition(segment.to); // gets the coordinates of the second intersection

        // finds the length of the street segment if no curve points
	if (curveCount == 0) {
            distance = find_distance_between_two_points(intStart, intEnd); // finds the distance between the two intersections, adds it to total
	} else { //finds the length of the street segment if there are curve points
		LatLon C1, C2; // initializes two lat/lon coordinates
		// loops through the number of curve points
		for (int i = 1; i< curveCount; i++) {
			C1 = getStreetSegmentCurvePoint(street_segment_id, i); // gets the lat/lon of the current point
			C2 = getStreetSegmentCurvePoint(street_segment_id, i - 1); // gets the lat/lon of the previous point
			distance += find_distance_between_two_points(C1, C2); // finds the distance between the two points, adds it to total
	    }
		C1 = getStreetSegmentCurvePoint(street_segment_id, 0); // gets the lat/lon of the first curve point
		C2 = getStreetSegmentCurvePoint(street_segment_id, curveCount - 1); // gets the lat/lon of the last curve point
		distance += find_distance_between_two_points(C1, intStart); // finds the distance between the first intersection and first curve point
		distance += find_distance_between_two_points(intEnd, C2); // finds the distance between the end intersection and last curve point
	}
    return distance; // returns the distance
} // done

//find the length of a whole street

double find_street_length(string street_name) /*TESTED*/ {
    double streetLength = 0.0; // creates + initializes variable to hold street length
    vector<unsigned> segments = find_street_street_segments(street_name); // get a vector of all the street segments

    // iterates through the vector of street segments
    for (vector<unsigned>::iterator it = segments.begin(); it != segments.end(); ++it) {
        streetLength += find_street_segment_length(*it); // calls "find_street_segment_length" on the vector element pointed at
    }
    return streetLength;
} // done

double find_street_length_from_id(unsigned int id) {
    vector<unsigned> segs = streetIDs[id];
    double length = 0.0;

    for (vector<unsigned>::iterator it = segs.begin(); it != segs.end(); ++it) {
        length += find_street_segment_length(*it); // calls "find_street_segment_length" on the vector element pointed at
    }
    return length;
}

//find the travel time to drive a street segment (time(minutes) = distance(Km)/speed_limit(Km/h)*60)

double find_segment_travel_time(unsigned street_segment_id) /*TESTED*/ {
    // I couldn't find the units explicitly stated for anything, but I assume that
    // the distances were given in KM and the speed limits would be KM/H
    // I'm guessing that the "60" there in the above comment converts it to minutes
    // But running the unit tests led me off by a factor of 1000
    // So I'm guessing the distance was originally in meters
    // in any case, it passes the unit tester now
    // the 0.060 takes into account that discrepancy of M -> KM conversion and HOUR -> MIN
    // I'd have defined it but honestly I have no idea why this is
    // If only people told me what units for distances/ speed limits were.. then
    // there wouldn't be use for this magic number and this humongous comment
    double distance = find_street_segment_length(street_segment_id); // variable for length of segment (km)
    double speedLimit = getStreetSegmentSpeedLimit(street_segment_id); // get the speed limit in km/h
    double time = distance / speedLimit * 0.060; // variable for time it takes to travel the segment, in minutes (*60)
    return time; // time is returned in minutes
}

//find the nearest point of interest (by name) to a given position

unsigned int find_closest_point_of_interest(LatLon my_position) {
	LatLon here = my_position;
	LatLon interesting; //placeholder for point of interest being examined in the loop
	unsigned long long poi = getNumberOfPointsOfInterest();
	double distance;
	double closest;
    unsigned int closestID;

	for (unsigned int i = 0; i < poi; i++) {
		interesting = getPointOfInterestPosition(i);
		distance = find_distance_between_two_points(here, interesting); //find the distance from here to interesting
		if (i == 0) { //This is the first point whose distance we are comparing. So far it is the closest.
			closest = distance;
			closestID = i;
		} else {
			if (distance < closest) {
				closest = distance;	//This is the new closest point of interest. The previous one is usurped!
				closestID = i; //Update the name of the new closest point.
			}
		}
	}
    return closestID;
} // NOT TESTED
