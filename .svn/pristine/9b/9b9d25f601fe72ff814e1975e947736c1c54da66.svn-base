#include "intersectInfo.h"
#include "gridParam.h"
#include "cityBlock.h"
#include <queue>
#include <set>
#include <algorithm>
#include <functional>

#include "m1.h"
#include "m2.h"
#include "m3.h"
#include "m4.h"

#define NO_EDGE -1
#define NO_TIME 0
#define TURN_OFFSET 0.25

#define OPEN_SET 0
#define CLOSED_SET 1
#define PENDING_SET -1

#define DEFAULT_DISTANCE 10000000

extern gridParam parameters;
extern vector<cityBlock> cityGrid;

using namespace std;

struct queueIntersect {
    unsigned int intersectID;
    double travelTime;

    queueIntersect(unsigned id, double time) {
        intersectID = id;
        travelTime = time;
    }

    queueIntersect() {
        intersectID = -1;
        travelTime = 1000000;
    }

    const bool operator>(const queueIntersect & rhs) const {
        return (travelTime > rhs.travelTime);
    }

    const bool operator<(const queueIntersect & rhs) const {
        //cout << "Operator < (less than): " << endl;
        return (travelTime < rhs.travelTime);
    }

    const bool operator==(const queueIntersect & rhs) const {
        return (intersectID == rhs.intersectID);
    }

};

struct intersectPath {
    unsigned int intersectionID;
    unsigned int reaching_segment;
    unsigned int prevIntersectionID;
    double timeToNode;
    int setAffinity;

    intersectPath(unsigned intersect, unsigned seg, unsigned prev, double travel, int affinity) {
        intersectionID = intersect;
        reaching_segment = seg;
        prevIntersectionID = prev;
        timeToNode = travel;
        setAffinity = affinity;
    }

    intersectPath() {
        intersectionID = -1;
        reaching_segment = -1;
        prevIntersectionID = -1;
        timeToNode = 1000000;
        setAffinity = 0;
    }

};


// /cad2/ece297s/public/m3/unit_tests$ - Unit Test Location

double compute_path_travel_time(const std::vector<unsigned>& path) {
    // Returns the time required to travel along the path specified. The path
    // is passed in as a vector of street segment ids, and this function can
    // assume the vector either forms a legal path or has size == 0.
    // The travel time is the sum of the length/speed-limit of each street
    // segment, plus 15 seconds per turn implied by the path. A turn occurs
    // when two consecutive street segments have different street names.

    double time = 0;
    for (unsigned int i = 0; i < path.size(); i++) {
        time += find_segment_travel_time(path[i]);

        if (i == 0)
            continue;

        if (getStreetSegmentStreetID(path[i]) != getStreetSegmentStreetID(path[i - 1])) {
            time += TURN_OFFSET;
        }
    }

    return time;

}

map<unsigned, intersectInfo> find_path_to_mult_intersections(unsigned int intersect_id_start, vector<unsigned> intersect_id_end) {

    vector<bool> found;
    int size = intersect_id_end.size();
    int count = 0;
    for (int i = 0; i < intersect_id_end.size(); i++) {
        found.push_back(false);
    }

    map<unsigned, intersectInfo> paths;

    //cout << "Start" << endl;

    int numberIntersects = getNumberOfIntersections();
    vector<intersectPath> pathMap(numberIntersects);

    multiset<queueIntersect> min_path_heap;

    pathMap[intersect_id_start] = intersectPath(intersect_id_start, NO_EDGE, NO_EDGE, NO_TIME, PENDING_SET);

    min_path_heap.insert(queueIntersect(intersect_id_start, NO_TIME));

    while (!min_path_heap.empty()) {
        //Get the node with the lowest travel time from the heap and then delete it from the vector
        multiset<queueIntersect>::iterator current = min_path_heap.begin();

        //Mark that this node has been visited
        pathMap [current->intersectID].setAffinity = CLOSED_SET;

        //Perform check for reaching destination
        for (int j = 0; j < size; j++) {
            if (current->intersectID == intersect_id_end[j]) {
                found[j] = true;
                count++;
            }
        }

        if (count == size) {
            break;
        }


        //Iterate through the street segments that are attached to the intersection
        vector<unsigned> segments = find_intersection_street_segments(current->intersectID);
        for (unsigned int i = 0; i < segments.size(); i++)
            //for (vector<unsigned>::iterator iter = segments.begin(); iter != segments.end(); iter++)
        {
            unsigned visitingID, segment;
            StreetSegmentEnds ends;
            double pathTime, heuristicTime;

            segment = segments[i];

            //Perform check for legality of street
            ends = getStreetSegmentEnds(segment);
            if ((current->intersectID == ends.to) && getStreetSegmentOneWay(segment)) {
                continue;
            }

            //Gets the appropriate intersection ID from the ends
            if (current->intersectID == ends.to) {
                visitingID = ends.from; //static_cast<unsigned int> (ends.from);
            } else {
                visitingID = ends.to; //static_cast<unsigned int>(ends.to);
            }

            //Performs check to see if the node has already been visited personally
            if (pathMap [visitingID].setAffinity == CLOSED_SET) {
                continue;
            }

            pathTime = pathMap[current->intersectID].timeToNode;
            pathTime += find_segment_travel_time(segment);
            //heuristicTime = find_distance_between_two_points(getIntersectionPosition(visitingID), getIntersectionPosition(intersect_id_end)) * 60 / (100000);
            //heuristicTime *= 1.5;

            //Implement a check to add offset for turns
            if (((pathMap[current->intersectID]).reaching_segment != NO_EDGE) && (getStreetSegmentStreetID(pathMap[current->intersectID].reaching_segment) != getStreetSegmentStreetID(segment))) {
                pathTime += TURN_OFFSET;
            }

            if ((pathMap[visitingID]).setAffinity == OPEN_SET) {
                pathMap[visitingID] = intersectPath(visitingID, segment, current->intersectID, pathTime, PENDING_SET);
                min_path_heap.insert(queueIntersect(visitingID, (pathTime /*+ heuristicTime*/)));
                continue;
            } else if (pathMap[visitingID].timeToNode <= pathTime) {
                continue;
            }

            pathMap[visitingID] = intersectPath(visitingID, segment, current->intersectID, pathTime, PENDING_SET);
            //multiset<queueIntersect>::iterator iter = min_path_heap.find(queueIntersect(visitingID, (pathTime + heuristicTime)));
            //min_path_heap.erase(iter);
            //min_path_heap.erase(queueIntersect(visitingID, (pathTime/* + heuristicTime*/)));
            min_path_heap.insert(queueIntersect(visitingID, (pathTime/* + heuristicTime*/)));
        }
        min_path_heap.erase(current);
    }

    //    if (count != size) {
    //
    //
    //        return path;
    //    }
    //
    for (int k = 0; k < size; k++) {

        vector<unsigned> path;
        intersectPath traverse = pathMap[intersect_id_end[k]];


        while (traverse.reaching_segment != NO_EDGE) {
            path.push_back(traverse.reaching_segment);
            //cout << traverse.reaching_segment << endl;
            traverse = pathMap[traverse.prevIntersectionID];
        }
        reverse(path.begin(), path.end());
        unsigned int id = intersect_id_end[k];
        double time = pathMap[id].timeToNode;
        paths[id] = intersectInfo(path, time);
    }

    //cout << "Done" << endl;

    return paths;

}

vector<unsigned> find_path_between_intersections(unsigned int intersect_id_start, unsigned int intersect_id_end) {
    bool found = false;
    vector <unsigned> path;

    if (intersect_id_start == intersect_id_end) {
        return path;
    }


    int numberIntersects = getNumberOfIntersections();
    vector<intersectPath> pathMap(numberIntersects);

    multiset<queueIntersect> min_path_heap;

    pathMap[intersect_id_start] = intersectPath(intersect_id_start, NO_EDGE, NO_EDGE, NO_TIME, PENDING_SET);

    min_path_heap.insert(queueIntersect(intersect_id_start, NO_TIME));

    while (!min_path_heap.empty()) {
        //Get the node with the lowest travel time from the heap and then delete it from the vector
        multiset<queueIntersect>::iterator current = min_path_heap.begin();

        //Mark that this node has been visited
        pathMap [current->intersectID].setAffinity = CLOSED_SET;

        //Perform check for reaching destination
        if (current->intersectID == intersect_id_end) {
            found = true;
            break;
        }

        //Iterate through the street segments that are attached to the intersection
        vector<unsigned> segments = find_intersection_street_segments(current->intersectID);
        for (unsigned int i = 0; i < segments.size(); i++)
            //for (vector<unsigned>::iterator iter = segments.begin(); iter != segments.end(); iter++)
        {
            unsigned visitingID, segment;
            StreetSegmentEnds ends;
            double pathTime, heuristicTime;

            segment = segments[i];

            //Perform check for legality of street
            ends = getStreetSegmentEnds(segment);
            if ((current->intersectID == ends.to) && getStreetSegmentOneWay(segment)) {
                continue;
            }

            //Gets the appropriate intersection ID from the ends
            if (current->intersectID == ends.to) {
                visitingID = ends.from; //static_cast<unsigned int> (ends.from);
            } else {
                visitingID = ends.to; //static_cast<unsigned int>(ends.to);
            }

            //Performs check to see if the node has already been visited personally
            if (pathMap [visitingID].setAffinity == CLOSED_SET) {
                continue;
            }

            pathTime = pathMap[current->intersectID].timeToNode;
            pathTime += find_segment_travel_time(segment);
            heuristicTime = find_distance_between_two_points(getIntersectionPosition(visitingID), getIntersectionPosition(intersect_id_end)) * 60 / (100000);
            //heuristicTime *= 1.5;

            //Implement a check to add offset for turns
            if (((pathMap[current->intersectID]).reaching_segment != NO_EDGE) && (getStreetSegmentStreetID(pathMap[current->intersectID].reaching_segment) != getStreetSegmentStreetID(segment))) {
                pathTime += TURN_OFFSET;
            }

            if ((pathMap[visitingID]).setAffinity == OPEN_SET) {
                pathMap[visitingID] = intersectPath(visitingID, segment, current->intersectID, pathTime, PENDING_SET);
                min_path_heap.insert(queueIntersect(visitingID, (pathTime + heuristicTime)));
                continue;
            } else if (pathMap[visitingID].timeToNode <= pathTime) {
                continue;
            }

            pathMap[visitingID] = intersectPath(visitingID, segment, current->intersectID, pathTime, PENDING_SET);
            //multiset<queueIntersect>::iterator iter = min_path_heap.find(queueIntersect(visitingID, (pathTime + heuristicTime)));
            //min_path_heap.erase(iter);
//            min_path_heap.erase(queueIntersect(visitingID, (pathTime + heuristicTime)));
            min_path_heap.insert(queueIntersect(visitingID, (pathTime + heuristicTime)));
        }
        min_path_heap.erase(current);
    }

    if (!found) {
        return path;
    }

    intersectPath traverse = pathMap[intersect_id_end];

    //int count = 1;
    while (traverse.reaching_segment != NO_EDGE) {
        path.push_back(traverse.reaching_segment);
        //cout << traverse.reaching_segment << endl;
        traverse = pathMap[traverse.prevIntersectionID];
        //count++;
    }

    //cout << count << endl;
    reverse(path.begin(), path.end());
    return path;

}


// /cad2/ece297s/public/m3/unit_tests$ - Unit Test Location

void drawPath(vector<unsigned> & path, unsigned int width, t_bound_box mapBounds) {
    setcolor(BLUE);
    setlinewidth(width);
    StreetSegmentEnds ends;
    unsigned int fromID, toID;
    LatLon from, to;
    unsigned int curveCount;
    float boundaryRight = mapBounds.right();
    float boundaryLeft = mapBounds.left();
    float boundaryTop = mapBounds.top();
    float boundaryBottom = mapBounds.bottom();

    int count = 1;


    // Iterate through all of the street segments
    // If a street segment is within the bounds of the visible world, it will be drawn. Else, it will be ignored
    // For no curve points, we can simply draw a straight line
    // For 1 curve point, draw to the curve point, then from it to the end
    // In all other cases, draw to the first one, between all existing pairs and then from the last one
    for (int Q = 0; Q < path.size(); Q++) {
        //cout << "Iteration: " << count << endl;
        count++;
        ends = getStreetSegmentEnds(path[Q]);
        from = getIntersectionPosition(ends.from);
        to = getIntersectionPosition(ends.to);
        if (((boundaryRight > from.lon) && (from.lon > boundaryLeft)) && (from.lat < boundaryTop) && (from.lat > boundaryBottom) || (boundaryRight > to.lon) && (to.lon > boundaryLeft) && (to.lat < boundaryTop) && (to.lat > boundaryBottom)) {

            curveCount = getStreetSegmentCurvePointCount(path[Q]);
            if (curveCount == 0) {
                fromID = ends.from;
                toID = ends.to;
                from = getIntersectionPosition(fromID);
                to = getIntersectionPosition(toID);
                drawline(from.lon, from.lat, to.lon, to.lat);
            } else if (curveCount == 1) {
                fromID = ends.from;
                from = getIntersectionPosition(fromID);
                to = getStreetSegmentCurvePoint(path[Q], 0);
                drawline(from.lon, from.lat, to.lon, to.lat);
                from = to;
                toID = ends.to;
                to = getIntersectionPosition(toID);
                drawline(from.lon, from.lat, to.lon, to.lat);
            } else {
                fromID = ends.from;
                from = getIntersectionPosition(fromID);
                to = getStreetSegmentCurvePoint(path[Q], 0);
                drawline(from.lon, from.lat, to.lon, to.lat);
                for (unsigned int i = 1; i < curveCount; i++) {
                    from = getStreetSegmentCurvePoint(path[Q], i);
                    to = getStreetSegmentCurvePoint(path[Q], i - 1);
                    drawline(from.lon, from.lat, to.lon, to.lat);
                }
                from = getStreetSegmentCurvePoint(path[Q], curveCount - 1);
                toID = ends.to;
                to = getIntersectionPosition(toID);
                drawline(from.lon, from.lat, to.lon, to.lat);
            }
        }
    }
}

// Returns the shortest travel time path (vector of street segments) from
// the start intersection to a point of interest with the specified name.
// If no such path exists, returns an empty (size == 0) vector.
vector<unsigned> find_path_to_point_of_interest(unsigned intersect_id_start, string point_of_interest_name) {
    LatLon position = getIntersectionPosition(intersect_id_start);
    double closest = DEFAULT_DISTANCE;
    double secondClosest = DEFAULT_DISTANCE;
    double checkDistance;
    unsigned int closestID, secondClosestID;

    cityBlock proximity;
    vector<unsigned> checklist, path, pathTwo;
    double timeOne, timeTwo;

    LatLon copy = position;
    LatLon origin;
    double xDistance, yDistance;

    vector<float> bounds = findCoordinates();
    origin.lon = bounds[0];
    origin.lat = bounds[1];

    copy.lon = bounds[0];
    yDistance = find_distance_between_two_points(origin, copy);
    copy.lat = bounds[1];
    copy.lon = position.lon;
    xDistance = find_distance_between_two_points(origin, copy);

    if ((xDistance <= parameters.getFirstLon()) && (yDistance >= parameters.getSecondLat())) {
        proximity = cityGrid[TOP_LEFT];
    } else if ((xDistance >= parameters.getSecondLon()) && (yDistance >= parameters.getSecondLat())) {
        proximity = cityGrid[TOP_RIGHT];
    } else if ((xDistance <= parameters.getFirstLon()) && (yDistance <= parameters.getFirstLat())) {
        proximity = cityGrid[BOTTOM_LEFT];
    } else if ((xDistance >= parameters.getSecondLon()) && (yDistance >= parameters.getFirstLat())) {
        proximity = cityGrid[BOTTOM_RIGHT];
    } else if ((xDistance <= parameters.getFirstLon()) && (yDistance <= parameters.getSecondLat()) && (yDistance >= parameters.getFirstLat())) {
        proximity = cityGrid[MID_LEFT];
    } else if ((xDistance >= parameters.getFirstLon()) && (yDistance >= parameters.getSecondLat()) && (xDistance <= parameters.getSecondLon())) {
        proximity = cityGrid[TOP_MID];
    } else if ((xDistance >= parameters.getSecondLon()) && (yDistance <= parameters.getSecondLat()) && (yDistance >= parameters.getFirstLat())) {
        proximity = cityGrid[MID_RIGHT];
    } else if ((xDistance >= parameters.getFirstLon()) && (yDistance <= parameters.getFirstLat()) && (xDistance <= parameters.getSecondLon())) {
        proximity = cityGrid[BOTTOM_MID];
    } else {
        proximity = cityGrid [CENTER];
    }

    checklist = proximity.getHere();
    vector<unsigned> near = proximity.getNear();
    vector<unsigned> far = proximity.getFar();
    unsigned int checkingID;
    string name;

    // Loop through the associated city block to see if we can find the POIs in proximity of the intersection
    for (unsigned int i = 0; i < checklist.size(); i++) {
        checkingID = checklist[i];
        name = getPointOfInterestName(checkingID);
        if (point_of_interest_name == name) {
            checkDistance = find_distance_between_two_points(position, getPointOfInterestPosition(checkingID));

            if (checkDistance < closest && checkDistance < secondClosest) {
                closestID = checkingID;
                closest = checkDistance;
                continue;
            }

            if (checkDistance > closest && checkDistance < secondClosest) {
                secondClosestID = checkingID;
                secondClosest = checkDistance;
                continue;
            }

        }

        // If we haven't found two POIs to compare, we must continue searching in adjacent blocks
        // We want to cycle through the near vector and visit each adjacent block separately
        //if(closest == DEFAULT_DISTANCE || secondClosest == DEFAULT_DISTANCE){
        int k = 0;

        while (k < near.size()) {

            proximity = cityGrid[near[k]];
            checklist = proximity.getHere();

            for (unsigned int z = 0; z < checklist.size(); z++) {
                checkingID = checklist[z];
                name = getPointOfInterestName(checkingID);
                if (point_of_interest_name == name) {
                    checkDistance = find_distance_between_two_points(position, getPointOfInterestPosition(checkingID));

                    if (checkDistance < closest && checkDistance < secondClosest) {
                        closestID = checkingID;
                        closest = checkDistance;
                        continue;
                    }

                    if (checkDistance > closest && checkDistance < secondClosest) {
                        secondClosestID = checkingID;
                        secondClosest = checkDistance;
                        continue;
                    }
                }
            }
            k++;
        }
        //}

        if (closest == DEFAULT_DISTANCE || secondClosest == DEFAULT_DISTANCE) {
            int y = 0;

            while (y < far.size()) {

                proximity = cityGrid[far[y]];
                checklist = proximity.getHere();

                for (unsigned int x = 0; x < checklist.size(); x++) {
                    checkingID = checklist[x];
                    name = getPointOfInterestName(checkingID);
                    if (point_of_interest_name == name) {
                        checkDistance = find_distance_between_two_points(position, getPointOfInterestPosition(checkingID));

                        if (checkDistance < closest && checkDistance < secondClosest) {
                            closestID = checkingID;
                            closest = checkDistance;
                            continue;
                        }

                        if (checkDistance > closest && checkDistance < secondClosest) {
                            secondClosestID = checkingID;
                            secondClosest = checkDistance;
                            continue;
                        }
                    }
                }
                y++;
            }
        }

    }

    LatLon posPls;

    if (closest == DEFAULT_DISTANCE && secondClosest == DEFAULT_DISTANCE) {
        return path;
    }

    if (secondClosest == DEFAULT_DISTANCE) {

        posPls = getPointOfInterestPosition(closestID);
        closestID = findClosestIntersection(posPls);
        return path = find_path_between_intersections(intersect_id_start, closestID);
    }

    posPls = getPointOfInterestPosition(closestID);
    closestID = findClosestIntersection(posPls);
    posPls = getPointOfInterestPosition(secondClosestID);
    secondClosestID = findClosestIntersection(posPls);

    path = find_path_between_intersections(intersect_id_start, closestID);
    pathTwo = find_path_between_intersections(intersect_id_start, secondClosestID);
    timeOne = compute_path_travel_time(path);
    timeTwo = compute_path_travel_time(pathTwo);

    if (timeOne <= timeTwo) {
        return path;
    } else {
        return pathTwo;
    }

    //    unsigned int pointOfInterestID = 0;
    //    vector<unsigned> matches;
    //    vector<LatLon> matchPos;
    //    unsigned int i;
    //
    //    for (i = 0; i < getNumberOfPointsOfInterest(); i++) {
    //        if (getPointOfInterestName(i) == point_of_interest_name) {
    //            matches.push_back(i);
    //            matchPos.push_back(getPointOfInterestPosition(i));
    //        }
    //    }
    //
    //    vector<double> dist;
    //    LatLon startPos = getIntersectionPosition(intersect_id_start);
    //
    //    for (i = 0; i < matches.size(); i++) {
    //       double heuristicTime = find_distance_between_two_points(startPos, matchPos[i]) * 60 / (100000);
    //       dist.push_back(heuristicTime);
    //    }
    //
    //    double minTime = 50000;
    //    vector<unsigned> closest;
    //	unsigned int index;
    //    unsigned int currClosest;
    //    double currDistance;
    //	for (i = 0; i < 3; i++) {
    //	    for (unsigned int j = 0; j < matches.size(); j++) {
    //	        if (i == 0) {
    //	            currClosest = matches[0];
    //	            currDistance = dist[0];
    //				index = 0;
    //	        } else {
    //	            if (dist[i] < currDistance) {
    //					currDistance = dist[i];
    //					currClosest = matches[i];
    //					index = i;
    //				}
    //		    }
    //		}
    //		closest.push_back(matches[index]);
    //		dist[index] = 1000000;
    //	}
    //
    //    vector<unsigned> path1 = find_path_between_intersections(intersect_id_start, closest[0]);
    //    double shortest = compute_path_travel_time(path1);
    //
    //    vector<unsigned> path2 = find_path_between_intersections(intersect_id_start, closest[1]);
    //    double shortest1 = compute_path_travel_time(path2);
    //
    //    vector<unsigned> path3 = find_path_between_intersections(intersect_id_start, closest[2]);
    //    double shortest2 = compute_path_travel_time(path3);
    //
    //    if((shortest > shortest1) && (shortest > shortest2)){
    //        return path1;
    //    }
    //    else if((shortest1 > shortest) && (shortest1 > shortest2)){
    //        return path2;
    //    }
    //    else if ((shortest2 > shortest) && (shortest2 > shortest1)){
    //        return path3;
    //    }
    //

    //return empty;
}

/*MITCH's FUNCTIONS - NOTED HERE FOR EASIER FOLDING*/
LatLon findCentre(void) {
    t_bound_box mapBounds = get_visible_world(); //find the bounds of the current map view

    LatLon bottomLeft, topRight;
    bottomLeft.lat = mapBounds.bottom();
    bottomLeft.lon = mapBounds.left();
    topRight.lat = mapBounds.top();
    topRight.lon = mapBounds.right();

    float avgLat = (bottomLeft.lat + topRight.lat) / 2; //Compute the average latitude of the view
    float avgLon = (bottomLeft.lon + topRight.lon) / 2; //Compute the average longitude of the view

    LatLon centre;
    centre.lat = avgLat;
    centre.lon = avgLon;

    return centre;
}

//Given an intersection id, does that intersection contain any street named unknown (i.e. id = 0)?

bool containsUnknown(unsigned id) {
    vector<unsigned> streetIDs = find_intersection_street_IDs(id);
    for (unsigned int i = 0; i < streetIDs.size(); i++) {
        if (!streetIDs[i]) {
            return true; //As soon as any streetID in the list is found to be 0, return true.
        }
    }
    return false; //If we survive to this point, then there are no unknown streets at this intersection.
}
/*
vector<unsigned> findClosestIntersectionsToPOI(unsigned POIid, unsigned n) {
    LatLon POIpos = getPointOfInterestPosition(POIid);
    return closestIntersectionList(POIpos, n);
}
*/

unsigned findClosestIntersection(LatLon here, vector<unsigned> exclude) {
    unsigned int id; //Once we find the closest intersection, its ID goes here
    LatLon position; //Holds the position of the intersection being tested.
    double distance; //The distance between the current intersection and POI
    double closest; //The currently closest intersection.
    bool already = false; //True if the current intersection is already in the list (i.e. exclude it).

    for (unsigned int i = 1; i < getNumberOfIntersections(); i++) { //We exclude id = 0 because those streets are unknown
        already = false; //Reset the exclude flag
        position = getIntersectionPosition(i);
        distance = find_distance_between_two_points(here, position);
        if (i == 1) { //If this is the first intersection, it is automatically the closest.
            closest = distance;
            id = 0;
        } else {
            for (unsigned int j = 0; j < exclude.size(); j++) {
                if (exclude[j] == i) {
                    already = true;
                }
            }
            if (distance < closest && !already && !containsUnknown(i)) { //This only happens if the intersection doesn't appear in exclude AND if the intersection doesn't contain any unknown streets.
                id = i; //USURPED!
                closest = distance; //Update the distance to the closest intersection
            }
        }
    }
    //When all of the dust has settled, return the ID of the closest intersection.
    return id;
}

//Use this function if we only want the geographically closest intersection.

unsigned findClosestIntersection(LatLon here) {
    unsigned int id; //Once we find the closest intersection, its ID goes here
    LatLon position; //Holds the position of the intersection being tested.
    double distance; //The distance between the current intersection and POI
    double closest; //The currently closest intersection.

    for (unsigned int i = 1; i < getNumberOfIntersections(); i++) { //We exclude id = 0 because those streets are unknown
        position = getIntersectionPosition(i);
        distance = find_distance_between_two_points(here, position);
        if (i == 1) { //If this is the first intersection, it is automatically the closest.
            closest = distance;
            id = 0;
        } else {
            if (distance < closest && !containsUnknown(i)) { //This only happens if the intersection doesn't appear in exclude.
                id = i; //USURPED!
                closest = distance; //Update the distance to the closest intersection
            }
        }
    }
    //When all of the dust has settled, return the ID of the closest intersection.
    return id;
}
/*
vector<unsigned> closestIntersectionList(LatLon here, unsigned int n) {
    vector<unsigned> list; //This will contain the n closest intersections.
    unsigned int id;
    for (unsigned int i = 0; i < n; i++) {
        id = findClosestIntersection(here, list);
        list.push_back(id);
    }
    return list;
}

	unsigned currentStreet, nextStreet; // street IDs
	string streetName, direction; // strings for the name and direction to be printed out
	StreetSegmentEnds segEnds2, segEnds3; // for finding the intersections
	// position coordinates - used to find left/ right turning
	LatLon intersectionOne = getIntersectionPosition(firstIntersection);
	LatLon intersectionTwo, intersectionThree;
	// check value for turn direction
	double checkTurn;
*/


/*FAN's FUNCTIONS - NOTED HERE FOR EASIER FOLDING*/
// Fan wrote this. It has yet to be tested, but Ivan will have to call it at some point to test it, right?
// This function will return a list of directions
// the code on determining left/ right is actually taken from the wikipedia article on graham scans
// this info can be found at http://en.wikipedia.org/wiki/Graham_scan
// another site that was looked at was for some complex hull algorithm which is actually pree cool
// the logic used here to determine the direction is http://www.geeksforgeeks.org/convex-hull-set-2-graham-scan/

void printDirections(vector <unsigned>& path, unsigned intersectionStart, unsigned intersectionEnd) {
    // instruction number
    unsigned instructionNum = 1;

    // all the intersection IDs
    unsigned firstIntersection = intersectionStart;
    unsigned secondIntersection, thirdIntersection;

    unsigned currentStreet, nextStreet;     // street IDs
    double legDistance = 0.0;               // Travel distance for each step of the route
    double legTime = 0.0;                   // Travel time for each step of the route
    int rLegDistance = 0;                   // Travel distance rounded to the nearest 100 m
    int rLegTime = 0.0;                     // Travel time rounded to the nearest minute
    bool km = false;                        // Should the distance be converted to and expressed in km?
    string streetName, direction;           // strings for the name and direction to be printed out
    StreetSegmentEnds segEnds2, segEnds3;   // for finding the intersections
    // position coordinates - used to find left/ right turning
    LatLon intersectionOne = getIntersectionPosition(firstIntersection);
    LatLon intersectionTwo, intersectionThree;
    // check value for turn direction
    double checkTurn;

    // loops through the vector of segment ids
    // if two segments are from different streets, then find turn direction
    // else continue thru loop
    cout << endl;
    cout << "Your travel directions for " << getIntersectionName(intersectionStart) << " to " << getIntersectionName(intersectionEnd) << "." << endl << endl
            << "Estimated travel time is " << compute_path_travel_time(path) << " minutes by car." << endl << endl;
    for (unsigned i = 0; i < path.size() - 1; i++) {
        currentStreet = getStreetSegmentStreetID(path[i]);  // get the street it belongs to
        nextStreet = getStreetSegmentStreetID(path[i + 1]); // get the street of the next segment to

        segEnds2 = getStreetSegmentEnds(path[i]);           // get the current segment's ends
        segEnds3 = getStreetSegmentEnds(path[i + 1]);       // get the next segment's ends

        //Add to the total distance/time so far on the current street
        legTime += find_segment_travel_time(path[i]);
        legDistance += find_street_segment_length(path[i]);

        // checking for the correct intersection ID
        // if the "to" value is the same as the previous intersection
        // change the next intersection to the "from" value
        // only matters really for the one-way streets
        secondIntersection = segEnds2.to; // want it to always be the following intersection, in case of 1-way
        if (secondIntersection == firstIntersection) {
            secondIntersection = segEnds2.from;
        }
        thirdIntersection = segEnds3.to; // want it to always be the following intersection, in case of 1-way
        if (thirdIntersection == secondIntersection) {
            thirdIntersection = segEnds3.from;
        }

        if (currentStreet != nextStreet) {
            //First, make the travel distance and time look nice.
            //Distance: Round to the nearest 0.1 km. If less than 1 km, express in m and round to the nearest 100 m.
            //Currently the distance is in metres.
            legDistance /= 100;                         //Now it is in hm
            rLegDistance = (int) legDistance;          //Cast it to an int (truncate decimals)
            double diff = legDistance - rLegDistance;  //Measure the error in truncation (always < 1)
            if (diff > 0.5) {
                rLegDistance++;                        //Round up (add 1 to the int) if 0.5 or greater
            }
            rLegDistance *= 100;                       //Convert it back to m

            double legDistance_km = 0.0;
            if (rLegDistance >= 1000) {                 //If greater than 1 km, express in km
                legDistance_km = (double) rLegDistance / 1000;
                km = true;
            }

            //Time: Round to the nearest min. If less than 1 min round up anyway.
            rLegTime = (int) legTime;                   //Cast the time to an int (truncate decimals)
            diff = legTime - rLegTime;                  //Measure the error in truncation (always < 1)
            if (diff > 0.5) {
                rLegTime++;                             //Round up (add 1 to the int) if 0.5 or greater
            }
            if (rLegTime == 0) {
                rLegTime = 1;                           //If less than 30 s, then round it up to 1 min
            }

            // only need to get positions if actually using
            intersectionTwo = getIntersectionPosition(secondIntersection);
            intersectionThree = getIntersectionPosition(thirdIntersection);
            // given in the form of an ordered triplet: (intersectionOne, intersectionTwo, intersectionThree)
            checkTurn = (intersectionTwo.lat - intersectionOne.lat) * (intersectionThree.lon - intersectionTwo.lon) -
                    (intersectionTwo.lon - intersectionOne.lon) * (intersectionThree.lat - intersectionTwo.lat);
            if (checkTurn == 0) {
                direction = "straight";
            } else if (checkTurn > 0) {
                direction = "right";
            } else {
                direction = "left";
            }
            // this algorithm was taken from wikipedia lmao
            streetName = getStreetName(nextStreet);
            if (!nextStreet) {
                streetName = "Sideroad";
            }
            int instrLength = 19; //Length of the instruction in chars. Starting with the constant parts, it is 4.
            if (direction == "right") {
                instrLength = 20; //"right" adds one char to the legnth.
            }
            instrLength += streetName.length();     //Add the length of the street name to the instruction length.
            if (instructionNum < 10) {
                cout << ' ';
            }
            cout << instructionNum << ". Turn " << direction << " onto " << streetName << ".";
            for (int i = 0; i < (75-instrLength); i++) {
                cout << ' ';
            }

            if (km) {
                cout << legDistance_km << " km, " << rLegTime << " min" << endl;
            } else {
                cout << rLegDistance << " m, " << rLegTime << " min" << endl;
            }

            instructionNum++;   // only increment if instruction actually made
            legDistance = 0.0;  //reset the leg distance to 0
            legTime = 0.0;      //reset the leg time to 0
            rLegDistance = 0;
            rLegTime = 0;
            km = false;
        }// if the segment's street == next segment's street, skip previous and loop

        // update the first intersection info
        firstIntersection = secondIntersection;
        intersectionOne = getIntersectionPosition(firstIntersection);
        // second and third intersection are set in function anyway, so won't set second to third here
    }
    cout << endl << "Have a nice trip!" << endl;
    cout << "\n" << endl; // adds extra space
}
