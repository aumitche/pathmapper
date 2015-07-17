#include "m1.h"
#include "m2.h"
#include "m3.h"
#include "m4.h"
#include "time.h"

#define TIME_LIMIT 30
#define TNAUGHT 1000000
#define OPEN_SET 0
#define CLOSED_SET 1
#define PENDING_SET -1
#define NO_EDGE -1
#define NO_TIME 0
#define TURN_OFFSET 0.25

using namespace std;

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


intersectInfo performTwoOpt(intersectInfo & itinerary, unordered_map <unsigned, map<unsigned, intersectInfo>> &preComputedIntersectionPaths, unordered_map <unsigned, map<unsigned, intersectInfo>> &preComputedDepotPaths, vector<unsigned> depot_locations) {

    itinerary.path.pop_back();
    itinerary.path.erase(itinerary.path.begin());

    unsigned int firstSwappingID, secondSwappingID;

    firstSwappingID = rand() % itinerary.path.size();
    secondSwappingID = rand() % itinerary.path.size();

    while (firstSwappingID == secondSwappingID) {
        secondSwappingID = rand() % itinerary.path.size();
    }

    unsigned int bottomBound = min(firstSwappingID, secondSwappingID);
    unsigned int topBound = max(firstSwappingID, secondSwappingID);

    intersectInfo firstIteration, secondIteration, thirdIteration;

    reverse(itinerary.path.begin() + bottomBound, itinerary.path.begin() + topBound + 1);
    vector<unsigned> beginningRoute(itinerary.path.begin(), itinerary.path.begin() + bottomBound);
    vector<unsigned> middleRoute(itinerary.path.begin() + bottomBound, itinerary.path.begin() + topBound + 1);
    vector<unsigned> endRoute(itinerary.path.begin() + topBound + 1, itinerary.path.end());

    //    cout << "Begin: ";
    //    for(vector<unsigned>::iterator iter = beginningRoute.begin(); iter != beginningRoute.end(); iter++){
    //        cout << (*iter) << " " ;
    //    }
    //    cout << endl;
    //    
    //    cout << "Middle: ";
    //    for(vector<unsigned>::iterator iter = middleRoute.begin(); iter != middleRoute.end(); iter++){
    //        cout << (*iter) << " " ;
    //    }
    //    cout << endl;
    //    
    //    cout << "End: ";
    //    for(vector<unsigned>::iterator iter = endRoute.begin(); iter != endRoute.end(); iter++){
    //        cout << (*iter) << " " ;
    //    }
    //    cout << endl;
    //    
    //    


    computeTwoOpt(firstIteration, beginningRoute, middleRoute, endRoute, preComputedIntersectionPaths, preComputedDepotPaths, depot_locations);
    computeTwoOpt(secondIteration, middleRoute, beginningRoute, endRoute, preComputedIntersectionPaths, preComputedDepotPaths, depot_locations);
    computeTwoOpt(thirdIteration, endRoute, beginningRoute, middleRoute, preComputedIntersectionPaths, preComputedDepotPaths, depot_locations);

    intersectInfo optimalIteration = findBest(firstIteration, secondIteration, thirdIteration);

    return optimalIteration;
}

void computeTwoOpt(intersectInfo & computedPath, vector<unsigned> firstSection, vector<unsigned> secondSection, vector<unsigned> thirdSection, unordered_map <unsigned, map<unsigned, intersectInfo>> &preComputedIntersectionPaths, unordered_map <unsigned, map<unsigned, intersectInfo>> &preComputedDepotPaths, vector<unsigned> depot_locations) {

    computedPath.path.insert(computedPath.path.end(), firstSection.begin(), firstSection.end());
    computedPath.path.insert(computedPath.path.end(), secondSection.begin(), secondSection.end());
    computedPath.path.insert(computedPath.path.end(), thirdSection.begin(), thirdSection.end());

    unsigned int startPoint, endPoint;
    startPoint = computedPath.path[0];
    endPoint = computedPath.path[computedPath.path.size() - 1];

    // Compute the path time of only the intersections

    double totalTravelTime = 0;

    for (int i = 0; i < (computedPath.path.size() - 1); i++) {
        totalTravelTime += ((preComputedIntersectionPaths[computedPath.path[i]])[computedPath.path[i + 1]]).time;
    }

    unsigned int closestStart, closestEnd;
    double startTime, endTime;
    double closestStartTime = TNAUGHT;
    double closestEndTime = TNAUGHT;

    // Locate closest start depo and closest end depo
    for (int i = 0; i < depot_locations.size(); i++) {
        startTime = ((preComputedDepotPaths[depot_locations[i]])[startPoint]).time;
        endTime = ((preComputedIntersectionPaths[endPoint])[depot_locations[i]]).time;

        if (startTime <= closestStartTime) {
            closestStartTime = startTime;
            closestStart = depot_locations[i];
        }

        if (endTime <= closestEndTime) {
            closestEndTime = endTime;
            closestEnd = depot_locations[i];
        }
    }

    totalTravelTime += closestStartTime + closestEndTime;
    computedPath.path.insert(computedPath.path.begin(), closestStart);
    computedPath.path.push_back(closestEnd);
    computedPath.time = totalTravelTime;

}

intersectInfo findBest(intersectInfo & first, intersectInfo & second, intersectInfo & third) {
    double path = first.time;
    double pathTwo = second.time;
    double pathThree = third.time;

    if (path <= pathTwo && path <= pathThree) {
        return first;
    } else if (pathTwo <= path && pathTwo <= pathThree) {
        return second;
    } else if (pathThree <= path && pathThree <= pathTwo) {
        return third;
    }

    return first;


}

//This is the attempt at two opt

vector<unsigned> traveling_salesman(vector<unsigned> intersections_to_traverse, vector<unsigned> depot_locations) {
    //In simulated annealing, we start with the list of destinations in random order.
    //Might as well just use the vector we're given

    clock_t startTime = clock();
    cout << "Start Pre Compute" << endl;

    unsigned int sizeI = intersections_to_traverse.size();
    unsigned int sizeD = depot_locations.size();

    vector<vector<unsigned>> results;
    //vector<bool> visited(sizeI, false);

    unordered_map <unsigned, map<unsigned, intersectInfo>> preComputedDepotPaths, preComputedIntersectPaths;
    for (int N = 0; N < intersections_to_traverse.size(); N++) {
        vector<unsigned> intersectCopy(intersections_to_traverse);
        intersectCopy.erase(intersectCopy.begin() + N);
        intersectCopy.insert(intersectCopy.end(), depot_locations.begin(), depot_locations.end());
        map<unsigned, intersectInfo> pathI = find_path_to_mult_intersections(intersections_to_traverse[N], intersectCopy);
        preComputedIntersectPaths[intersections_to_traverse[N]] = pathI;
    }


    for (int M = 0; M < depot_locations.size(); M++) {
        map<unsigned, intersectInfo> pathD = find_path_to_mult_intersections(depot_locations[M], intersections_to_traverse);
        preComputedDepotPaths[depot_locations[M]] = pathD;
    }

    cout << "Done Pre Compute" << endl;

    double shortestDistance;
    double time = 0;
    vector<unsigned> path, totalPath;
    unsigned int next, trueFlag;

    if ((intersections_to_traverse.size() == 0) || (depot_locations.size() == 0)) { // check for nothing
        return totalPath;
    }

    if (sizeD == 1 && sizeI == 1) {
        return (find_path_between_intersections(depot_locations[0], intersections_to_traverse[0]));
    }


    unsigned int iterate = 1;
    unsigned current = intersections_to_traverse[0];
    next = intersections_to_traverse[1];

    int countAll = 0;
    //visited[0] = true;
    //path.push_back(depot_locations[0]);
    //path.push_back(current);

    //We want multiple starts from every depot location
    while (countAll < sizeD) {
        // Until we have started from every depot, start from the next one
        current = depot_locations[countAll];
        path.push_back(depot_locations[countAll]);
        next = intersections_to_traverse[0];
        vector<bool> visited(sizeI, false);
        bool done = false;

        // First node is the depot we are starting from
        // Find what the closest one should be
        while (!done) {


            shortestDistance = 10000000;
            unsigned int count = 0;

            // Check every intersection that we can
            for (unsigned int i = 0; i < sizeI; i++) {

                //If we've visited the intersection, we don't want to go back
                if (visited[i]) {
                    count++;
                    continue;
                }

                if (i == 0) {
                    shortestDistance = find_distance_between_two_points(getIntersectionPosition(current), getIntersectionPosition(intersections_to_traverse[i]));
                    next = intersections_to_traverse[i];
                    trueFlag = i;
                }

                if (shortestDistance < find_distance_between_two_points(getIntersectionPosition(current), getIntersectionPosition(intersections_to_traverse[i]))) {
                    continue;
                }

                shortestDistance = find_distance_between_two_points(getIntersectionPosition(current), getIntersectionPosition(intersections_to_traverse[i]));
                next = intersections_to_traverse[i];
                trueFlag = i;
            }

            if (count == visited.size()) {
                done = true;
            }
            // Mark the intersection we end up visiting
            // Add the path required to the total path
            // Set the current node to the next node
            visited[trueFlag] = true;

            //        if(count == 0){
            //            path.push_back(current);
            //            path.push_back(next);
            //        }
            path.push_back(next);

            //totalPath.insert(totalPath.end(), path.begin(), path.end());
            current = next;

            // If all visited we are done with intersections
            iterate++;

        }
        //Remove the excess address
        path.pop_back();

        shortestDistance = TNAUGHT;
        //Find the closest distance to the end intersection
        for (int i = 0; i < sizeD; i++) {
            if (i == 0) {
                shortestDistance = find_distance_between_two_points(getIntersectionPosition(depot_locations[i]), getIntersectionPosition(next));
                current = depot_locations[i];
                continue;
            }

            if (shortestDistance < find_distance_between_two_points(getIntersectionPosition(depot_locations[i]), getIntersectionPosition(next))) {
                continue;
            }

            shortestDistance = find_distance_between_two_points(getIntersectionPosition(depot_locations[i]), getIntersectionPosition(next));
            current = depot_locations[i];
        }
        path.push_back(current); // Push back the closest depo
        results.push_back(path);

        // Clear the vector for the next iteration
        path.erase(path.begin(), path.end());
        countAll++;
    }

    double shortestTime = TNAUGHT;
    unsigned int shortest;
    for (int j = 0; j < (results.size()); j++) {
        //Start at the shortest found path
        path = results[j];
        time = 0;
        for (int i = 0; i < (path.size() - 1); i++) {
            if (i == 0) {
                time += preComputedDepotPaths[path[i]][path[i + 1]].time;
            } else {
                time += preComputedIntersectPaths[path[i]][path[i + 1]].time;
            }
        }

        // If better path, save its time and where to find it
        if (time < shortestTime) {
            shortestTime = time;
            shortest = j;
        }

        path.erase(path.begin(), path.end());
    }

    path = results[shortest];
    intersectInfo itinerary(path, shortestTime);
    float timeSecs;
    int count = 0;
    // Loop while time has not expired
    do {
        //        count++;
        itinerary = performTwoOpt(itinerary, preComputedIntersectPaths, preComputedDepotPaths, depot_locations);
        if (itinerary.time <= shortestTime) {
            path = itinerary.path;
        }
        clock_t currentTime = clock();
        timeSecs = ((float) (currentTime - startTime)) / CLOCKS_PER_SEC;
    } while (timeSecs < TIME_LIMIT * 0.9 /*|| count < 20*/);

    for (int i = 0; i < (/*itinerary.*/path.size() - 1); i++) {

        if (i == 0) {
            totalPath.insert(totalPath.end(), preComputedDepotPaths[path[i]][path[i + 1]].path.begin(), preComputedDepotPaths[path[i]][path[i + 1]].path.end());
        } else {
            totalPath.insert(totalPath.end(), preComputedIntersectPaths[path[i]][path[i + 1]].path.begin(), preComputedIntersectPaths[path[i]][path[i + 1]].path.end());
        }
    }
    printDirections(totalPath, totalPath[0], totalPath[totalPath.size() - 1]);
    return totalPath;
}
