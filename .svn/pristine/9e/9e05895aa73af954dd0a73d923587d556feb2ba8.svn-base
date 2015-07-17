#include <iostream>
#include <utility>
#include <cstdlib>
#include <vector>
#include <string>
#include <time.h>
#include <cmath>
#include <readline/readline.h>
#include <readline/history.h>
#include <boost/unordered_map.hpp>

#include <queue>
#include <set>
#include <algorithm>
#include <functional>

#include "easygl_constants.h"
#include "graphics.h"
#include "IntersectionGraph.h"
#include "StreetsDatabaseAPI.h"
#include "LatLon.h"

#include "Feature.h"
#include "Highway.h"
#include "POI.h"
#include "Street.h"
#include "intersectInfo.h"
#include "gridParam.h"
#include "cityBlock.h"

using namespace std;

// Returns a path (route) between the start intersection and the end
// intersection, if one exists. If no path exists, this routine returns
// an empty (size == 0) vector. If more than one path exists, the path
// with the shortest travel time is returned. The path is returned as a vector
// of street segment ids; traversing these street segments, in the given order,
// would take one from the start to the end intersection.
vector<unsigned> find_path_between_intersections(unsigned intersect_id_start, unsigned intersect_id_end);

map<unsigned, intersectInfo> find_path_to_mult_intersections(unsigned int intersect_id_start, vector<unsigned> intersect_id_end);

//vector<unsigned> bfsTraceBack(unsigned destID, vector<intersectPath> & pathMap);
//bool bfsFindPath (unsigned intersect_id_start, unsigned intersect_id_end, vector<intersectPath> & pathMap);
void drawPath(vector<unsigned> & path, unsigned int width, t_bound_box mapBounds);

// Returns the time required to travel along the path specified. The path
// is passed in as a vector of street segment ids, and this function can
// assume the vector either forms a legal path or has size == 0.
// The travel time is the sum of the length/speed-limit of each street
// segment, plus 15 seconds per turn implied by the path. A turn occurs
// when two consecutive street segments have different street names.
double compute_path_travel_time(const std::vector<unsigned>& path);

double estimated_travel_time(unsigned start_intersection_id, unsigned end_intersection_id);


// Returns the shortest travel time path (vector of street segments) from
// the start intersection to a point of interest with the specified name.
// If no such path exists, returns an empty (size == 0) vector.
vector<unsigned> find_path_to_point_of_interest (unsigned intersect_id_start, string point_of_interest_name);

unsigned findClosestIntersection (LatLon here, vector<unsigned> exclude);
unsigned findClosestIntersection (LatLon here);
LatLon findCentre(void);
vector<unsigned> findClosestIntersectionsToPOI (unsigned POIid, unsigned n);
vector<unsigned> closestIntersectionList(LatLon here, unsigned int n);

void printDirections(vector <unsigned> &path, unsigned intersectionStart, unsigned intersectionEnd);

void input(void);
void search(string query);
