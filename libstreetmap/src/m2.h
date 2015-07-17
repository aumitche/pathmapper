#pragma once
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
using namespace std;
// Draws the map whose at map_path; this should be a .bin file.

extern vector<Street> sideroad;
extern vector<Street> highway;
extern vector<Street> mainRoad;
extern vector<Highway> filteredHighways;

extern vector<Feature> naturalClosed;
extern vector<Feature> naturalOpen;
extern vector<Feature> waterClosed;
extern vector<Feature> waterOpen;

void draw_map(string map_name, vector<Street> sideroad_, vector<Street> mainroad_, vector<Street> highway_);
void sortStreets (vector<Street> streets);
unsigned int find_num_segs(unsigned streetID);
//double averageSpeed (unsigned streetID);
vector<float> findCoordinates ();
void filterHighways(vector<Street> & highway);
LatLon findFirstInstanceOfHighway (string highwayName, vector<unsigned> & segs, t_bound_box mapBounds, double & angle, int count, int zoomLevel, vector<LatLon> & takenSpace, double distance);
void drawHighways (vector<Highway> & highways, t_bound_box mapBounds, int zoomLevel, double distance);
bool checkProximity(LatLon position, vector<LatLon> & takenPositions, int zoomLevel, double distance);

// finds all the coordinates in each feature and returns them in a vector of t_points
vector<t_point> featureCoordCount(unsigned featNum, unsigned featurePointCount);

// returns if the feature is open or closed
bool isFeatClosed (unsigned numPoints, vector<t_point> featureCoords);

// sorts the features according to attribute and if open/ closed
void sortFeatures();
// draws all natural features
void drawFeatures(int zoom);
// functions related to features
void drawOpenFeatures (vector<Feature> & features, int zoom);
void drawClosedFeatures (vector<Feature> & features);
vector<t_point> featureCoordCount(unsigned featNum, unsigned featurePointCount);
bool isFeatClosed (unsigned numPoints, vector<t_point> featureCoords);
void sortFeatures();
void drawOpenFeatures (vector<Feature> & features, int zoom);
void drawClosedFeatures (vector<Feature> & features);

void draw_screen(void);

// functions related to streets
unsigned int find_num_segs(unsigned streetID);

vector<float> findCoordinates ();
double findAngle(LatLon & from, LatLon & to);
LatLon findFirstInstance (string streetName, t_bound_box mapBounds, double & angle, int zoomLevel, vector<LatLon> & takenSpace, double distance);
void drawNames (vector<Street>  & streets, t_bound_box mapBounds, int zoomLevel, double distance);
void drawVector(vector<Street> & streets, int colour, unsigned int width, t_bound_box mapBounds);
//void drawHighways (vector<Street> & streets, t_bound_box mapBounds);

// all the POI sorting and drawing things
void POIon(void (*drawscreen) (void));
void sortPOI();
void drawPOI(double radius, int zoomLevel, t_bound_box mapBounds);

// Does route drawing things
void drawRouteEnds();
void drawRoute(unsigned int width, t_bound_box mapBounds);
void locationSelect (float x, float y, t_event_buttonPressed button_info);

// button and toggle things
void printHelp(void (*drawscreen) (void));

// drawing search things
void toggleSearch(void (*drawscreen) (void));
void drawSearch();

// drawing for TSP
void pathOn(void (*drawscreen) (void));
void depotOn(void (*drawscreen) (void));
void delivOn(void (*drawscreen) (void));
void pointSelect (float x, float y, t_event_buttonPressed button_info);
void drawPoints();

void toggleClear(void (*drawscreen) (void));
void drawTSP(unsigned int width, t_bound_box mapBounds);
void runTSP(void (*drawscreen) (void));
