#define SIDEROADLENGTH 2000 // length is in metres
#define MAINROADLENGTH 5000
#define PI 3.1415926535897932384626433

//Zoom levels are based on the length of the diagonal of the screen, measured in metres.
#define ZOOM0 50000
#define ZOOM1 40000
#define ZOOM2 30000
#define ZOOM3 20000
#define ZOOM4 10000
#define ZOOM5 7500
#define ZOOM6 5000
#define ZOOM7 2500
#define ZOOM8 1500
#define ZOOM9 1000

#include "m1.h"
#include "m2.h"
#include "m3.h"
#include "m4.h"

using namespace std;

/*----------------------Global Variable things---------------------------*/

    // vectors for different streets
        vector<Street> sideroad;
        vector<Street> highway;
        vector<Street> mainroad;
        vector<Highway> filteredHighways;
        vector<LatLon> takenSpaceStreet;
    //vector<LatLon> takenSpacePOI;

    // vectors to hold the different types of natural features
    // note that "land" and "natural" are both lumped into the same vector, as they'll be drawn green
        vector<Feature> naturalClosed;
        vector<Feature> naturalOpen;
        vector<Feature> waterClosed;
        vector<Feature> waterOpen;
        vector<Feature> landClosed;

    // vector of points of interest
        vector<POI> poiVector;

    // variables for streets
        double avgSpeed;
        string name;
        double length;
        double numSegs;

    // variables for features
        unsigned featNum, numPoints;
        vector<t_point> featCoords;
        bool isClosed, isWater, isLand;
        string featAttr;

    // variables for POI
        bool drawFlag = false;

    // route drawing things
        unsigned locationNum = 0;
        LatLon intersectionArray[2];
        unsigned intersectionIDArray[2];
        bool printAlready = false;
        vector<unsigned> path;

    // TSP things
        bool drawDepotFlag = false;
        bool drawDropoffFlag = false;
        bool drawPathFlag = false;
        bool alreadyAdded = false;

        vector<unsigned> depot;
        vector<unsigned> dropoff;

    // variables for search bar
        bool drawSearchFlag = false;
        bool selectMessage;
        string displayName;

        bool clearPoints = false;
        bool drawTSPalready = false;
        vector<unsigned> TSPpath;
        bool drawTSPflag = false;


// Draws the map whose at map_path; this should be a .bin file.
void draw_map(string map_name, vector<Street> sideroad_, vector<Street> mainroad_, vector<Street> highway_) {
    sideroad = sideroad_;
    mainroad = mainroad_;
    highway = highway_;
    init_graphics("PathMapper 3.0", LIGHTGREY);
    create_button("Window", "HELP", printHelp);
    create_button("HELP", "Show Search", toggleSearch);
    create_button("HELP", "Turn POI On", POIon);
    create_button("HELP", "Path", pathOn);
    create_button("Path", "Depot", depotOn);
    create_button("Depot", "Dropoff", delivOn);
    create_button("Dropoff", "Clear Points", toggleClear);
    create_button("Dropoff", "Find TSP", runTSP);
    create_button("Window", "---", NULL);
    vector<float> bounds = findCoordinates();
    set_visible_world(bounds[0], bounds[1], bounds[2], bounds[3]);
    update_message("");
    event_loop(pointSelect, NULL, NULL, draw_screen);

    // this creates a stupid fucking button
    //create_button("Proceed", "---", NULL);

    close_graphics();
}

/*-----------------------------Streets-----------------------------------*/

    // Sorts the streets
    // Vector is previously created in the load_map function in m1.cpp
    // Takes all the streets, get the information on id, length, segments and average speed and creates an object of type Street
    // Sorts the streets into different vectors for layers (to be used at different zoom levels)

    // Finds the number of street segments that are associated with a given street
    unsigned int find_num_segs(unsigned streetID) {
        vector<unsigned> segs = find_street_street_segments_from_id(streetID);
        return segs.size();
    }

    // Given an intersection graph or map, returns the extremities of lat/lon coordinates for any intersection
    vector<float> findCoordinates(void) {
        vector<float> coordinates;
        float maxLon = 0;
        float minLon = 0;
        float maxLat = 0;
        float minLat = 0;
        LatLon checkingCoordinates;

        // For every intersection, check if its lat/lon points are the current max or min
        // If they are neither, simply iterate to the next intersection
        for (unsigned int i = 0; i < getNumberOfIntersections(); i++) {
            checkingCoordinates = getIntersectionPosition(i);

        // initializes the first set of coordinates to the first point
        if (i == 0){
          maxLon = checkingCoordinates.lon;
          minLon = checkingCoordinates.lon;
          maxLat = checkingCoordinates.lat;
          minLat = checkingCoordinates.lat;
            }

        if (checkingCoordinates.lon > maxLon) {
                maxLon = checkingCoordinates.lon;
            } else if (checkingCoordinates.lon < minLon) {
                minLon = checkingCoordinates.lon;
            }

            if (checkingCoordinates.lat > maxLat) {
                maxLat = checkingCoordinates.lat;
            } else if (checkingCoordinates.lat < minLat) {
                minLat = checkingCoordinates.lat;
            }
        }

        for (unsigned int i = 0; i < getFeatureCount(); i++){
            for (unsigned int j = 0; j < getFeaturePointCount(i); j++){
                checkingCoordinates = getFeaturePoint(i, j);
                if (checkingCoordinates.lon > maxLon) {
                    maxLon = checkingCoordinates.lon;
                }
                if (checkingCoordinates.lon < minLon) {
                    minLon = checkingCoordinates.lon;
                }

                if (checkingCoordinates.lat > maxLat) {
                    maxLat = checkingCoordinates.lat;
                }
                if (checkingCoordinates.lat < minLat) {
                    minLat = checkingCoordinates.lat;
                }
            }
        }
        //Place the variables in a vector in the order they would be called to generate the initial box
        coordinates.push_back(minLon);
        coordinates.push_back(minLat);
        coordinates.push_back(maxLon);
        coordinates.push_back(maxLat);

        return coordinates;
    }

/*------------------------Natural Features-------------------------------*/

    // find coordinates
    // returns a vector of all the coordinates of a feature given the number of the feature
    vector<t_point> featureCoordCount(unsigned featNum, unsigned featurePointCount) {
        // create a vector of t_point coordinates
        vector<t_point> featureCoords;
        t_point coord; // creates a t_point called coord
        // gets all the points of the feature and puts it into the vector
        for (unsigned i = 0; i < featurePointCount; i++) {
            // gets the lat and lon and sticks it into the x, y coordinate as a t_point
            coord.x = getFeaturePoint(featNum, i).lon;
            coord.y = getFeaturePoint(featNum, i).lat;
            // adds coord to the vector of coordinates
            featureCoords.push_back(coord);
        }
        return featureCoords;
    }

    // if first and last coordinates are same, then closed = true
    bool isFeatClosed (unsigned numPoints, vector<t_point> featureCoords) {
        t_point firstCoord = featureCoords[0];
        t_point lastCoord = featureCoords[numPoints-1];
        bool isClosed = true;
        // not considered closed if coordinates are different or there is only one coordinate
        if ((firstCoord.x != lastCoord.x) || (firstCoord.y != lastCoord.y) || (featureCoords.size() == 1)){
            isClosed = false;
        }
        return isClosed;
    }

    // Go through all of the features and create a new class holding all of its info
    // Then sort into water or land

    void sortFeatures() {
        unsigned totalFeatures = getFeatureCount();
        for (unsigned i = 0; i < totalFeatures; i++){
            featNum = i; // index of the feature in overall list
            numPoints = getFeaturePointCount(i); // number of coordinates
            featCoords = featureCoordCount(i, numPoints); // the actual coordinates, returned as a vector of t_points
            isClosed = isFeatClosed(numPoints, featCoords); // checks if closed or open feature

            // for feature attribute and type, check if of type "water" first
            isWater = true; // initially set to true
            featAttr = getFeatureAttribute(i, "water"); // check if water
            if (featAttr == "") { // Not water
                featAttr = getFeatureAttribute(i, "waterway"); // check if waterway

                if (featAttr == "") { // Not waterway
                    isWater = false; // if returns empty string, then not water type

                    // the rest test if it's natural or land
                    featAttr = getFeatureAttribute(i, "natural"); // check if natural
                    // for some reason, "water" is under "natural" and not "water"/ "waterway"
                    if (featAttr == "water") {
                        isWater = true;
                    } // Is water (what the fuck)
                    else if (featAttr == "") { // if it's nothing, just assume land lmfao
                        featAttr = getFeatureAttribute(i, "land"); // check if natural
                        if (featAttr == "island"){
                            isLand = true;
                        }
                    } // Not natural
                } // Not waterway
            } // Not water

            // throw all the above information into the constructor of Feature objects
            Feature feature(featNum, isWater, isClosed, featAttr, numPoints, featCoords);

            // sort them into five vectors, based on open/ closed and water/ non-water and land lmao
            // land
            if (isLand && isClosed && !isWater) {
                landClosed.push_back(feature);
            }
            // open, non-water
            else if (!isLand && !isWater && !isClosed){
                naturalOpen.push_back(feature);
            }
            // closed, non-water
            else if (!isLand && !isWater && isClosed){
                naturalClosed.push_back(feature);
            }
            // closed, water
            else if (!isLand && isWater && isClosed){
                waterClosed.push_back(feature);
            }
            // open, water
            else if (!isLand && isWater && !isClosed){
                waterOpen.push_back(feature);
            }
        }
    }

/*---------------------------Draws things--------------------------------*/
    void draw_screen(void) {
        clearscreen();
        setfontsize(10);
        vector<float> bounds = findCoordinates();
        t_bound_box mapBounds = get_visible_world();
        double distance; //A measure of the diagonal of the screen, in metres.

        vector<float> testBounds = findCoordinates();

        LatLon bottomLeft, topRight;
        bottomLeft.lat = mapBounds.bottom();
        bottomLeft.lon = mapBounds.left();
        topRight.lat = mapBounds.top();
        topRight.lon = mapBounds.right();


        distance = find_distance_between_two_points(bottomLeft, topRight);
        drawFeatures(distance);

        if (distance > ZOOM0){
            // update_message("Zoom level 0");
            setcolor(MAPBASE);
            fillrect(testBounds[0], testBounds[1], testBounds[2], testBounds[3]);
            drawFeatures(distance);
            drawVector(mainroad, MAINROAD, 2, mapBounds);
            drawVector(highway, HIGHWAY, 5, mapBounds);
            drawRoute(5, mapBounds);
            drawTSP(5, mapBounds);
        } else if (distance > ZOOM1) {
            // update_message("Zoom level 1");
            setcolor(MAPBASE);
            fillrect(testBounds[0], testBounds[1], testBounds[2], testBounds[3]);
            drawFeatures(distance);
            drawVector(mainroad, MAINROAD, 3, mapBounds);
            drawVector(highway, HIGHWAY, 6, mapBounds);
            drawRoute(6, mapBounds);
            drawTSP(6, mapBounds);
        } else if (distance > ZOOM2) {
            // update_message("Zoom level 2");
            setcolor(MAPBASE);
            fillrect(testBounds[0], testBounds[1], testBounds[2], testBounds[3]);
            drawFeatures(distance);
            drawVector(mainroad, MAINROAD, 5, mapBounds);
            drawVector(highway, HIGHWAY, 8, mapBounds);
            drawHighways(filteredHighways, mapBounds, ZOOM2, distance);
            takenSpaceStreet.clear();
            drawRoute(8, mapBounds);
            drawTSP(8, mapBounds);
        } else if (distance > ZOOM3) {
            // update_message("Zoom level 3");
            setcolor(MAPBASE);
            fillrect(testBounds[0], testBounds[1], testBounds[2], testBounds[3]);
            drawFeatures(distance);
            drawVector(sideroad, BACKGROUND, 1, mapBounds);
            drawVector(mainroad, MAINROAD, 5, mapBounds);
            drawVector(highway, HIGHWAY, 8, mapBounds);
            drawRoute(8, mapBounds);
            drawTSP(8, mapBounds);
        } else if (distance > ZOOM4) {
            // update_message("Zoom level 4");
            setcolor(MAPBASE);
            fillrect(testBounds[0], testBounds[1], testBounds[2], testBounds[3]);
            drawFeatures(distance);
            drawVector(sideroad, BACKGROUND, 1, mapBounds);
            drawVector(mainroad, MAINROAD, 7, mapBounds);
            drawVector(highway, HIGHWAY, 10, mapBounds);
            drawNames(mainroad, mapBounds, ZOOM4, distance);
            drawRoute(10, mapBounds);
            drawTSP(10, mapBounds);
            takenSpaceStreet.clear();
        } else if (distance > ZOOM5) {
            // update_message("Zoom level 5");
            setcolor(MAPBASE);
            fillrect(testBounds[0], testBounds[1], testBounds[2], testBounds[3]);
            drawFeatures(distance);
            drawVector(sideroad, BACKGROUND, 3, mapBounds);
            drawVector(mainroad, MAINROAD, 7, mapBounds);
            drawVector(highway, HIGHWAY, 10, mapBounds);
            drawNames(mainroad, mapBounds, ZOOM5, distance);
            drawPOI(distance, ZOOM5, mapBounds);
            takenSpaceStreet.clear();
            //takenSpacePOI.clear();
            drawRoute(10, mapBounds);
            drawTSP(10, mapBounds);
        } else if (distance > ZOOM6) {
            // update_message("Zoom level 6");
            setcolor(MAPBASE);
            fillrect(testBounds[0], testBounds[1], testBounds[2], testBounds[3]);
            drawFeatures(distance);
            drawVector(sideroad, BACKGROUND, 4, mapBounds);
            drawVector(mainroad, MAINROAD, 8, mapBounds);
            drawVector(highway, HIGHWAY, 11, mapBounds);
            drawNames(mainroad, mapBounds, ZOOM6, distance);
            drawPOI(distance, ZOOM6, mapBounds);
            takenSpaceStreet.clear();
            //takenSpacePOI.clear();
            drawRoute(11, mapBounds);
            drawTSP(11, mapBounds);
        } else if (distance > ZOOM7) {
            // update_message("Zoom level 7");
            setcolor(MAPBASE);
            fillrect(testBounds[0], testBounds[1], testBounds[2], testBounds[3]);
            drawFeatures(distance);
            drawVector(sideroad, BACKGROUND, 5, mapBounds);
            drawVector(mainroad, MAINROAD, 9, mapBounds);
            drawVector(highway, HIGHWAY, 12, mapBounds);
            drawNames(mainroad, mapBounds,ZOOM7, distance);
            drawNames(sideroad, mapBounds, ZOOM7, distance);
            drawPOI(distance, ZOOM7, mapBounds);
            takenSpaceStreet.clear();
            //takenSpacePOI.clear();
            drawRoute(12, mapBounds);
            drawTSP(12, mapBounds);
        } else if (distance > ZOOM8) {
            // update_message("Zoom level 8");
            setcolor(MAPBASE);
            fillrect(testBounds[0], testBounds[1], testBounds[2], testBounds[3]);
            drawFeatures(distance);
            drawVector(sideroad, BACKGROUND, 7, mapBounds);
            drawVector(mainroad, MAINROAD, 10, mapBounds);
            drawVector(highway, HIGHWAY, 13, mapBounds);
            drawNames(mainroad, mapBounds, ZOOM8, distance);
            drawNames(sideroad, mapBounds, ZOOM8, distance);
            drawPOI(distance, ZOOM8, mapBounds);
            takenSpaceStreet.clear();
            //takenSpacePOI.clear();
            drawRoute(13, mapBounds);
            drawTSP(13, mapBounds);
        } else {
            //update_message("Zoom level 9+");
            setcolor(MAPBASE);
            fillrect(testBounds[0], testBounds[1], testBounds[2], testBounds[3]);
            drawFeatures(distance);
            drawVector(sideroad, BACKGROUND, 8, mapBounds);
            drawVector(mainroad, MAINROAD, 10, mapBounds);
            drawVector(highway, HIGHWAY, 13, mapBounds);
            drawNames(mainroad, mapBounds, ZOOM9, distance);
            drawNames(sideroad, mapBounds, ZOOM9, distance);
            drawPOI(distance, ZOOM9, mapBounds);
            takenSpaceStreet.clear();
            //takenSpacePOI.clear();
            drawRoute(13, mapBounds);
            drawTSP(13, mapBounds);
        }
        drawSearch(); ////// fan fan fan
        drawPoints();
        LatLon centre = findCentre();
    //    update_message("(" + centre.lat + "," + centre.lon + ")");
        //cout << "Centred at (" << centre.lat << "," << centre.lon << ")" << endl;
        string message;
        string selectMode;
        if (selectMessage){
            if (drawPathFlag){
                selectMode = "Currently in: path selection mode. You have selected: ";
            }
            else if (drawDepotFlag){
                selectMode = "Currently in: depot selection mode. You have selected: ";
            }
            else if (drawDropoffFlag){
                selectMode = "Currently in: dropoff selection mode. You have selected: ";
            }
            message = selectMode + displayName;
            selectMessage = !selectMessage;
        }
        else {
            unsigned closestIntersection = findClosestIntersection(centre);
            message = "Centred near " + getIntersectionName(closestIntersection);
        }
        update_message(message);
    }

    // Uses the vector of Streets to draw all the streets that should appear on the map
    void drawVector (vector<Street>& streets, int colour, unsigned int width, t_bound_box mapBounds) {
        // Declare variables needed and set the colour and width to the appropriate values
        setcolor(colour);
        setlinewidth(width);
        StreetSegmentEnds ends;
        unsigned int fromID, toID;
        LatLon from, to;
        unsigned int curveCount;
        float boundaryRight = mapBounds.right();
        float boundaryLeft = mapBounds.left();
        float boundaryTop = mapBounds.top();
        float boundaryBottom = mapBounds.bottom();

        // Begins looping through the street segments to see where and if it can print them.
        for (unsigned int i = 0; i < streets.size(); i++) {
            vector<unsigned> segs = find_street_street_segments((streets[i].getName()));

            // Iterate through all of the street segments for the street being considered
            // If a street segment is within the bounds of the visible world, it will be drawn. Else, it will be ignored
            // For no curve points, we can simply draw a straight line
            // For 1 curve point, draw to the curve point, then from it to the end
            // In all other cases, draw to the first one, between all existing pairs and then from the last one
            for(vector<unsigned>::iterator iter = segs.begin(); iter < segs.end(); iter++) {
                ends  = getStreetSegmentEnds(*iter);
                from = getIntersectionPosition(ends.from);
                to = getIntersectionPosition(ends.to);
                if (((( boundaryRight > from.lon) && (from.lon > boundaryLeft)) && (from.lat < boundaryTop) && (from.lat > boundaryBottom)) || ((boundaryRight > to.lon) && (to.lon > boundaryLeft) &&  (to.lat < boundaryTop) && (to.lat > boundaryBottom))) {
                    curveCount = getStreetSegmentCurvePointCount(*iter);
                    if(curveCount == 0) {
                        fromID = ends.from;
                        toID = ends.to;
                        from = getIntersectionPosition(fromID);
                        to = getIntersectionPosition(toID);
                        drawline(from.lon, from.lat, to.lon, to.lat);
                    } else if (curveCount == 1) {
                        fromID = ends.from;
                        from = getIntersectionPosition(fromID);
                        to = getStreetSegmentCurvePoint(*iter, 0);
                        drawline(from.lon, from.lat, to.lon, to.lat);

                        from = to;
                        toID = ends.to;
                        to = getIntersectionPosition(toID);
                        drawline(from.lon, from.lat, to.lon, to.lat);
                    } else {
                        fromID = ends.from;
                        from = getIntersectionPosition(fromID);
                        to = getStreetSegmentCurvePoint(*iter, 0);
                        drawline(from.lon, from.lat, to.lon, to.lat);
                        for(unsigned int i = 1; i < curveCount; i++) {
                            from = getStreetSegmentCurvePoint(*iter, i);
                            to = getStreetSegmentCurvePoint(*iter, i-1);
                            drawline(from.lon, from.lat, to.lon, to.lat);
                        }
                        from = getStreetSegmentCurvePoint(*iter, curveCount-1);
                        toID = ends.to;
                        to = getIntersectionPosition(toID);
                        drawline(from.lon, from.lat, to.lon, to.lat);
                    }
                }
            }
        }
    }

    // Finds the angle at which the name of a street is to be written at
    double findAngle (LatLon & from, LatLon & to){
        double deltaX;
        double deltaY;
        double angle;
        deltaX = to.lon - from.lon;
        deltaY = to.lat - from.lat;

        angle = atan(deltaY/deltaX);
        angle = angle * (180/PI);

        // Since atan returns a value from -90 to 90 degrees and draw functions only accept values from
        // 0 to 360 degrees, we must rectify negative values.
        if (angle < 0){
            angle += 360;
            return angle;
        }
        return angle;
    }

    // draws all the features
    void drawFeatures(int zoom){
        // changes colour to green to draw non-water
        setcolor(NATURE);

        // draws all the closed natural
        drawClosedFeatures (naturalClosed);

        // draws all the open natural
        drawOpenFeatures (naturalOpen, zoom);

        setcolor(WATER);
        // draws all the open water
        drawOpenFeatures (waterOpen, zoom);

        // draws all the closed water by looping through the waterClosed vector
        drawClosedFeatures (waterClosed);

        // changes colour to white for islands
        setcolor(NATURE);
        drawClosedFeatures (landClosed);
    }

    // draws all the open features
    // use the line segment draw function
    void drawOpenFeatures (vector<Feature> & features, int zoom){
        string type;

        // need a vector of coordinates
        vector<t_point> vectorPoints;

        // loops through the vector of features
        for (unsigned i = 0; i < features.size(); i++){
            //First find out the zoom level and feature type to set line width
            type = getFeatureAttribute(i, "waterway");
            if (zoom >= ZOOM0) {
                setlinewidth(1);
            } else if (zoom >= ZOOM2) {
                if (type == "river") {
                    setlinewidth(1);
                } else if (type == "stream") {
                    setlinewidth(1);
                }
            } else if (zoom >= ZOOM4) {
                if (type == "river") {
                    setlinewidth(6);
                } else if (type == "stream") {
                    setlinewidth(5);
                } else if (type == "boatyard") {
                    setlinewidth(3);
                }
            } else if (zoom >= ZOOM6) {
                if (type == "river") {
                    setlinewidth(9);
                } else if (type == "stream") {
                    setlinewidth(7);
                } else if (type == "boatyard") {
                    setlinewidth(5);
                }
            } else if (zoom >= ZOOM8) {
                if (type == "river") {
                    setlinewidth(12);
                } else if (type == "stream") {
                    setlinewidth(10);
                } else if (type == "boatyard") {
                    setlinewidth(5);
                }
            }

            // gets the vector of coordinates for current feature
            vectorPoints = features[i].getCoords();
            if (vectorPoints.size() < 2) {
                fillrect(vectorPoints[0].x, vectorPoints[0].y, vectorPoints[0].x+0.0000000001, vectorPoints[0].x+0.0000000001);
            }
            // loops through the coordinates of the feature
            for (unsigned j = 0; j < vectorPoints.size()-1; j++){
                // passes in the individual x, y coords for each t_point
                //cout << vectorPoints[j].x << ", " << vectorPoints[j].y << endl;
                drawline(vectorPoints[j].x, vectorPoints[j].y, vectorPoints[j+1].x, vectorPoints[j+1].y);
            }
        }
    }

    void drawClosedFeatures (vector<Feature> & features){
        // need a vector of coordinates and the size
        vector<t_point> vectorPoints;
        unsigned vectorSize;

        // draws all the closed features by looping through the vector
        // use the polygon draw function

        for(vector<Feature>::iterator iter = features.begin(); iter < features.end(); iter++){
            vectorSize = (*iter).getPoints();
            vectorPoints = (*iter).getCoords();
            fillpoly(&vectorPoints[0], vectorSize -1);
        }
    }

    // Draws the names of all the streets given that they are not too close to each other
    void drawNames (vector<Street> & streets, t_bound_box mapBounds, int zoomLevel, double distance){
        LatLon position;

        string printName, extractName, appendName;
        double angle;
        for (vector<Street>::iterator iter = streets.begin(); iter < streets.end(); iter++) {
            position = findFirstInstance((*iter).getName(), mapBounds, angle, zoomLevel, takenSpaceStreet, distance);
            // Handle the case for a position being out of bounds or being too close to another name
            if (position.lat == 0 && position.lon == 0) {
                continue;
            }

            // Since we do not want to deal with excessively long strings (as they are both inconvenient and
            // interfere with drawing of other names) we limit the words for a street name to 3 words
            extractName = (*iter).getName();
            stringstream nameLimiter;
            nameLimiter<< extractName;
            for (int i = 0; i < 3; i++) {
                nameLimiter >> appendName;
                printName.append(appendName);
                printName.append(" ");
                appendName.clear();
            }

            //Lower the text size at the highest zoom levels to accomodate the display of more information
            if (zoomLevel < ZOOM7) {
                settextattrs(8, angle);
            } else {
                settextattrs(10,angle);
            }

            // Draw the name and clear strings used for parsing
            setcolor(BLACK);
            t_bound_box textBox(position.lon - 0.01, position.lat - 0.01, position.lon + 0.01, position.lat + 0.01);
            t_point centerPoints(position.lon, position.lat);
            drawtext_in(textBox, printName);
            printName.clear();
            appendName.clear();
            extractName.clear();
        }
    }

    // Finds the first instance (starting from the middle of the road) that is within the map
    // bounds and is not too close to any other written name
    LatLon findFirstInstance (string streetName, t_bound_box mapBounds, double & angle, int zoomLevel, vector<LatLon> & takenSpace, double distance) {
        LatLon fromCoord, toCoord, avgCoord;
        unsigned int delay;
        unsigned int i = 0;
        bool available = false;
        // Filter out all unknown roads (as this is not useful information) and roads that associate
        // with highways to prevent cluttering at highway intersections
        if(streetName.find("unknown") != -1 || streetName.find("Ramp") != -1 || streetName.find("Exit") != -1 ){
            avgCoord.lon = 0;
            avgCoord.lat = 0;
            return avgCoord;
        }

        // Set the boundaries of the map
        float boundaryRight = mapBounds.right();
        float boundaryLeft = mapBounds.left();
        float boundaryTop = mapBounds.top();
        float boundaryBottom = mapBounds.bottom();
        StreetSegmentEnds ends;

        // Get the number of segs and set the delay to be halfway through the street segments
        vector<unsigned> segIDs = find_street_street_segments(streetName);
        delay = (segIDs.size()/2);


        for(vector<unsigned>::iterator iter = segIDs.begin(); iter < segIDs.end(); iter++){

            if(i < delay){
                i++;
                continue;
            }

            available = false;

            ends = getStreetSegmentEnds(*iter);
            fromCoord = getIntersectionPosition(ends.from);
            toCoord = getIntersectionPosition(ends.to);
            unsigned int curveCount = getStreetSegmentCurvePointCount(*iter);
            unsigned int minCurveCount = *(segIDs.begin());

            // If there are no curve points on this street segment, then check if the names are within
            // bounds. If they are, check for the proximity of the
            if(curveCount == 0){
                if (( boundaryRight > fromCoord.lon) && (fromCoord.lon > boundaryLeft) && (fromCoord.lat < boundaryTop) && (fromCoord.lat > boundaryBottom) &&
                (boundaryRight > toCoord.lon) && (toCoord.lon > boundaryLeft) &&  (toCoord.lat < boundaryTop) && (toCoord.lat > boundaryBottom)  ){
                avgCoord.lon = (fromCoord.lon + toCoord.lon)/2;
                avgCoord.lat = (fromCoord.lat + toCoord.lat)/2;
                available = checkProximity(avgCoord, takenSpace, zoomLevel, distance);

                if(available == true){
                    angle = findAngle (fromCoord, toCoord);
                    takenSpaceStreet.push_back(avgCoord);
                    return avgCoord;
                } else {
                    avgCoord.lat = 0;
                    avgCoord.lon = 0;
                    return avgCoord;
                }


                }
            } else if (curveCount < minCurveCount){
                minCurveCount = curveCount;
                if (( boundaryRight > fromCoord.lon) && (fromCoord.lon > boundaryLeft) && (fromCoord.lat < boundaryTop) && (fromCoord.lat > boundaryBottom) &&
                (boundaryRight > toCoord.lon) && (toCoord.lon > boundaryLeft) &&  (toCoord.lat < boundaryTop) && (toCoord.lat > boundaryBottom)  ){
                avgCoord.lon = (fromCoord.lon + toCoord.lon)/2;
                avgCoord.lat = (fromCoord.lat + toCoord.lat)/2;
                available = checkProximity(avgCoord, takenSpace, zoomLevel, distance);

                if(available == true){
                    angle = findAngle (fromCoord, toCoord);
                    takenSpaceStreet.push_back(avgCoord);
                    return avgCoord;
                } else {
                    avgCoord.lat = 0;
                    avgCoord.lon = 0;
                    return avgCoord;
                }
                }
            }
        }
        avgCoord.lat = 0;
        avgCoord.lon = 0;
        return avgCoord;
    }

    void filterHighways (vector<Street> & highway){
        Highway current;
        vector<string> alreadyProcessed;
        vector<unsigned> segments;
        string nameConsidered, toPushBack;
        int pos1, pos2;

        bool numericalValue = false;
        bool nameExists = false;
        char theLetter;

        for(vector<Street>::iterator iterS = highway.begin(); iterS < highway.end(); iterS++){
            nameConsidered = (*iterS).getName();
            numericalValue = false;
            nameExists = false;

            for(string::iterator iterN = nameConsidered.begin(); iterN < nameConsidered.end(); iterN++){

                if (('0' < (*iterN)) && (*iterN) < '9'){
                    numericalValue = true;
                    theLetter = (*iterN);
                    break;
                }
        }
            if(numericalValue == true){
                pos1 = nameConsidered.find(theLetter);
                nameConsidered.erase(0, pos1);
                pos2 = nameConsidered.find_first_of(" ");
                toPushBack = nameConsidered.substr(0, pos2);


                // check if this number exists within the filteredHighways vector.
                // if so, simply add its street segments and street id to that Highway in the vector.
                for(vector<Highway>::iterator iterH = filteredHighways.begin(); iterH < filteredHighways.end(); iterH++){

                    if((*iterH).getName() == toPushBack){
                        nameExists = true;
                        (*iterH).setStreetID((*iterS).getID());
                        segments = find_street_street_segments((*iterS).getName());
                        (*iterH).setSegID(segments);
                        break;
                    }
                }

                if(nameExists ==  false){
                    alreadyProcessed.push_back(toPushBack);
                    current.setName(toPushBack);
                    current.setStreetID((*iterS).getID());
                    segments = find_street_street_segments((*iterS).getName());
                    current.setSegID(segments);
                    filteredHighways.push_back(current);
                }
            } else {

                current.setName(nameConsidered);
                current.setStreetID( (*iterS).getID() );
                segments = find_street_street_segments((*iterS).getName());
                current.setSegID(segments);
                filteredHighways.push_back(current);
            }
        }
    }

    void drawHighways (vector<Highway> & highways, t_bound_box mapBounds, int zoomLevel, double distance) {
        LatLon position;
        double angle;
        int count = 0;
        vector<unsigned> segs;
        for(vector<Highway>::iterator iter = filteredHighways.begin(); iter < filteredHighways.end(); iter++){
            segs = (*iter).getSegID();
            sort(segs.begin(), segs.end());
            position = findFirstInstanceOfHighway((*iter).getName(), segs, mapBounds, angle, count, zoomLevel, takenSpaceStreet, distance);
            count += 30;
            if (position.lat == 0 && position.lon == 0){
                continue;
            }
                settextattrs(10, angle);
                setcolor(BLACK);
                t_bound_box textBox(position.lon - 0.1, position.lat - 0.1, position.lon + 0.1, position.lat + 0.1);
                t_point centerPoints(position.lon, position.lat);
                drawtext(centerPoints, (*iter).getName(), textBox );


        }
    }

    // Finds the first instance that a highway name can be written, given a set delay
    LatLon findFirstInstanceOfHighway (string highwayName, vector<unsigned> & segs, t_bound_box mapBounds, double & angle, int count, int zoomLevel, vector<LatLon> & takenSpace, double distance) {
        LatLon fromCoord, toCoord, avgCoord;
        if(highwayName.find("Ramp") != -1 || highwayName.find("Exit") != -1 || highwayName.find("Lane") != -1 || highwayName.find("lane") != -1){
            avgCoord.lon = 0;
            avgCoord.lat = 0;
            return avgCoord;
        }

        bool available = true;
        float boundaryRight = mapBounds.right();
        float boundaryLeft = mapBounds.left();
        float boundaryTop = mapBounds.top();
        float boundaryBottom = mapBounds.bottom();
        StreetSegmentEnds ends;
        int firstConsider = 0;

        for(vector<unsigned>::iterator iter = segs.begin(); iter < segs.end(); iter++) {
            // Possesses a delay for the highways, ensuring that the names will not be situated
            // too close to each other
            if(firstConsider < count){
                firstConsider++;
                continue;
            }

            ends = getStreetSegmentEnds(*iter);
            fromCoord = getIntersectionPosition(ends.from);
            toCoord = getIntersectionPosition(ends.to);
            unsigned int curveCount = getStreetSegmentCurvePointCount(*iter);
            unsigned int minCurveCount = *(segs.begin());

            //Employ same methodology as findFirstInstance
            // See above code for more detailed comments
            if(curveCount == 0) {
                if (( boundaryRight > fromCoord.lon) && (fromCoord.lon > boundaryLeft) && (fromCoord.lat < boundaryTop) && (fromCoord.lat > boundaryBottom) &&
                (boundaryRight > toCoord.lon) && (toCoord.lon > boundaryLeft) &&  (toCoord.lat < boundaryTop) && (toCoord.lat > boundaryBottom)  ){
                avgCoord.lon = (fromCoord.lon + toCoord.lon)/2;
                avgCoord.lat = (fromCoord.lat + toCoord.lat)/2;
                available = checkProximity(avgCoord, takenSpace, zoomLevel, distance);

                if(available == true){
                    angle = findAngle (fromCoord, toCoord);
                    takenSpace.push_back(avgCoord);
                    return avgCoord;
                } else {
                    avgCoord.lat = 0;
                    avgCoord.lon = 0;
                    return avgCoord;
                }
            }

            } else if (curveCount < minCurveCount){
                minCurveCount = curveCount;
                if (( boundaryRight > fromCoord.lon) && (fromCoord.lon > boundaryLeft) && (fromCoord.lat < boundaryTop) && (fromCoord.lat > boundaryBottom) &&
                (boundaryRight > toCoord.lon) && (toCoord.lon > boundaryLeft) &&  (toCoord.lat < boundaryTop) && (toCoord.lat > boundaryBottom)  ){
                avgCoord.lon = (fromCoord.lon + toCoord.lon)/2;
                avgCoord.lat = (fromCoord.lat + toCoord.lat)/2;
                available = checkProximity(avgCoord, takenSpace, zoomLevel, distance);

                if(available == true){
                    angle = findAngle (fromCoord, toCoord);
                    takenSpace.push_back(avgCoord);
                    return avgCoord;
                } else {
                    avgCoord.lat = 0;
                    avgCoord.lon = 0;
                    return avgCoord;
                }
            }



        }
    }
        avgCoord.lat = 0;
        avgCoord.lon = 0;
        return avgCoord;
    }

    // Checks that there is enough space to write a name based on the zoom level
    // Called for both POI and for streets
    bool checkProximity (LatLon position, vector<LatLon> & takenSpace, int zoomLevel, double distance){
        double radius;
        bool writePossible = true;
        switch(zoomLevel){
            case ZOOM0:
                    radius = 16000;
                    for(vector<LatLon>::iterator iter = takenSpace.begin(); iter < takenSpace.end(); iter++){
                        if( find_distance_between_two_points((*iter), position) < radius){
                            writePossible = false;
                            break;
                        }
                    }
                    break;

                case ZOOM1:
                    radius = 12000;
                    for(vector<LatLon>::iterator iter = takenSpace.begin(); iter < takenSpace.end(); iter++){
                        if( find_distance_between_two_points((*iter), position) < radius ){
                            writePossible = false;
                            break;
                        }
                    }
                    break;

                case ZOOM2:
                    radius = 8000;
                    for(vector<LatLon>::iterator iter = takenSpace.begin(); iter < takenSpace.end(); iter++){
                        if( find_distance_between_two_points((*iter), position) < radius){
                            writePossible = false;
                            break;
                        }
                    }
                    break;

                case ZOOM3:
                    radius = 4000;
                    for(vector<LatLon>::iterator iter = takenSpace.begin(); iter < takenSpace.end(); iter++){
                        if( find_distance_between_two_points((*iter), position) < radius ){
                            writePossible = false;
                            break;
                        }
                    }
                    break;

                case (ZOOM4):
                    radius = 1500;
                    for(vector<LatLon>::iterator iter = takenSpace.begin(); iter < takenSpace.end(); iter++){
                        if(  find_distance_between_two_points((*iter), position) < radius){
                            writePossible = false;
                            break;
                        }
                    }
                    break;

                case (ZOOM5):
                    radius = 1000;
                    for(vector<LatLon>::iterator iter = takenSpace.begin(); iter < takenSpace.end(); iter++){
                        if( find_distance_between_two_points((*iter), position) < radius ){
                            writePossible = false;
                            break;
                        }
                    }
                    break;

                case (ZOOM6):
                    radius = 500;
                    for(vector<LatLon>::iterator iter = takenSpace.begin(); iter < takenSpace.end(); iter++){
                        if( find_distance_between_two_points((*iter), position) < radius ){
                            writePossible = false;
                            break;
                        }
                    }
                    break;

                case (ZOOM7):
                    radius = 400;
                    for(vector<LatLon>::iterator iter = takenSpace.begin(); iter < takenSpace.end(); iter++){
                        if( find_distance_between_two_points((*iter), position) < radius){
                            writePossible = false;
                            break;
                        }
                    }
                    break;

                case (ZOOM8):
                    radius = 200;
                    for(vector<LatLon>::iterator iter = takenSpace.begin(); iter < takenSpace.end(); iter++){
                        if(find_distance_between_two_points((*iter), position) < radius){
                            writePossible = false;
                            break;
                        }
                    }
                    break;

                case (ZOOM9):
                    radius = 100;
                    //radius = distance / 100;
                    for(vector<LatLon>::iterator iter = takenSpace.begin(); iter < takenSpace.end(); iter++){
                        if( find_distance_between_two_points((*iter), position) < radius){
                            writePossible = false;
                            break;
                        }
                    }
                    break;

            default:
                return writePossible;


        }
        return writePossible;
    }

/*------------------------Points of Interest-----------------------------*/
    void POIon(void (*drawscreen) (void)){
        // toggles the flag for draw
        drawFlag = !drawFlag;
        if (drawFlag == false) {
            change_button_text("Turn POI Off", "Turn POI On");
        }
        else {
            change_button_text("Turn POI On", "Turn POI Off");
        }
        drawscreen();
    }

    void sortPOI(){
        unsigned long long totalPOI = getNumberOfPointsOfInterest();
        string namePOI;
        t_point coordsPOI;
        for (unsigned i = 0; i < totalPOI; i++){
            namePOI = getPointOfInterestName(i);
           coordsPOI.x = getPointOfInterestPosition(i).lon;
            coordsPOI.y = getPointOfInterestPosition(i).lat;

            POI pointsOfInterest(namePOI, coordsPOI);
            poiVector.push_back(pointsOfInterest);
        }
    }

    void drawPOI(double distance, int zoomLevel, t_bound_box mapBounds) {
        // only draws POI if the drawFlag is raised (true)
        // doesn't do anything otherwise (hopefully)
        if (drawFlag){
            bool available = true;
            LatLon poiPos;
            double radius = distance / 50000000; //This will be used to draw the points of interest.
            // loops through the vector and draws all of the points of vector

            float boundaryRight = mapBounds.right();
            float boundaryLeft = mapBounds.left();
            float boundaryTop = mapBounds.top();
            float boundaryBottom = mapBounds.bottom();

            for (unsigned i = 0; i < poiVector.size(); i++) {
                available = true;
                t_point poiCoords = poiVector[i].getCoordsPOI();
                poiPos.lon = poiCoords.x;
                poiPos.lat = poiCoords.y;

                if( ( boundaryRight > poiPos.lon) && (poiPos.lon > boundaryLeft) && (poiPos.lat < boundaryTop) && (poiPos.lat > boundaryBottom)){

                    setcolor(RED); // all the points of interest show up as red
                    fillarc(poiCoords, radius, 0, 360);

                    if(zoomLevel <= ZOOM7){
                        available = checkProximity(poiPos, takenSpaceStreet, zoomLevel, distance);

                        if(available == true){
                            takenSpaceStreet.push_back(poiPos);
                            settextattrs(8, 0);
                            setcolor(BLACK);
                            t_bound_box textBox(poiCoords.x - 0.1, poiCoords.y - 0.1, poiCoords.x + 0.1, poiCoords.y + 0.1);
                            t_point centerPoints(poiCoords.x, poiCoords.y);
                            drawtext(centerPoints, poiVector[i].getNamePOI(), textBox );
                        }
                    }
                }
            }
        }
    }

/*-------------------------Draw Route Things!----------------------------*/
    // allows the user to select things
    void locationSelect (float x, float y, t_event_buttonPressed button_info) {
        // creates a LatLon variable and stores the x,y values at LatLon.
        // stores this into glbal locationArray
        LatLon location;
        location.lat = y;
        location.lon = x;

        // gets the closest intersection to the point that was clicked, gets its position, stores it
        unsigned intersectionID = findClosestIntersection(location);

        // get the thing
        if (locationNum < 2) {
            cout << endl << "You have selected: " << getIntersectionName(intersectionID) << endl;
        }

        // check if the thing is correct!
        if ((locationNum == 1) && (intersectionID == intersectionIDArray[0])){
            cout << endl << "Error: " << getIntersectionName(intersectionID) << " has already been selected. Please try again." << endl;
        }
        else {
            intersectionIDArray[locationNum] = intersectionID; // stores the intersection ID here
            LatLon intersectionPos = getIntersectionPosition(intersectionID);
            intersectionArray[locationNum] = intersectionPos; // stores the position LatLon of the intersection

            locationNum++; // increment things (highest is 2, where it resets itself)
            draw_screen(); // draws the screen with the route
        }
    }

    void drawRouteEnds() {
        vector<float> bounds = findCoordinates();
            t_bound_box mapBounds = get_visible_world();
            double distance; //A measure of the diagonal of the screen, in metres.

            vector<float> testBounds = findCoordinates();

            LatLon bottomLeft, topRight;
            bottomLeft.lat = mapBounds.bottom();
            bottomLeft.lon = mapBounds.left();
            topRight.lat = mapBounds.top();
            topRight.lon = mapBounds.right();


            distance = find_distance_between_two_points(bottomLeft, topRight);

        double radius = distance / 13500000;
            t_point centrePoint; // the fillarc function accepts sacrifices in t_point so need to convert
            centrePoint.x = intersectionArray[0].lon;
            centrePoint.y = intersectionArray[0].lat;
        if (locationNum == 1) {
            setcolor(START); // bright green colour!
            fillarc(centrePoint, radius, 0, 360); // draws a circle
        }
        else if (locationNum == 2) {
            setcolor(STOP); // bright red colour!
            fillrect(intersectionArray[1].lon-radius, intersectionArray[1].lat-radius, intersectionArray[1].lon+radius, intersectionArray[1].lat+radius); // if the user is sloppy while clicking then they get ugly things
            setcolor(START); // bright green colour!
            fillarc(centrePoint, radius, 0, 360); // draws a circle
        }
    }

    void drawRoute(unsigned int width, t_bound_box mapBounds){
        if (drawPathFlag){ // only draws if the flag is set
            if (locationNum < 3) { // only draws if there are fewer than 3 points selected
                if (locationNum == 2) { // only draws if there are two locations selected, 0 and 1
                    path = find_path_between_intersections(intersectionIDArray[0], intersectionIDArray[1]); // the [0] has start and [1] has end
                    if (path.size() == 0) { // path is empty, draw nothing and clear the screen of location markers
                        cout << endl << "Error: no valid path found. Please try a different location. Locations will be cleared." << endl;
                        locationNum = 0; // resets things
                        draw_screen(); // redraws the screen ayyy
                        printAlready = false; // resets the flag
                    }
                    else { // only draw a route if the vector isn't empty...
                        if (!printAlready){ // find the path to draw if not already done
                            // Ivan's function to find shortest path and things
                            drawPath(path, width, mapBounds); // this calls Ivan's draw functions
                            printDirections(path, intersectionIDArray[0], intersectionIDArray[1]); // this prints out directions to terminal
                            printAlready = true; // changes the flag to signified already printed
                        }
                        else { // just draw the path (without recalculating the path)
                            // this calls Ivan's draw functions
                            drawPath(path, width, mapBounds);
                        }
                    }
                }
                drawRouteEnds(); // draws these on top
            }
            else { // on third selection, clears things
                locationNum = 0; // resets things
                draw_screen(); // redraws the screen ayyy
                printAlready = false; // resets the flag
            }
        }
    }

    void printHelp(void (*drawscreen) (void)){
        cout << endl<<
                "Click anywhere on the map to set a start and end destination." << endl <<
                "The program will automatically set it to the closest intersection." << endl <<
                "\n Click once to set your start destination (a green marker will show)" << endl <<
                "and click again to set your end destination (a red marker will show) " << endl <<
                "The shortest path will be drawn between the two points, and directions will be printed here." << endl <<
                "Click again anywhere on the screen to remove all previous selections."
                << endl;
        drawscreen();
    }

/*-------------------------Draws TSP Things!-----------------------------*/
    // can only have one type of selection at a time, either single path or TSP (depot or dropoff)

    void toggleClear(void (*drawscreen) (void)){
        // deletes all the depots and the path
        depot.clear();
        dropoff.clear();
        cout << "cleared others" << endl;
        TSPpath.clear();
        cout << "cleared TSP" << endl;
        drawTSPalready = false;
        draw_screen();
    }

    void pathOn(void (*drawscreen) (void)){
        drawPathFlag = !drawPathFlag;
        // toggles other selection flags to off
        drawDropoffFlag = false;
        drawDepotFlag = false;
        if (!drawPathFlag){
            locationNum = 0; // clears to 0
        }
        draw_screen();
    }

    void depotOn(void (*drawscreen) (void)){ // toggles the depot location selection flag
        drawDepotFlag = !drawDepotFlag;
        // toggles other selection flags to off
        drawPathFlag = false;
        drawDropoffFlag = false;
        locationNum = 0; // sets the single path number to 0
        //cout << endl << "DEPOT FLAG ONNNN: " << drawDepotFlag << endl;
        draw_screen();
    }

    void delivOn(void (*drawscreen) (void)){ // toggles the dropoff location selection flag
        drawDropoffFlag = !drawDropoffFlag;
        // toggles other selection flags to off
        drawPathFlag = false;
        drawDepotFlag = false;
        locationNum = 0; // sets the single path number to 0
        //cout << endl << "DELIV FLAG ONNNN: " << drawDropoffFlag <<  endl;
        draw_screen();
    }

    void pointSelect (float x, float y, t_event_buttonPressed button_info) {
        // creates a LatLon variable and stores the x,y values at LatLon.
        // stores this into global locationArray
        selectMessage = true;
        alreadyAdded = false;
        LatLon location;
        location.lat = y;
        location.lon = x;

        // gets the closest intersection to the point that was clicked, gets its position, stores it
        unsigned intersectionID = findClosestIntersection(location);
        displayName = getIntersectionName(intersectionID);

        cout << endl << "You have selected: " << displayName << endl;

        if (drawPathFlag){
            // check if the thing is correct!
            if ((locationNum == 1) && (intersectionID == intersectionIDArray[0])){
                cout << endl << "Error: " << displayName << " has already been selected. Please try again." << endl;
            }
            else {
                intersectionIDArray[locationNum] = intersectionID; // stores the intersection ID here
                LatLon intersectionPos = getIntersectionPosition(intersectionID);
                intersectionArray[locationNum] = intersectionPos; // stores the position LatLon of the intersection

                locationNum++; // increment things (highest is 2, where it resets itself)
            }
        } // end drawpathflag check

        else if (drawDepotFlag){ // if that's set to on
            if (depot.size() == 0) { // nothing there
                depot.push_back(intersectionID);
                cout << endl << displayName << " has been added to the list of depots." << endl;
                //cout << endl << "depot: " << depot.size() << endl;
            }
            else {
                for (unsigned i = 0; i < depot.size(); i++){ // iterate to find duplicates
                    if (intersectionID == depot[i]) { // already selected
                        cout << endl << "Error: " << displayName << " has already been selected. Please try again." << endl;
                        alreadyAdded = true;
                        break;
                    }
                } // end for loop (duplicate check)
                if (!alreadyAdded){ // not yet selected, add to vector
                    depot.push_back(intersectionID);
                    cout << endl << displayName << " has been added to the list of depots." << endl;
                    //cout << endl << "depot: " << depot.size() << endl;
                }
            }
        }

        else if (drawDropoffFlag){ // if that's set to on
            if (dropoff.size() == 0) { // empty vector, just push_back
                dropoff.push_back(intersectionID);
                cout << endl << displayName << " has been added to the list of dropoff points" << endl;
                //cout << endl << "dropoff: " << dropoff.size() << endl;
            }
            else { // things exist, must check if already selected
                for (unsigned i = 0; i < dropoff.size(); i++){ // iterate to find duplicates
                    if (intersectionID == dropoff[i]) { // if exists, exit exit exit!
                        cout << endl << "Error: " << displayName << " has already been selected. Please try again." << endl;
                        alreadyAdded = true;
                        break;
                    }
                } // end for loop (duplicate check)
                if (!alreadyAdded){
                   dropoff.push_back(intersectionID);
                   cout << endl << displayName << " has been added to the list of dropoff points" << endl;
                   //cout << endl << "dropoff: " << dropoff.size() << endl;
                }
            } // end dropoff size check
        }

        draw_screen(); // draws the screen with the route
    }

    void drawPoints() {
         vector<float> bounds = findCoordinates();
            t_bound_box mapBounds = get_visible_world();
            double distance; //A measure of the diagonal of the screen, in metres.

            vector<float> testBounds = findCoordinates();

            LatLon bottomLeft, topRight;
            bottomLeft.lat = mapBounds.bottom();
            bottomLeft.lon = mapBounds.left();
            topRight.lat = mapBounds.top();
            topRight.lon = mapBounds.right();

            distance = find_distance_between_two_points(bottomLeft, topRight);

            double radius = distance / 13500000;
            t_point centrePoint; // the fillarc function accepts sacrifices in t_point so need to convert

        setcolor(MAGENTA);
        for (unsigned i = 0; i < depot.size(); i++){
            //cout << "draw something i'm giving up on u" << endl;
            LatLon intersectionPos = getIntersectionPosition(depot[i]);
            centrePoint.x = intersectionPos.lon;
            centrePoint.y = intersectionPos.lat;
            fillarc(centrePoint, radius, 0, 360); // draws a circle
        }

        setcolor(YELLOW);
        for (unsigned i = 0; i < dropoff.size(); i++){
            //cout << "draw something i'm giving up on u" << endl;
            LatLon intersectionPos = getIntersectionPosition(dropoff[i]);
            centrePoint.x = intersectionPos.lon;
            centrePoint.y = intersectionPos.lat;
            fillarc(centrePoint, radius, 0, 360); // draws a circle
        }
    }

    void drawTSP(unsigned int width, t_bound_box mapBounds){
        if (drawTSPflag){
            if (!drawTSPalready){
                TSPpath = traveling_salesman(dropoff, depot); // these are global variables containing the path
                if (TSPpath.size() != 0){ // only if TSPpath returns a valid path
                    drawTSPalready = true; // toggles it to true
                }
            }
            drawPath(TSPpath, width, mapBounds); // this calls Ivan's draw functions
            drawPoints(); // draw these on top of the selected things
        }
    }

    // draws TSP when told to
    void runTSP(void (*drawscreen) (void)){
        drawTSPflag = true;
        draw_screen();
    }
/*--------------------------Draws a sidebar!-----------------------------*/
    void toggleSearch(void (*drawscreen) (void)){
        // toggles the flag for drawing sidebars
        drawSearchFlag = !drawSearchFlag;
        if (drawSearchFlag == false) {
            change_button_text("Hide Search", "Show Search");
        }
        else {
            change_button_text("Show Search", "Hide Search");
        }
        drawscreen();
    }

    void drawSearch() {
        if (drawSearchFlag){
            vector<float> bounds = findCoordinates();
            t_bound_box mapBounds = get_visible_world();
            double distance; //A measure of the diagonal of the screen, in metres.

            vector<float> testBounds = findCoordinates();

            LatLon bottomLeft, topRight;
            bottomLeft.lat = mapBounds.bottom();
            bottomLeft.lon = mapBounds.left();
            topRight.lat = mapBounds.top();
            topRight.lon = mapBounds.right();

            distance = mapBounds.right() - mapBounds.left();
            double width = distance/5.0;
            double offset = distance/30.0;
            double buttonwidth = width/5.0;
            double buttonheight = offset/2.0;

            // top left point
            t_point topLeftPoint;
            topLeftPoint.x = mapBounds.left()+offset;
            topLeftPoint.y = mapBounds.top()-offset;

            // top right point
            t_point topRightPoint = topLeftPoint;
            topRightPoint.x += width;

            // bottom left point
            t_point bottomLeftPoint = topLeftPoint;
            bottomLeftPoint.y -= buttonheight*5.5;

            // bottom right point
            t_point bottomRightPoint = bottomLeftPoint;
            bottomRightPoint.x += width;

            setcolor(WHITE);
            fillrect(bottomLeftPoint, topRightPoint);
            setcolor(CYAN);
            setlinewidth(4);
            drawline(topRightPoint, topLeftPoint);
            drawline(bottomRightPoint, bottomLeftPoint);
            drawline(topRightPoint, bottomRightPoint);
            drawline(topLeftPoint, bottomLeftPoint);

            // draw two buttons

            t_point buttonAddBLP, buttonAddTRP;
            t_point buttonSearchBLP, buttonSearchTRP;

            //  _________________________________________
            // |                                        |
            // |       _________       __________       |
            // |       |  Add  |       | Search |       |
            // |                                        |
            // the first button is add, the second is search
            buttonAddBLP.x = topLeftPoint.x + buttonwidth;
            buttonAddBLP.y = topLeftPoint.y - 5*buttonheight;

            buttonAddTRP.x = topLeftPoint.x + 2*buttonwidth;
            buttonAddTRP.y = topLeftPoint.y - 4*buttonheight;

            // add second button

            buttonSearchTRP = buttonAddTRP;
            buttonSearchBLP = buttonAddBLP;

            buttonSearchBLP.x += 2*buttonwidth;
            buttonSearchTRP.x += 2*buttonwidth;

            // draws the two buttons
            setcolor(GREEN);
            fillrect(buttonAddBLP, buttonAddTRP);

            setcolor(MAGENTA);
            fillrect(buttonSearchBLP, buttonSearchTRP);
        }
    }

    // void searchButton(float x, float y, t_event_buttonPressed button_info) {
    //     if
    // }
