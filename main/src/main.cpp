#include "StreetsDatabaseAPI.h"
#include "m1.h"
#include "m2.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include "m3.h"
#include "m4.h"

using namespace std;


int main(int argc, char* argv[]) {
	string city = argv[1];
	string bin = "/cad2/ece297s/public/maps/" + city + ".bin";
	load_map(bin);

    close_map();



    return 0;
}
