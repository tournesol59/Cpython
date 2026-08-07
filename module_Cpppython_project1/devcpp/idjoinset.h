#include <iostream>
#include <cstring>
using namespace std;
//Definition enum
enum Insurance {
   WENDINS,  // carte week end
   FREQINS, // frequence plus
   BIGINS   // grand voyageur
};

typedef struct JoinSet JOIN_SET;

struct JoinSet {
   // !!! int id; //!!!
   char emailaddr[40];
   int insurance;  // enum
   char refnum[40];
   char departure[40];
   char destination[40];
   char spare[40];  // spare field for future, eg price
   int spareref;  // spare refnum for future, eg category dependent on price
};
typedef JOIN_SET* JoinSet_Ptr;

JoinSet_Ptr create_joinset(const char* email, int insurance, const char* refnum, const char* departure, const char* destination, const char* spare, int spareref);


