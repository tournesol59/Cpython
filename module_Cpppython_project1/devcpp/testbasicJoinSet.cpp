#include <iostream>
#include <string>
#include <cstring>
#include <vector>

#include "idjoinset.h"

using namespace std;

int main() {
   const char* email = "sausan@aol.com";
  // Insurance insurance = WENDINS;
   const int insurance = 1;
   const char* refnum = "BF6A4N";
   const char* departure = "Toulouse";
   const char* destination = "Bayonne";
   const char* spare = "210";
   const int spareref = 1;
   
   JOIN_SET* js_p = create_joinset(email, insurance, refnum, departure, destination, spare, spareref);
   vector<JOIN_SET> joinset_vec;
   joinset_vec.push_back(*js_p);
   
   for (vector<JOIN_SET>::iterator it=joinset_vec.begin(); it != joinset_vec.end(); ++it) {
       cout << "The first element of our join-result record list\n";
       //printf("%s \t %d \t %s \t %s \t %s \t %s \t %s \n", *it->emailaddr, *it->insurance, *it->refnum, *it->departure, *it->destination, *it->spare, *it->spareref);   
       printf("%s \n", it->emailaddr); 
   }

   //joinset_vec.remove();
   return 0;
}
