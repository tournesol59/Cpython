#include "idjoinset.h"

using namespace std;

JoinSet_Ptr create_joinset(const char* email, int insurance, const char* refnum, const char* departure, const char* destination, const char* spare, int spareref)
{
   JoinSet_Ptr js_p = (JOIN_SET*) malloc(1*sizeof(JOIN_SET));   
   if (js_p != NULL) 
   {
      strcpy(js_p->emailaddr, email);
      js_p->insurance = insurance;
      strcpy(js_p->refnum, refnum);
      strcpy(js_p->departure, departure);
      strcpy(js_p->destination, destination);
      strcpy(js_p->spare, spare);
      js_p->spareref = spareref;
      
      return js_p;
   }
   else return NULL;   
}
