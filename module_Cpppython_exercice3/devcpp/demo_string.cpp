#include <iostream>
#include <cstring>
#include <string>
#include <vector>
#include <sstream>

using namespace std;
vector<int> parseInts(string str) {
    // Complete this function
    istringstream iss(str);
    vector<int> res;
    int x;
    string temp;
    char delim = '-'; // whatever you want
    while (getline(iss, temp, delim)) {
        if (istringstream(temp) >> x) { // or std::stoi(), std::strtol(), etc
            res.push_back(x);
            cout << x;
        }
    }
    return res;
}

int main()   {
   string mystr("1-2-3");
   vector<int> myvect=parseInts(mystr);
   cout << endl;
   return 0;
}