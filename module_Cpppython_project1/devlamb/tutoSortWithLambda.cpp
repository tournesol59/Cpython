#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;
// routine de base de lalgorithm quicksort, un element pris au hasard
// sans expression lambda pour commencer

int qsort_core(std::vector<int> vec, int middle) 
{
   // une structure de donnees vec qui est modifiable  
  std::vector<int>::iterator it1 = vec.begin();
  std::vector<int>::iterator it2 = vec.begin();
  std::vector<int>::iterator it3 = vec.begin();
  // place the iterator it2 at the middle
  for (int i=0; i < middle; i++) {
    ++it2;
  }
  // swap vec(middle) and it1
   int temp=*it1;
   *it1 = *it2;
   *it2 = temp;
      // now the middle-tobe element is placed at the beginning of the vector table,
      // it2 must be init again
   it2 = vec.begin();
   it1 = std::next(vec.begin());
   it3 = std::next(it2);
  // while loop
  int Nv = vec.size();
  int i=0;
  while (i < Nv-1) {
       if (*it2 < *it1) {
         //nothing but advance it2
         ++it1; i=i+1;
       }
       else if (*it2 > *it1) {
         // it2 must be pushed one time but it2++ and it1 must be swaped before
            // if it2 and it1=next(it2) then a simple swap is sufficient
            if (it1 == next(it2)) {
                 temp = *it1;
                 *it1= *it2;  
                 *it2 = temp;
                 it2 = it1;
                 ++it1; i=i+1;
            }
            //else if it1> next(it2), a triple swap must be done, except if it1 is already the lastelemt
            else if (i<Nv-2) {
                 temp = *it1;
                 *it1 = *it3;
                 *it3 = temp;
                 
                 temp=*it3;
                 *it3 = *it2;
                 *it2 = temp;
                 ++it2;
                 it3=std::next(it2);
                 ++it1; i=i+1;
            } 
            // else: at the end of the table a simple swap is sufficient
            else {
              temp = *it1;
              *it1 = *it2;
              *it2 = temp;
              it2 = it1; i=i+1;
              // ++it1; // already at the end
            }
       }
  }
  return 0;
}


int main() {

  std::cout <<"demo prix a vector of string with the same prefix using lambdas\n";
  std:vector<string> v({"9", "12", "13", "14", "18", "42", "75"});

  std::string prefix("hi");
  
  std::for_each(v.begin(), v.end(), [&prefix](const auto &a) {
       std::cout << prefix + a << std::endl;
  });     
  
  std::cout << "demo of a vector of MyClass{int, string} be sorted using lambda functions \n";

  std::vector<std::pair<int,string>> myvec;
  myvec.push_back(std::pair<int, string>({1,"abc"}));
  myvec.push_back(std::pair<int, string>({3,"aef"}));
  myvec.push_back(std::pair<int, string>({6,"aaz"}));
  myvec.push_back(std::pair<int, string>({5,"efg"}));
  myvec.push_back(std::pair<int, string>({2,"dac"}));
  myvec.push_back(std::pair<int, string>({4,"czu"}));

    // now define and call my qsort_core function
  std::vector<int> vec({-1,-2,-6,-4,5,1});
  
  std::cout << "Before sort core" << std::endl;
  for (auto it=vec.begin(); it != vec.end(); ++it) {
       std::cout << *it << std::endl;
  } 
  qsort_core(vec, 2);
  
  std::cout << "After sort core" << std::endl;
  for (auto it=vec.begin(); it != vec.end(); ++it) {
       std::cout << *it << std::endl;
  }    
     
  // now define and call the lambda
  //std::vector<int> vec({-1,-2,-6,-4,5,1});

  std::sort(vec.begin(), vec.end(),
              [](const auto& a, const auto& b) {
                return a < b; });
                
  for (auto it=vec.begin(); it != vec.end(); ++it) {
       std::cout << *it << std::endl;
  }
 /*
  std::cout << "demo of classic error"
  std::vector<int> kk;
  kk.push_back(2);
  std::function<int(int)> foo = std::function<int(int)> foo{[kk](intx) mutable
  {
    kk.push_back(1);
    return kk[0]+1;
  }};
  */
  
  
  return 0;
}
