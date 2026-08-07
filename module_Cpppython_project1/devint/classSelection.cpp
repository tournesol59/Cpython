#include "classSelection.h" // d'abord ne faire que la decl de la classe de base

#include <iostream>
#include <cstring>

using namespace std;
/*
// equivaut a la classe Component du Design Pattern Composite
class Menu {
   public:
      // le destructeur d une interface doit etre declare virtuel
      virtual ~Component() {}
   
   public:
      // Parcours recursif de l arbre   
      virtual void ls_r() = 0; 
      // virtuelle pure doit etre redefini ainsi par son implementation - Selection / Operation
       
      virtual void scanf_d() = 0; // equivaut a "ajouter selection" car va prendre un entier comme propriete

      //virtual void scanf_s() = 0;  // APRES
   protected:
      id select_;
      string name_;
      
};
*/

// equivaut a la classe Leaf du design pattern
//class Operation : public Menu {
//   public:
//};


void Operation :: ls_r()
     {
         cout << "\t\t selected: " << select_ << "name: " << name_ << endl;
     }     

void Operation :: scan_d() 
      {
         cout << "Entrez votre selection: 1=faire une recherche sur email, 2=faire une recherche sur refnum" << endl;  
         cin >> select_;
         cout << "Vous avez entre: " << select_ << endl;
         switch (select_) {
            case 1:
              cout << "pointe vers recherche par Email" << endl;
              break;
            case 2:
              cout << "pointe vers recherche par Refnum" << endl;
              break;
         }
            
      }
      
// equivaut a la classe Composite du design pattern
//class Selection : public Menu {
//   public:
//};
Selection::Selection(int dselect, string dname) : select_(dselect), name_(dname)
{
}

void Selection :: ls_r()
      {
         cout << "selected: " << select_ << "name: " << name_ << endl;
         //for (auto &menu_ptr : menu_ptrs_) 
         //{
         // fred: menu ptr est un unique ptr genere de base, 
         //   iterateur sur la liste de smart pointer suivante:
           //menu_ptr->ls_r();     
         //}
         menu_ptr->ls_r();  

      }
       
void Selection :: scan_d() 
      {
         cout << "Entrez votre selection: 1=faire une recherche, 2=faire une modification" << endl;  
         cin >> select_;
         cout << "Vous avez entre: " << select_ << endl;
         
         
         std::unique_ptr<Menu> menu_ptr = std::make_unique<Menu>(new Selection({1,"name"})); // cause a compilation error !!!
         //menu_ptrs_.push_back(q);
         //menu_ptr=q;
         
      }