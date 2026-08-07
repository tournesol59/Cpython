#include <iostream>
#include <cstring>
#include <list>
#include <locale>
#include <memory>

using namespace std;

// equivaut a la classe Component du Design Pattern Composite
class Menu {
   public:
      // le destructeur d une interface doit etre declare virtuel
      virtual ~Menu() {}
   
   public:
      // Parcours recursif de l arbre   
      virtual void ls_r() = 0; 
      // virtuelle pure doit etre redefini ainsi par son implementation - Selection / Operation
       
      virtual void scan_d() = 0; // equivaut a "ajouter selection" car va prendre un entier comme propriete

      //virtual void scan_s() = 0;  // APRES
   protected:
      int select_;
      string name_;
      
};

class Operation : public Menu 
  {
  virtual void ls_r() override;
  
  virtual void scan_d() override;
     
  protected:
      int select_;
      string name_;
  };
  
class Selection : public Menu 
  {
  Selection(int dselect, string dname);
  virtual void ls_r() override;
  
  virtual void scan_d() override;
  
  protected:
      int select_;
      string name_;
     //std::list<std::unique_ptr<Menu>> menu_ptrs_; // on teste d'abord plus petit: ci-dessous
     std::unique_ptr<Menu> menu_ptr;
  };



