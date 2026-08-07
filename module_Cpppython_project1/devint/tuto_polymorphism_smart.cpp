#include <iostream>
#include <memory>
#include <assert.h>

using namespace std;
/*** la classe abstraite de base ***/
class Object
{
   unsigned long int new_handle(void);
   
   protected:
      unsigned long int h;
   
   public:
      Object(void);  // le constructeur
      virtual ~Object(void);  // le destructeur virtuel, oblige pour polymorphism
      virtual void print(void) =0; // fonction virtuelle pure
      unsigned long int handle(void) const;  // Fonction renvoyant le numero d'identifiant de l'objet
    
};

// cette fonction n est appelable que par la classe Objet
unsigned long int Object::new_handle(void)
{
   static unsigned long int hc=0;
   return hc=hc+1;    
}

// le constructeur de Object doit etre appele par les classe derivees:
Object::Object(void) {
   h = new_handle();
   return;   
}

Object::~Object(void)
{
   return;   
}

unsigned long int Object::handle(void) const
{
   return h;  // renvoie le numero de l'objet   
}

/*** la classe sac ***/
class Bag : public Object // la classe sac. Elle herite de Object
{
   struct BagList {
      unique_ptr<Object> ptr;
      unique_ptr<BagList> next;  
   };
   
   unique_ptr<BagList> head;
      
   public:
      Bag(void);
      ~Bag(void);
      void print();
      bool has(int h);
      void add(std::unique_ptr<Object> ob);
      //void remove(Object &o);
};

Bag::Bag() : Object()
{
   return;
}

Bag::~Bag() 
{
   while (head) {
      auto next= std::move(head->next);
      head = std::move(next);  // smart pointer shall dispose automatically
   }
}
 
void Bag::print() 
{
   cout << "Sac numero " << handle() << "." << endl;
   std::unique_ptr<BagList> temp = std::move(head->next); 
   while (temp) { // il y a un petit mais: si head est null
      //temp->print(); // a tester marche pas, le pointeur est sur BagList, pas sur bag
      cout << "Sac numero " << temp->ptr->handle() << "." << endl;
      temp = std::move(temp->next);
   }
   return;
}

bool Bag::has(int h)
{
   if (handle()==h)
      return true;
      std::unique_ptr<BagList> temp = std::move(head->next);
      while (temp) {
         if (temp->ptr->handle()==h) 
            return true;
         temp = std::move(temp->next); 
      }
      return false;
}

void Bag::add(std::unique_ptr<Object> ob) 
{
   head = std::unique_ptr<BagList> (new BagList{ std::move(ob), std::move(head)});  
   return;
}


/*** pour utiliser une instance Object, il faut definir une classe non abstraite  ***/
class MonObjet : public Object
{
    /*  Définir la méthode print() pour l'objet...  */
    public: 
       void print(void) override;
};

void MonObjet::print(void)
{
   cout << "Mon objet: " << handle()  << endl;    
}


// a function consuming a unique_ptr can take it by value of by rvalue reference
std::unique_ptr<MonObjet> pass_through(std::unique_ptr<MonObjet> p) {
   p->print();
   return p;
}



Bag monSac;

int main() {
       std::cout << "1) Unique ownership semantics demo\n";
   {
   std::unique_ptr<MonObjet> a =  std::make_unique<MonObjet>();
   std::unique_ptr<MonObjet> qa = pass_through(std::move(a));
   
   assert(!a);
   //b,c;

   }
      
   std::cout << "\n" "2) Runtime polymorphism demo\n";
    {
        // Create a derived resource and point to it via base type
        std::unique_ptr<Object> pa =  std::make_unique<MonObjet>();
        // Dynamic dispatch works as expected
        pa->print();

    }
     
    std::cout << "\n" "3) Runtime polymorphism on List\n";         
    {
        // Create a derived resource and point to it via base type
        std::unique_ptr<Object> pa =  std::make_unique<MonObjet>();
        monSac.add(pa);
   //monSac.add(b);
   //monSac.add(c);
   
   //   monSac.print();
    }

   return 0;   
}

