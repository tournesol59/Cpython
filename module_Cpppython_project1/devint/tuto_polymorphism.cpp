#include <iostream>

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
   //protected:
   //   int value;
   //Object *ptr_next; // non ce n'est pas une liste chaine directement: elle contient une declaration
      // de liste qui contient l'objet de base.
      // En resume Bag peut contenir une liste qui contient Bag lui meme mais via la classe abstract
      //dont elle herite
      
   struct BagList
   {
      BagList *next;
      Object *ptr;
   };
   
   BagList *head;
   
   public:
      Bag(void);
      ~Bag(void);
      void print(void);  // Fonction d affichage du sac
      bool has(unsigned long int) const;
      void add(Object &ob);
    //  void remove(Object &ob);
       
};
 /* ma version  
Bag::Bag(void) 
{
   BagList list = {
      NULL, new Object() };
      
   head = &list;
}
*/
//corrige
Bag::Bag(void) : Object()
{
    return; //ne fait rien d'autre que d'appeler Object::Object
}

// corrige
Bag::~Bag(void)
{
   Bag::BagList *temp = head; // detruit le liste d'objet
   while (temp != NULL)
   {
      temp=temp->next;
      delete head;
      head = temp;   
   }
   return;
}


void Bag::print(void) {
   Bag::BagList *temp=head;
   cout << "Sac numero " << handle() << "." << endl;
   cout << "   Contenu:" << endl;
   if (temp==NULL) {
      cout << "empty bag\n"; 
      return;
   }
   while (temp!=NULL) {
      cout << "Sac numero " << temp->ptr->handle() << "." << endl; //rem: en fait on pourrait appeler temp->print()
      temp = temp->next;      
   }
   return;
}
  
bool Bag::has(unsigned long int h) const 
{
   Bag::BagList *temp = head;
   if (temp==NULL) {
      cout << "empty bag\n";
      return false;   
   }
   while (temp!=NULL) {
      if (temp->ptr->handle() == h)
         return true;
   }
   return false;
}

void Bag::add(Object &ob) 
{
   Bag::BagList *temp = new BagList;
      temp->ptr=&ob;
      temp->next=head;
      head=temp;
   
   return;
}

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


Bag monSac;

int main() {
   MonObjet a,b,c;
   
   monSac.add(a);
   monSac.add(b);
   monSac.add(c);
   
   monSac.print();
   return 0;   
}

