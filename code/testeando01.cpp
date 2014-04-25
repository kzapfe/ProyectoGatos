//Probar el hpp de estados coherentes

#include "CoherentState01.hpp"
#include <iostream>
#include <fstream>

using namespace std;


int main(){
  //No mames que bien trabajas despues de hacer ejercicio 
  // Y SIN INTERNET. Recuerda eso SIEMPRE.

  CoherentState MedioGato;
  simplectic xuno;
  xuno.q=2;
  xuno.p=1.5;
  CoherentState OtroMedioGato(xuno);

  //cout<<OtroMedioGato.centre.q<<endl;
  //cout<<OtroMedioGato.centre.p<<endl;
  
  

  double poq;
  simplectic x;
  double testeando;
  int resol=500; 

  for(int i=-resol; i< resol; i++){
    for(int j=-resol; j< resol; j++){
    
      x.q=(double)i/100.0;
      x.p=(double)j/100.0;

      testeando=OtroMedioGato.CentreRepresentation(x);
      cout<<x.q<<"\t"<< x.p << "\t" << testeando << endl;
    }
    cout<<endl;

  }
  return 0;
 
}

