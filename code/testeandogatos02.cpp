//Probemos hacer gato estados cuadridimensionales.

#include <iostream>
#include <fstream>

#include "CatStates01.hpp"


using namespace std;


int main(){
  //No mames que bien trabajas despues de hacer ejercicio 
  // Y SIN INTERNET. Recuerda eso SIEMPRE.

  
  simplectic x,y;
  int resol=400;

  double DobleWignerFunction;
  gsl_complex DobleWeylFunction;
 

  x.q=0.5;
  x.p=0.5;
  y.q=-0.5;
  y.p=-0.5;

  CoherentState GatoX(x), GatoY(y);
 
  CatState GatoTotal(GatoX,GatoY);

  ofstream DobleWigner;
  DobleWigner.open("Testing2dimCatsCentre01.dat");
  ofstream DobleWeyl;
  DobleWeyl.open("Testing2dimCatsChords01.dat");

  simplectic puntoeval1,puntoeval2;

  
  for(int i=-resol; i< resol; i++){
    for(int j=-resol; j< resol; j++){
    
//  for(int i=12; i< 14; i++){
//    for(int j=12; j< 14; j++){
    
  
      puntoeval1.q=3*(double)i/resol;
      puntoeval1.p=3*(double)j/resol;
      
      DobleWignerFunction=GatoTotal.CentreRepresentation(puntoeval1);
      
      DobleWeylFunction=GatoTotal.ChordRepresentation(puntoeval1);

      DobleWigner<<puntoeval1.q<<"\t"<< puntoeval1.p << "\t"  << 
	"\t" << DobleWignerFunction 
		<< endl;
   
      DobleWeyl<<puntoeval1.q<<"\t"<< puntoeval1.p << "\t"  << 
	"\t" << GSL_REAL(DobleWeylFunction) <<
	"\t" << GSL_IMAG(DobleWeylFunction) 
	       << endl;
   
      


    }
    DobleWigner <<endl;
    DobleWeyl <<endl;

  }
  
  
  return 0;
 
}

