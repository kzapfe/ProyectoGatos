//Probar el hpp de estados coherentes

#include "CoherentState01.hpp"
#include <iostream>
#include <fstream>

using namespace std;


int main(){
  //No mames que bien trabajas despues de hacer ejercicio 
  // Y SIN INTERNET. Recuerda eso SIEMPRE.

  CoherentState MedioGato;
  simplectic xuno, xdos;
  xuno.q=2;
  xuno.p=1.5;
  CoherentState OtroMedioGato(xuno);
  CoherentState Otroestado(OtroMedioGato);

  swap(xuno, xdos);
  simplectic xtres;
  xtres=xdos;

  //cout<<OtroMedioGato.centre.q<<endl;
  //cout<<OtroMedioGato.centre.p<<endl;
  
  ofstream cuerdasuno, cuerdasdos;
  cuerdasuno.open("CoherentCuerdas01.dat");
  cuerdasdos.open("CoherentCuerdas02.dat");

  double poq;
  simplectic x;
  gsl_complex testuno, testdos;
  int resol=500; 

  for(int i=-resol; i< resol; i++){
    for(int j=-resol; j< resol; j++){
    
      x.q=(double)i/100.0;
      x.p=(double)j/100.0;

      testuno=MedioGato.ChordRepresentation(x);
      testdos=OtroMedioGato.ChordRepresentation(x);

      cuerdasuno<<x.q<<"\t"<< x.p << "\t"  << 
	"\t" << GSL_REAL(testuno) << "\t" << GSL_IMAG(testuno) 
		<< endl;
   
      cuerdasdos<<x.q<<"\t"<< x.p << "\t"  << 
	"\t" << GSL_REAL(testdos) << "\t" << GSL_IMAG(testdos) 
		<< endl;
   
      


    }
    cuerdasuno <<endl;
    cuerdasdos <<endl;

  }
  return 0;
 
}

