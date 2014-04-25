//Probar el hpp de estados coherentes

#include "CatStates01.hpp"
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
  CatState Gato1;
  CatState Gato2(MedioGato, OtroMedioGato);


  //cout<<OtroMedioGato.centre.q<<endl;
  //cout<<OtroMedioGato.centre.p<<endl;
  
  simplectic yuno(1,0), ydos(0,0);
  simplectic zuno, zdos;

  cout<<"Antes"<<endl;
  cout<<yuno.q<<"\t"<<yuno.p<<endl;
  cout<<ydos.q<<"\t"<<ydos.p<<endl;
  
  ydos=yuno;

  cout<<"Despues"<<endl;
  cout<<yuno.q<<"\t"<<yuno.p<<endl;
  cout<<ydos.q<<"\t"<<ydos.p<<endl;

  zuno=yuno;
  zdos.q=0.5;
  zdos.p=-0.5;

  yuno=zuno+zdos;

  cout<<"Despues"<<endl;
  cout<<zuno.q<<"\t"<<zuno.p<<endl;
  cout<<zdos.q<<"\t"<<zdos.p<<endl;
  cout<<yuno.q<<"\t"<<yuno.p<<endl;

  CoherentState MedioGato3(yuno);


  CatState OtroGato(MedioGato3, OtroMedioGato);

  ofstream cuerdasuno, cuerdasdos;
  cuerdasuno.open("GatosWeyl01.dat");
  cuerdasdos.open("GatosWelr02.dat");

  double poq;
  simplectic x;
  gsl_complex testuno, testdos;
  int resol=500; 

  for(int i=-resol; i< resol; i++){
    for(int j=-resol; j< resol; j++){
    
      x.q=(double)i/100.0;
      x.p=(double)j/100.0;

      testuno=Gato2.InterferenceTermChordRepresentation(x);
      testdos=OtroGato.InterferenceTermChordRepresentation(x);

      cuerdasuno<<x.q<<"\t"<< x.p << "\t"  << 
	"\t" << GSL_REAL(testuno) 
		<< endl;
   
      cuerdasdos<<x.q<<"\t"<< x.p << "\t"  << 
	"\t" << GSL_IMAG(testdos) 
		<< endl;
   
      


    }
    cuerdasuno <<endl;
    cuerdasdos <<endl;

  }

  
  return 0;
 
}

