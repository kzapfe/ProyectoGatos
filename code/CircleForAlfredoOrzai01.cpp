//Probemos hacer gato estados bidimiensionales aleatoreos sobre la camada de un o.a.

#include <iostream>
#include <fstream>
#include <vector>


#include "CatStates01.hpp"

using namespace std;

 
int main(){
  //No mames que bien trabajas despues de hacer ejercicio 
  // Y SIN INTERNET. Recuerda eso SIEMPRE.

  vector   <simplectic> x,y;
  //El numero de estados coherentes. El numero de gatos es 
  // la combinacion de pares posibles.

  const double radiointeres=8.*sqrt(hbar);

  int resol=200;
  double energiaoscilador=10;
  double Qmay=4.*sqrt(hbar);
  
  double DobleWignerFunction;
  gsl_complex DobleWeylFunction;
  
  CoherentState CentroUno, CentroDos ;
  

  CentroUno.SetCentre(Qmay, 0.00);
  CentroDos.SetCentre(-Qmay, 0.00);
  
  CatState  Gatocentral(CentroUno, CentroDos, 1.0, -1.00);
  
  ofstream Centros;
  Centros.open("CentrosGaussianas.dat");

  
  //ahora toca calcular el valor de las chivas en la funcion:
  //Corte en el plano de las q y en el de las y.

  ofstream CorteEnX;
  
  CorteEnX.open("WignerPunchCorteX.dat");  
  
  double auxxq, auxxp;


  for(int n=-resol; n<resol; n++){
    for(int m=-resol; m<resol; m++){
      //partes chatas
      
      double WignerFunctionCoherent=0.000;
      double WignerFunctionInterference=0.000;
      double WignerTotal=0.000;

      auxxq=radiointeres*(double)n/(double)resol;
      auxxp=radiointeres*(double)m/(double)resol;      
      simplectic falsax(auxxq, auxxp);

      //primero las puntas gaussianas

	
	WignerFunctionCoherent=CentroUno.CentreRepresentation(falsax)
	+CentroDos.CentreRepresentation(falsax);
           
       	WignerFunctionInterference=
	  	  Gatocentral.InterferenceTermCentreRepresentation(falsax);
      
	
	WignerTotal=WignerFunctionInterference+WignerFunctionCoherent;

	CorteEnX<<auxxq<<"\t"<<auxxp<<"\t"<<WignerFunctionCoherent<<"\t"
	      <<WignerFunctionInterference<<"\t"<<WignerTotal
	      <<endl;
      
      
    }
    
	CorteEnX<<endl;
    
  }
  
  CorteEnX.close();
  

  
  return 0;
 
}

