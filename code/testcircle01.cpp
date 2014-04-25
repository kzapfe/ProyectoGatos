//Probemos hacer gato estados bidimiensionales aleatoreos sobre la camada de un o.a.

#include <iostream>
#include <fstream>
#include <vector>
#include <armadillo>


#include "CatStates01.hpp"
#include "PopulateSphericalShells01.hpp"


using namespace std;
using namespace arma;
 

int main(){
  //No mames que bien trabajas despues de hacer ejercicio 
  // Y SIN INTERNET. Recuerda eso SIEMPRE.

  vector   <simplectic> x,y;
  //El numero de estados coherentes. El numero de gatos es 
  // la combinacion de pares posibles.
  const int maximumgauss=70;
  const double radiointeres=14;
  int resol=600;

  mat centros;
  centros=Populate2BallShell(2000, maximumgauss);
  //Los centros, otra vez para variar.
 

  double DobleWignerFunction;
  gsl_complex DobleWeylFunction;
 
  //iniocializar el estado poblado
  //parece ser que vector no funciona aqui bien
  CoherentState *CentroX ;
  CatState * Gatos;

  CentroX=new CoherentState[maximumgauss]; 
  Gatos=new CatState[maximumgauss*(maximumgauss-1)/2];

  for(int i=0; i<maximumgauss; i++){

    CentroX[i].SetCentre(centros(i,0), centros(i,1));
  
  };

  ofstream Centros;
  Centros.open("CentrosGaussianas.dat");

  Centros<<centros<<endl;

  int cuentainterferencias=0;

  for(int i=0; i<maximumgauss-1; i++){
    for(int j=i+1; j<maximumgauss; j++){
      Gatos[cuentainterferencias].SetGaussians(CentroX[i], CentroX[j]);
      cuentainterferencias++;

    }

  }

  //ahora toca calcular el valor de las chivas en la funcion:
  //Corte en el plano de las q y en el de las y.

  ofstream CorteEnX;
  
  CorteEnX.open("WignerCorteX.dat");  
  
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
      for(int i=0; i<maximumgauss; i++){
	
	WignerFunctionCoherent+=CentroX[i].CentreRepresentation(falsax);
      
      }

      WignerFunctionCoherent/=(double)maximumgauss;



      for(int i=0; i<cuentainterferencias; i++){
	WignerFunctionInterference+=
	  Gatos[i].InterferenceTermCentreRepresentation(falsax);

      }
      WignerFunctionInterference/=(double)cuentainterferencias;
  
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

