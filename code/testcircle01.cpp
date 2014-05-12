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
  const int maximumgauss=7;
  const double radiointeres=16;
  int resol=300;

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

  ofstream Centros, Cuerdas;
  Centros.open("CentrosCirculo7.dat");

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
  
  CorteEnX.open("WignerCirculo7X.dat");  
  
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

      //Luego las interferencias

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
  
  ofstream CorteEnXhi;
  
  CorteEnXhi.open("WeylCirculo7Xhi.dat");  

 
  for(int n=-resol; n<resol; n++){
    for(int m=-resol; m<resol; m++){
      //partes chatas
      
      
      gsl_complex WeylFunctionCoherent;
      gsl_complex WeylFunctionInterference;
      gsl_complex WeylFunctionTotal;
	 
      auxxq=radiointeres*(double)n/(double)resol;
      auxxp=radiointeres*(double)m/(double)resol;      
      simplectic falsax(auxxq, auxxp);

      //primero las puntas gaussianas
      for(int i=0; i<maximumgauss; i++){
	
	WeylFunctionCoherent=
	  gsl_complex_add(WeylFunctionCoherent,
			  CentroX[i].ChordRepresentation(falsax));
      	
      }


      for(int i=0; i<cuentainterferencias; i++){
	//We may need sintactic Suger here
	WeylFunctionInterference=
	  gsl_complex_add(WeylFunctionInterference,
			  Gatos[i].InterferenceTermChordRepresentation(falsax));	
      }
      
      
	  //Esto esta siendo usado suponiendo que en el valor dela funcion todo este bien normalizado
      WeylFunctionCoherent=gsl_complex_div_real(WeylFunctionCoherent,
						(double)maximumgauss);
      WeylFunctionInterference=gsl_complex_div_real(WeylFunctionInterference,
						    (double)cuentainterferencias);
      WeylFunctionTotal=gsl_complex_add(WeylFunctionCoherent, 
					    WeylFunctionInterference);
      
      CorteEnXhi<<falsax.q<<"\t"<<falsax.p<<"\t"		  
		<<GSL_REAL(WeylFunctionCoherent)<<"\t"
		<<GSL_IMAG(WeylFunctionCoherent)<<"\t"
		<<GSL_REAL(WeylFunctionInterference)<<"\t"
		<<GSL_IMAG(WeylFunctionInterference)<<"\t"
		<<GSL_REAL(WeylFunctionTotal)<<"\t"
		<<GSL_IMAG(WeylFunctionTotal)<<endl;
      
	   	  
      
    }
    
    CorteEnXhi<<endl;
    
  }
  
        
      
  
  CorteEnXhi.close();
 
  
  return 0;
 
}

