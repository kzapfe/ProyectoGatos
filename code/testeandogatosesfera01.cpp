//Probemos hacer gato estados cuadridimensionales.

#include <iostream>
#include <fstream>
#include <vector>
#include <armadillo>

#include <ctime> //para ver cuanto se tarda en cada punto

//#include <gsl/gsl_rng.h> //Inicializa el gsl random
//#include <gsl/gsl_randist.h> //las distribuciones de gsl


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
  const int maximumgauss=200;
  int resol=100;
  double magnituddeinteres=3.0;

  mat centros;
  centros=Populate4BallShell(100, maximumgauss);
  //Los centros, otra vez para variar.
 


  double DobleWignerFunction;
  gsl_complex DobleWeylFunction;
 

  //iniocilizar el estado poblado
  //parece ser que vector no funciona aqui bien
  CoherentState *CentroX, *CentroY;
  CatState *GatosX, *GatosY;

  CentroX=new CoherentState[maximumgauss];
  CentroY=new CoherentState[maximumgauss];

  GatosX=new CatState[maximumgauss*(maximumgauss-1)/2];
  GatosY=new CatState[maximumgauss*(maximumgauss-1)/2];

  clock_t empieza,final;
  double tiemponecesario=0.0;
  empieza=clock();

  for(int i=0; i<maximumgauss; i++){
    //Selecciona centros en los subespacios X y Y
    CentroX[i].SetCentre(centros(i,0), centros(i,1));
    CentroY[i].SetCentre(centros(i,2), centros(i,3));


  };

  int cuentainterferencias=0;

  for(int i=0; i<maximumgauss-1; i++){
    for(int j=i+1; j<maximumgauss; j++){
      //las interferencias son Indepèndientes
      GatosX[cuentainterferencias].SetGaussians(CentroX[i], CentroX[j]);
      GatosY[cuentainterferencias].SetGaussians(CentroY[i], CentroY[j]);      
      cuentainterferencias++;      
    }
  }

  final=clock();

  tiemponecesario=double(final-empieza)/CLOCKS_PER_SEC;
  cout<<"Me tarde "<<tiemponecesario<< "s en llenar de gatos la 4 esfera. "<<endl;

  //ahora toca calcular el valor de las chivas en la funcion:
  //Corte en el plano de las q y en el de las y.

  ofstream CorteEnQ, CorteEnY;
  
  CorteEnQ.open("WignerCorteQ.dat");
  CorteEnY.open("WignerCorteY.dat");
  
  double auxxq, auxxp, auxyq, auxyp;
 
  empieza=clock();

  for(int n=-resol; n<resol; n++){
    for(int m=-resol; m<resol; m++){
      //partes chatas
      
      double WignerFunctionCoherent=0.000;
      double WignerFunctionInterference=0.000;
      double WignerTotal=0.000;


      auxxq=magnituddeinteres*(double)n/(double)resol;     
      auxxp=magnituddeinteres*(double)m/(double)resol;      

      
      simplectic falsax(auxxq, auxxp), falsay(0.00,0.00);
      
      //primero las puntas gaussianas
      //Resulta que si no usas la cabezota, no parece ser que las
      //dos cortes se puedan ver juntos.
      for(int i=0; i<maximumgauss; i++){
	
	WignerFunctionCoherent+=CentroX[i].CentreRepresentation(falsax)*
	  CentroY[i].CentreRepresentation(falsay);
      
      }

      for(int i=0; i<cuentainterferencias; i++){
	WignerFunctionInterference+=
	  GatosX[i].InterferenceTermCentreRepresentation(falsax)
	  *GatosY[i].InterferenceTermCentreRepresentation(falsax);

      }



      WignerFunctionCoherent/=(double)maximumgauss;
      WignerFunctionInterference/=(double)cuentainterferencias;
      WignerTotal=WignerFunctionInterference+WignerFunctionCoherent;


      CorteEnQ<<auxxq<<"\t"<<auxxp<<"\t"
	      <<auxyq<<"\t"<<auxyp<<"\t"
	      <<WignerFunctionCoherent<<"\t"
	      <<WignerFunctionInterference<<"\t"
	      <<WignerTotal
	      <<endl;
      
      
    }
    
    CorteEnQ<<endl;
    
  }
  
  CorteEnQ.close();
  final=clock();

  tiemponecesario=double(final-empieza)/CLOCKS_PER_SEC;
  cout<<"Me tarde "<<tiemponecesario
      << "s en calcular el primer corte en Wigner. "<<endl;


  
  return 0;
 
}

