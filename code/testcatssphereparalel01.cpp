//Probemos hacer gato estados cuadridimensionales.

#include <iostream>
#include <fstream>
#include <vector>
#include <armadillo>

#include <ctime> //para ver cuanto se tarda en cada punto
#include <omp.h> //para ver cuanto se tarda en cada punto



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
  const int maximumgauss=400;
  int resol=200;
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

  ofstream CorteEnQ, CorteEnY, CorteEnP;
  
  CorteEnQ.open("WignerCorteQ4S02.dat");
  CorteEnY.open("WignerCorteY4S02.dat");
  CorteEnP.open("WignerCorteP4S02.dat");
  
  
#pragma omp parallel num_threads(3)
  {

  
    int hilo;

    hilo=omp_get_thread_num();
    if(hilo==0){ 
      cout<<"Hilo 0, haciendo Corte en Q"<<endl;

      double auxxq, auxxp, auxyq, auxyp;
    
      empieza=clock();
    
      for(int n=-resol; n<resol; n++){
	for(int m=-resol; m<resol; m++){
	  //partes chatas
	  
	  double WignerFunctionCoherent=0.000;
	  double WignerFunctionInterference=0.000;
	  double WignerTotal=0.000;
	  	
	  auxxq=magnituddeinteres*(double)n/(double)resol;     
	  auxyq=magnituddeinteres*(double)m/(double)resol;      
           
	  simplectic falsax(auxxq, 0.00), falsay(auxyq,0.00);
      
	  //primero las puntas gaussianas	  
	  for(int i=0; i<maximumgauss; i++){	
	    WignerFunctionCoherent+=CentroX[i].CentreRepresentation(falsax)*
	      CentroY[i].CentreRepresentation(falsay);	
	  }
	  //Luego las interferencias
	  for(int i=0; i<cuentainterferencias; i++){
	    WignerFunctionInterference+=
	      GatosX[i].InterferenceTermCentreRepresentation(falsax)
	      *GatosY[i].InterferenceTermCentreRepresentation(falsay);
	    
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
	
      } //Termina el dibujo (sobre el segundo ciclo de resoluciones)
  
      CorteEnQ.close();
      final=clock();
      
      tiemponecesario=double(final-empieza)/CLOCKS_PER_SEC;
      cout<<"Me tarde "<<tiemponecesario
	<< "s en calcular el Q-corte en Wigner. "<<endl;
      //Aqui termina su trabajo el hilo 0

    }else if(hilo==1){
      //Comienza su trabajo el hilo 1

      cout<<"Hilo 0, haciendo Corte en Y"<<endl;

      double auxxq, auxxp, auxyq, auxyp;
    
      empieza=clock();
    
      for(int n=-resol; n<resol; n++){
	for(int m=-resol; m<resol; m++){
	  //partes chatas
	  
	  double WignerFunctionCoherent=0.000;
	  double WignerFunctionInterference=0.000;
	  double WignerTotal=0.000;
	  	
	  auxyq=magnituddeinteres*(double)n/(double)resol;     
	  auxyp=magnituddeinteres*(double)m/(double)resol;      
           
	  simplectic falsax(0.000, 0.000), falsay(auxyq,auxyp);
      
	  //primero las puntas gaussianas	  
	  for(int i=0; i<maximumgauss; i++){	
	    WignerFunctionCoherent+=CentroX[i].CentreRepresentation(falsax)*
	      CentroY[i].CentreRepresentation(falsay);	
	  }
	  //Luego las interferencias
	  for(int i=0; i<cuentainterferencias; i++){
	    WignerFunctionInterference+=
	      GatosX[i].InterferenceTermCentreRepresentation(falsax)
	      *GatosY[i].InterferenceTermCentreRepresentation(falsay);
	    
	  }
      
        
	  WignerFunctionCoherent/=(double)maximumgauss;
	  WignerFunctionInterference/=(double)cuentainterferencias;
	  WignerTotal=WignerFunctionInterference+WignerFunctionCoherent;
	  
	  CorteEnY<<auxxq<<"\t"<<auxxp<<"\t"
		  <<auxyq<<"\t"<<auxyp<<"\t"
		  <<WignerFunctionCoherent<<"\t"
		  <<WignerFunctionInterference<<"\t"
		  <<WignerTotal
		  <<endl;
      
	}
    
	CorteEnY<<endl;
	
      } //Termina el dibujo (sobre el segundo ciclo de resoluciones)
  
      CorteEnY.close();
      final=clock();
      
      tiemponecesario=double(final-empieza)/CLOCKS_PER_SEC;
      cout<<"Me tarde "<<tiemponecesario
	<< "s en calcular el  corte Y en Wigner. "<<endl;

    }else if(hilo==2){
      cout<<"Hilo 2, haciendo Corte en P"<<endl;

      double auxxq, auxxp, auxyq, auxyp;
    
      empieza=clock();
    
      for(int n=-resol; n<resol; n++){
	for(int m=-resol; m<resol; m++){
	  //partes chatas
	  
	  double WignerFunctionCoherent=0.000;
	  double WignerFunctionInterference=0.000;
	  double WignerTotal=0.000;
	  	
	  auxxp=magnituddeinteres*(double)n/(double)resol;     
	  auxyp=magnituddeinteres*(double)m/(double)resol;      
           
	  simplectic falsax(0.000, auxxp), falsay(0.00,auxyp);
      
	  //primero las puntas gaussianas	  
	  for(int i=0; i<maximumgauss; i++){	
	    WignerFunctionCoherent+=CentroX[i].CentreRepresentation(falsax)*
	      CentroY[i].CentreRepresentation(falsay);	
	  }
	  //Luego las interferencias
	  for(int i=0; i<cuentainterferencias; i++){
	    WignerFunctionInterference+=
	      GatosX[i].InterferenceTermCentreRepresentation(falsax)
	      *GatosY[i].InterferenceTermCentreRepresentation(falsay);
	    
	  }
      
        
	  WignerFunctionCoherent/=(double)maximumgauss;
	  WignerFunctionInterference/=(double)cuentainterferencias;
	  WignerTotal=WignerFunctionInterference+WignerFunctionCoherent;
	  
	  CorteEnP<<auxxq<<"\t"<<auxxp<<"\t"
		  <<auxyq<<"\t"<<auxyp<<"\t"
		  <<WignerFunctionCoherent<<"\t"
		  <<WignerFunctionInterference<<"\t"
		  <<WignerTotal
		  <<endl;
      
	}
    
	CorteEnP<<endl;
	
      } //Termina el dibujo (sobre el segundo ciclo de resoluciones)
  
      CorteEnP.close();
      final=clock();
      
      tiemponecesario=double(final-empieza)/CLOCKS_PER_SEC;
      cout<<"Me tarde "<<tiemponecesario
	<< "s en calcular el corte P en Wigner. "<<endl;

      
    }

  } //Cierra el pragma paralel omp loop

  
  return 0;
 
}

