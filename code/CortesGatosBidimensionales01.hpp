//Veamos si puedes ordanar un poco tu programa sacando el corte afuera del main

/*Rutinas que obtienen cortes  
  bidimensionales de Wigner y Weyl 
y los escriben al disco duro*/

#include <iostream>
#include <fstream>
//#include "CatStates01.hpp"
//#include "PopulateSphericalShells01.hpp"

using namespace std;

void WignerSection(CoherentState CenX[], CoherentState CenY[],
		   CatState CatX[], CatState CatY[], 
		   int maxgauss, int maxcats,
		   bool selectx, bool selecty,
		   string nombrearchivo){

  ofstream Archivo;
  //Funcion chingona que convierte un string a c_str que funciona para open.
  Archivo.open(nombrearchivo.c_str()); 
  Archivo<<"#Testing" <<endl;
  
  /*Tabla Booleana */
  //No te interesan dos cortes: los cruzados
  /* 0 0 -> qx qy 
     0 1 -> qx px
     1 0 -> qy py
     1 1 -> px py */
    
  
   double auxxq=0.00, auxxp=0.00, auxyq=0.00, auxyp=0.00;
   double coordenadauno, coordenadados;
   simplectic falsax(auxxq, auxxp), falsay(auxyq,auxyp);
  
   clock_t empieza,final;
   double tiemponecesario=0.0;
   
    if((!selectx)&&(!selecty))
	    { //Corte en Q	 
	      Archivo<<"# Corte en qx, qy"<<endl;
	    }
	  else if((selectx)&&(!selecty))
	    { //Corte en X	      
	      Archivo<<"# Corte en qx, px"<<endl;
	    }
	  else if ((!selectx)&&(selecty))
	    { // Corte en Y
	      Archivo<<"# Corte en qy,py"<<endl;
	    }
	  else if((selectx)&&(selecty))
	    { //Corte en P     
	      Archivo<<"# Corte en px, py"<<endl;
	    }
	 


   empieza=clock();
    
      for(int n=-resol; n<resol; n++){
	for(int m=-resol; m<resol; m++){
	  //partes chatas
	  
	  double WignerFunctionCoherent=0.000;
	  double WignerFunctionInterference=0.000;
	  double WignerTotal=0.000;
	  	
	  coordenadauno=magnituddeinteres*(double)n/(double)resol;     
	  coordenadados=magnituddeinteres*(double)m/(double)resol;      
           
	  if((!selectx)&&(!selecty))
	    { //Corte en Q
	      falsax.q=coordenadauno;
	      falsay.q=coordenadados;
	      
	    }
	  else if((selectx)&&(!selecty))
	    { //Corte en X
	      falsax.q=coordenadauno;
	      falsax.p=coordenadados;
	      
	    }
	  else if ((!selectx)&&(selecty))
	    { // Corte en Y
	      falsay.q=coordenadauno;
	      falsay.p=coordenadados;
	      
	    }
	  else if((selectx)&&(selecty))
	    { //Corte en P
	      falsax.p=coordenadauno;
	      falsay.p=coordenadados;
	      
	    }
	  
	  //lets get to business
      
	  //primero las puntas gaussianas	  
	  for(int i=0; i<maxgauss; i++){	
	    WignerFunctionCoherent+=CenX[i].CentreRepresentation(falsax)*
	      CenY[i].CentreRepresentation(falsay);	
	  }
	  //Luego las interferencias
	  for(int i=0; i<maxcats; i++){
	    WignerFunctionInterference+=
	      CatX[i].InterferenceTermCentreRepresentation(falsax)
	      *CatY[i].InterferenceTermCentreRepresentation(falsay);
	    
	  }
      
        
	  WignerFunctionCoherent/=(double)maxgauss;
	  WignerFunctionInterference/=(double)maxcats;
	  WignerTotal=WignerFunctionInterference+WignerFunctionCoherent;
	  
	  Archivo<<coordenadauno<<"\t"<<coordenadados<<"\t"		  
		  <<WignerFunctionCoherent<<"\t"
		  <<WignerFunctionInterference<<"\t"
		  <<WignerTotal
		  <<endl;
      
	}
    
	Archivo<<endl;
	
      } //Termina el dibujo (sobre el segundo ciclo de resoluciones)
  
      Archivo.close();
      final=clock();
      
      tiemponecesario=double(final-empieza)/CLOCKS_PER_SEC;
      cout<<"Me tarde "<<tiemponecesario
	  << "s en calcular el corte "<< 
	selectx << selecty<<
	" en Wigner. "<<endl;
      //Aqui termina su trabajo el hilo 0



      return;

};

void WeylSection(CoherentState CenX[], CoherentState CenY[],
		   CatState CatX[], CatState CatY[], 
		   int maxgauss, int maxcats,
		   bool selectx, bool selecty,
		   string nombrearchivo){

  ofstream Archivo;
  //Funcion chingona que convierte un string a c_str que funciona para open.
  Archivo.open(nombrearchivo.c_str()); 
  Archivo<<"#Testing" <<endl;
  
  //mu es la SFT de x, xhi de y
  /*Tabla Booleana */
  //No te interesan dos cortes: los cruzados
  /* 0 0 -> muq xhiq 
     0 1 -> muq mup
     1 0 -> xhiq xhip
     1 1 -> mup xhip */
    
   double auxmuq=0.00, auxmup=0.00, auxxhiq=0.00, auxxhip=0.00;
   double coordenadauno, coordenadados;
   simplectic mu(auxmuq, auxmup), xhi(auxxhiq,auxxhip);
  
   clock_t empieza,final;
   double tiemponecesario=0.0;
   
    if((!selectx)&&(!selecty))
	    { //Corte en Q	 
	      Archivo<<"# Corte en muq, xhiq"<<endl;
	    }
	  else if((selectx)&&(!selecty))
	    { //Corte en X	      
	      Archivo<<"# Corte en muq, mup"<<endl;
	    }
	  else if ((!selectx)&&(selecty))
	    { // Corte en Y
	      Archivo<<"# Corte en xhiq,xhip"<<endl;
	    }
	  else if((selectx)&&(selecty))
	    { //Corte en P     
	      Archivo<<"# Corte en mup, xhip"<<endl;
	    }
	 


   empieza=clock();
    
      for(int n=-resol; n<resol; n++){
	for(int m=-resol; m<resol; m++){
	  //partes chatas
	  
	  gsl_complex WeylFunctionCoherent;
	  gsl_complex WeylFunctionInterference;
	  gsl_complex WeylFunctionTotal;
	 
	  //El corte de Weyl tiene la escala AL REVEZ , aprox
	  coordenadauno=(1./magnituddeinteres)*(double)n/(double)resol;     
	  coordenadados=(1./magnituddeinteres)*(double)m/(double)resol;      
           
	  if((!selectx)&&(!selecty))
	    { //Corte en Q
	      mu.q=coordenadauno;
	      xhi.q=coordenadados;
	      
	    }
	  else if((selectx)&&(!selecty))
	    { //Corte en Mu
	      mu.q=coordenadauno;
	      mu.p=coordenadados;
	      
	    }
	  else if ((!selectx)&&(selecty))
	    { // Corte en Xhi
	      xhi.q=coordenadauno;
	      xhi.p=coordenadados;
	      
	    }
	  else if((selectx)&&(selecty))
	    { //Corte en P
	      mu.p=coordenadauno;
	      xhi.p=coordenadados;
	      
	    }
	  
	  //lets get to business
      
	  //primero las puntas gaussianas	  
	 
	  for(int i=0; i<maxgauss; i++){	
	    WeylFunctionCoherent=
	      gsl_complex_add(WeylFunctionCoherent,					
			      gsl_complex_mul(
					      CenX[i].ChordRepresentation(mu),
					      CenY[i].ChordRepresentation(xhi)));
	  }
	  //Luego las interferencias
	  for(int i=0; i<maxcats; i++){
	    //We may need sintactic Suger here
	    WeylFunctionInterference=
	      gsl_complex_add(WeylFunctionInterference,
			      gsl_complex_mul(
					      CatX[i].InterferenceTermChordRepresentation(mu),
					      CatY[i].InterferenceTermChordRepresentation(xhi)));
	    
	  }
      
        
	  WeylFunctionCoherent=gsl_complex_div_real(WeylFunctionCoherent,
						    (double)maxgauss);
	  WeylFunctionInterference=gsl_complex_div_real(WeylFunctionInterference,
							(double)maxcats);
	  WeylFunctionTotal=gsl_complex_add(WeylFunctionCoherent, 
					    WeylFunctionInterference);
	  
	  Archivo<<coordenadauno<<"\t"<<coordenadados<<"\t"		  
		 <<GSL_REAL(WeylFunctionCoherent)<<"\t"
		 <<GSL_IMAG(WeylFunctionCoherent)<<"\t"
		 <<GSL_REAL(WeylFunctionInterference)<<"\t"
		 <<GSL_IMAG(WeylFunctionInterference)<<"\t"
		 <<GSL_REAL(WeylFunctionTotal)<<"\t"
		 <<GSL_IMAG(WeylFunctionTotal)<<endl;

	   
      
	}
    
	Archivo<<endl;
	
      } //Termina el dibujo (sobre el segundo ciclo de resoluciones)
  
      Archivo.close();
      final=clock();
      
      tiemponecesario=double(final-empieza)/CLOCKS_PER_SEC;
      cout<<"Me tarde "<<tiemponecesario
	  << "s en calcular el corte "<< 
	selectx << selecty<<
	" en Wigner. "<<endl;
      //Aqui termina su trabajo el hilo 0



      return;

};
