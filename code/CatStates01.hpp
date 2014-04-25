
#ifndef __CatState__
#define __CatState__

#include "CoherentState01.hpp"

class CatState{
  /*Clase de Estados Gatos. 
    Requiere dos estados coherentes para funcionar */

public: 

 
  gsl_complex PesoUno, PesoDos;
  CoherentState  Uno, Dos;
  simplectic PuntaUno, PuntaDos;

  CatState(){
    //Todo se va a valores por omision
     //Por omision vamos a hacer uno siempre los 
    //pesos de las dos mitades coherentes
    GSL_SET_COMPLEX(&PesoUno, 1.0, 0.0);
    GSL_SET_COMPLEX(&PesoDos, 1.0, 0.0);
  }

  CatState(CoherentState & EstadoUno, CoherentState & EstadoDos){
    //Por omision vamos a hacer uno siempre los 
    //pesos de las dos mitades coherentes
    GSL_SET_COMPLEX(&PesoUno, 1.0, 0.0);
    GSL_SET_COMPLEX(&PesoDos, 1.0, 0.0);
    
    this->Uno=EstadoUno;
    this->Dos=EstadoDos;

  }


  CatState(CoherentState & EstadoUno, CoherentState & EstadoDos, 
	   double contribucionuno, double contribuciondos){
    GSL_SET_COMPLEX(&PesoUno, contribucionuno, 0.0);
    GSL_SET_COMPLEX(&PesoDos, contribuciondos, 0.0);
    this->Uno=EstadoUno;
    this->Dos=EstadoDos;
    PesoUno=gsl_complex_rect(contribucionuno,0.00);
    PesoDos=gsl_complex_rect(contribuciondos,0.00);
  }

  CatState(CoherentState & EstadoUno, CoherentState & EstadoDos, 
	   gsl_complex contribucionuno, gsl_complex contribuciondos){
    this->PesoUno=contribucionuno;
    this->PesoDos=contribuciondos;
    this->Uno=EstadoUno;
    this->Dos=EstadoDos;

  }

  void SetGaussians(CoherentState &EstadoUno, CoherentState &EstadoDos){
    this->Uno=EstadoUno;
    this->Dos=EstadoDos;
  };
    

  friend void swap(CatState& Primero, CatState& Segundo){
    //intercambia dos CoherentStates, depende de <algorithm>
    using std::swap;
    swap(Primero.Uno, Segundo.Uno);    
    swap(Primero.Dos, Segundo.Dos);    
    swap(Primero.PesoUno, Segundo.PesoUno);    
    swap(Primero.PesoDos, Segundo.PesoDos);    
  };

  CatState& operator=(CatState RightHandSide){
    swap(*this, RightHandSide);
    return *this;
  }


  double InterferenceTermCentreRepresentation(simplectic& punto){
    /*la funcion devuelve SOLAMENTE el termino de interferencia
      de la funcion de Wigner correspondiente al punto,
    la Notacion sigue el articulo de Nicacio y Raul Vallejos.
    La transformacion metaplectica se reduce a la unidad.*/

    double magnitud;
    simplectic zeta, eta;
    //la magnitud lleva un factor que no aparece en el art de Raul.
    // ya que tu consideras masas y vibraciones arbitrarias.
    double FactorCorrectivopormasayomega;
    double exponentegaussiana, faseinterferencia;
    double result;

    zeta=Uno.centre-Dos.centre;
    eta=(Uno.centre+Dos.centre);
    eta.q=eta.q/2.0;
    eta.p=eta.p/2.0;    

    exponentegaussiana=-((punto.p-eta.p)*(punto.p-eta.p)+
			(punto.q-eta.q)*(punto.q-eta.q))/hbar;
    

    faseinterferencia=punto.simplecticproduct(zeta)/hbar;


    magnitud=gsl_complex_abs(gsl_complex_mul(PesoUno,PesoDos));
    
    FactorCorrectivopormasayomega=sqrt(Uno.masa*Uno.omega*
				       Dos.masa*Dos.omega)/(pi*hbar);

   
    result=2*FactorCorrectivopormasayomega*
      magnitud*
      cos(faseinterferencia)*
      exp(exponentegaussiana);
      
    return result;

  };    


  
  gsl_complex InterferenceTermChordRepresentation(simplectic& cuerda){
    /*la funcion devuelve SOLAMENTE los terminos: 
      (son dos, sumados) de interferencia
      de la funcion de Weyl correspondiente al punto,
    la Notacion sigue el articulo de Nicacio y Raul Vallejos.
    La transformacion metaplectica se reduce a la unidad.*/

    double magnitud;
    simplectic zeta, eta;
    simplectic cuerdaefectiva;
    //la magnitud lleva un termino que no aparece en el art de Raul.
    // ya que tu consideras masas y vibraciones arbitrarias.
    double FactorCorrectivopormasayomega;
    double exponentegaussiana, faseinterferencia;


    //Recuerda que esto es un coseno: la transformada
    //se debe a dos terminos que suman con fases opuestas
    //en su version de cuerdas estos se separan!
    gsl_complex medioresultadouno, medioresultadodos;

    gsl_complex result;

    zeta=Uno.centre-Dos.centre;
    eta=(Uno.centre+Dos.centre);
    eta.q=eta.q/2.0;
    eta.p=eta.p/2.0;    


    //Primer termino
    cuerdaefectiva=cuerda-zeta;

    exponentegaussiana=-((cuerdaefectiva.p)*(cuerdaefectiva.p)+
			 (cuerdaefectiva.q)*(cuerdaefectiva.q))/(4.0*hbar);
    
    faseinterferencia=-cuerdaefectiva.simplecticproduct(eta)/hbar;
  

    medioresultadouno=gsl_complex_polar(exp(exponentegaussiana), faseinterferencia);

    
    //segundo termino
    cuerdaefectiva=cuerda+zeta;

    exponentegaussiana=-((cuerdaefectiva.p)*(cuerdaefectiva.p)+
			 (cuerdaefectiva.q)*(cuerdaefectiva.q))/(4.0*hbar);
  
    faseinterferencia=-cuerdaefectiva.simplecticproduct(eta)/hbar;
  
    medioresultadodos=gsl_complex_polar(exp(exponentegaussiana), faseinterferencia);

    magnitud=gsl_complex_abs(gsl_complex_mul(PesoUno,PesoDos));
    
    FactorCorrectivopormasayomega=sqrt(Uno.masa*Uno.omega*
				       Dos.masa*Dos.omega)/(hbar*pi);
    

    result=gsl_complex_add(medioresultadouno,medioresultadodos);

    result=gsl_complex_mul_real(result, magnitud*FactorCorrectivopormasayomega);


    return result;

  };    


  double CentreRepresentation( simplectic & punto){
    //Devuelve la representacion en Centros completa evaluada en el punto.
    double result;
    

    result=InterferenceTermCentreRepresentation(punto)
      +Uno.CentreRepresentation(punto)
      +Dos.CentreRepresentation(punto);

    return result;



    
  };

  
  gsl_complex ChordRepresentation( simplectic & punto){
    //Devuelve la representacion en Centros completa evaluada en el punto.
    gsl_complex parteaditiva;
    gsl_complex result;
    

    parteaditiva=gsl_complex_add(Uno.ChordRepresentation(punto),
			   Dos.ChordRepresentation(punto));

    result=gsl_complex_add(parteaditiva,
			   InterferenceTermChordRepresentation(punto));
	


    return result;



    
  };
  

  
};


#endif 
