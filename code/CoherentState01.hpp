/*Rutinas y clases y metodos para hacer estados coherentes y superposiciones */
/* Karel Zapfe */

#ifndef __CoherentState__
#define __CoherentState__


#include "QuantumConstants.hpp"
#include "simplectic01.hpp"
#include <algorithm>
#include <cmath>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>


//Piensa donde poner hbarra, depende del estado GLOBAL.



class CoherentState{

public:


  simplectic centre;
  
  //En principio, la definicion del estado coherente depende de
  //la aproximacion cuadratica de la Func Hamiltoniana. 
  //asi que aparecen los siguientes valores.

  double masa, omega;

   CoherentState(){
    this->centre.q=0;
    this->centre.p=0;
    this->masa=1;
    this->omega=1;
  };

  
   CoherentState(simplectic & desplazamiento){
    this->centre.q=desplazamiento.q;
    this->centre.p=desplazamiento.p;
    this->masa=1;
    this->omega=1;
  };


  
  CoherentState(double q, double p){
    this->centre.q=q;
    this->centre.p=p;
    this->masa=1;
    this->omega=1;
  };


  CoherentState(CoherentState & Original){
    //Constructor de Copias.
    this->centre.q=Original.centre.q;
    this->centre.p=Original.centre.p;
    this->masa=Original.masa;
    this->omega=Original.omega;
  };


void  SetCentre(simplectic &x){
    centre.q=x.q;
    centre.p=x.p;

  };

  
void  SetCentre(double q, double p){
    centre.q=q;
    centre.p=p;

  };


  friend void swap(CoherentState& Uno, CoherentState& Dos){
    //intercambia dos CoherentStates, depende de <algorithm>
    using std::swap;
    swap(Uno.centre, Dos.centre);    
    swap(Uno.masa, Dos.masa);    
    swap(Uno.omega, Dos.omega);    
  };

  CoherentState& operator=(CoherentState RightHandSide){
    swap(*this, RightHandSide);
    return *this;
  }


  
   gsl_complex q_RepresentationHeisen(double xq){
     /*Representacion q en Heisenberg pictures. La parte temporal esta
       en la observable. */
     double Normalizador;
     double exponentereal;
     double exponenteimag;
     double Magnitud;
     gsl_complex result;

     Normalizador=pow( masa*omega/pi/hbar, 0.25);
     exponentereal=- omega*masa/(2*hbar)*(xq-centre.q)*(xq-centre.q);
     exponenteimag=centre.p*xq/hbar;
     
     Magnitud=Normalizador*exp(exponentereal);  
     result=gsl_complex_polar(Magnitud, exponenteimag);
     
     return result;
     
   };
   

   gsl_complex p_RepresentationHeisen(double xp){
     /*Representacion p en Heisenberg pictures. La parte temporal esta
       en la observable. Tu obtuviste esto, asi que esperemos este vien*/
     double Normalizador;
     double exponentereal;
     double exponenteimag;
     double Magnitud;
     gsl_complex result;

     Normalizador=pow( masa*omega/pi/hbar, 0.25);
     exponentereal=-(xp-centre.p)*(xp-centre.p)/(masa*omega*2.0*hbar);
     exponenteimag=-centre.q*xp/hbar;

     Magnitud=Normalizador*exp(exponentereal);  
     result=gsl_complex_polar(Magnitud, exponenteimag);
     
     return result;
     
   };


  double CentreRepresentation(simplectic & point){
    //Funcion de Wigner a t=0 (a la Heisenberg maybe) 
    //de un estado coherente.
    //no metaplectic transforms here please.
    /* Aqui point es el punto de evaluacion, un elemento del espacio
       de centros. */
    
    double Normalizador;
    double exponente;
    double result;

    Normalizador=masa*omega/(pi*hbar);
    exponente=-((centre.q-point.q)*(centre.q-point.q)*omega*masa
		+(centre.p-point.p)*(centre.p-point.p)/omega*masa)/hbar;
    result=Normalizador*exp(exponente);
    
    return result;

  };

  double CentreRepresentation(double q, double p){
    simplectic puntoaux(q, p);
    double result=CentreRepresentation(puntoaux);
    return result;
  };


  gsl_complex ChordRepresentation(simplectic & point){
    //Funcion de Wigner a t=0 (a la Heisenberg maybe) 
    //de un estado coherente.
    //no metaplectic transforms here please.
    /* Aqui point es el punto de evaluacion, un elemento del espacio
       de cuerdas. */


    double Normalizador;
    double exponentereal;
    double exponenteimag;
    gsl_complex exponente;
    gsl_complex result;

    Normalizador=masa*omega/(pi*hbar);
    exponentereal=(-1.0/(4.0*hbar))*
      ((point.p*point.p)/(masa*omega)+(point.q*point.q)*(masa*omega));
    exponenteimag=-1/hbar*(centre.simplecticproduct(point));

    exponente=gsl_complex_rect(exponentereal, exponenteimag);
    result=gsl_complex_mul_real(gsl_complex_exp(exponente), Normalizador);

    return result;

  };

  gsl_complex ChordRepresentation(double q, double p){
    simplectic puntoaux(q, p);
    gsl_complex result= ChordRepresentation(puntoaux);
    return result;
  };

  
};



#endif
