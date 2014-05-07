DOCUMENTACION PARA EL PAQUETE DE GATOS
======================================

Primera Parte: Estados Coherentes.
----------------------------------------

La primera parte del paquete consiste
en un archivo de objetos llamado CoherentState*.hpp
Es necesario incluir este archivo con el 
indicativo 

	#include "CoherentState*.hpp"

El paquete depende de `gsl`, por lo que
es necesario pasar la opción de compilación `-lgsl`
y de preferencia optimizar. También usa
el paquete `simplectic01.hpp` que debe encontrarse
en la ruta de compilación.


El asterisco representa la versión.

Para el uso de estados coherentes en semiclásica
usé las convenciones Ozorio de Almeida, Raul Vallejos
y otras. Un estado coherente es simplemente
un estado base del oscilador armónico fuera del
origen. Un estado coherente también se le conoce
como estado Gaussiano por su representación
de Wigner (centros). 

Para simplificar el asunto hablare de un 
estado coherente usando la siguiente notación:
El estado base de un oscilador armónico normalizado
será denotado por

|0>

El operador de Traslación por una cuerda \mu
aplicado a |0> resultara en el estado coherente
centrado en (mu_q, mu_p):

T_\mu |0> = |\mu>

De momento el paquete que he preparado solamente puede
hacer representaciones estáticas. He escogido
las más usuales para empezar. El estado asume
que tenemos en general una masa _m_ y una 
frecuencia oscilatoria _w_

#### Constructores

El constructor por default `CoherentState()`
produce un estado coherente centrado en el origen.
Por omisión, _m_=1, _w_=1.

El constructor con un argumento `simplectic` da
un estado coherente centrado en el argumento,
es decir, `CoherentState(x)` con x de tipo `simplectic`
da un estado coherente |x>.
La variante de este constructor es usar
dos argumentos tipo double, `CoherentState(q,p)`,
con el resultado esperado.
Finalmente podemos llamar al constructor
con otro estado coherente como argumento,
resultando en uno cuyas cualidades
son las  del original:
 
	CohentState OtroEstado;
	CoherentState(OtroEstado) X;


####Metodos para accesar las propiedades.

Por filosofia, todas las propiedades
son publicas, asi que se pueden acceder
directamente del objeto con un punto y el nombre
de la propiedad. Las propiedades son:
	simplectic centre;
	double masa;
	double omega;

Además hay otros metodos para acceder a las
propiedades, cuyo nombre revela su funcion.

	SetCentre(simplectic);
	SetCentre(double, double);
	SetMasa(double);
	SetOmega(double);

Y tenemos un operador de `=` (asignación). 





#### Estado Coherente: Representación en _q_

La representación en el espacio de posiciones
o de las _q_ es una funcion de amplitud compleja.
Dado que mi código todavía no cuenta con
propagadores, siempre asumimos el tiempo igual a cero.
Con ello, la representacion