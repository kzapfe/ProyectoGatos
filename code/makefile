LIBS=-lgsl
FUENTECOM=CoherentState01.hpp simplectic01.hpp QuantumConstants.hpp 
FUENTESESFERA=CatStates01.hpp PopulateSphericalShells01.hpp

TestCentre01.x: testeando01.cpp $(FUENTECOM)
	g++ $(LIBS) $<  -o $@


TestChords01.x: testeando02.cpp $(FUENTECOM)
	g++ $(LIBS) $<  -o $@


TestCats01.x: testeandogatos01.cpp $(FUENTECOM) CatStates01.hpp 
	g++ $(LIBS) $<  -o $@


TestCats02.x: testeandogatos02.cpp $(FUENTECOM) CatStates01.hpp 
	g++ $(LIBS) $<  -o $@

TestCircle01.x: testcircle01.cpp $(FUENTECOM)  $(FUENTESESFERA)
	g++ $(LIBS) -larmadillo  $<  -o $@

TestSphere01.x: testeandogatosesfera01.cpp $(FUENTECOM)  $(FUENTESESFERA)
	g++ $(LIBS) -larmadillo  $<  -o $@


TestSphereParalel01.x: testcatssphereparalel01.cpp $(FUENTECOM)  $(FUENTESESFERA)
	g++ $(LIBS) -fopenmp -larmadillo  $<  -o $@

#Not really testing anymore
SphereParalel01.x: catspheresectionparalel01.cpp CortesGatosBidimensionales01.hpp \
$(FUENTECOM)  $(FUENTESESFERA)
	g++ $(LIBS) -fopenmp -larmadillo  $<  -o $@	

CatGenParalel01.x: catsgeneralsectionparalel01.cpp CortesGatosBidimensionales01.hpp \
$(FUENTECOM)  $(FUENTESESFERA)
	g++ $(LIBS) -fopenmp -larmadillo  $<  -o $@	

CatGenParalel02.x: catsgeneralsectionparalel02.cpp CortesGatosBidimensionales02.hpp \
$(FUENTECOM)  $(FUENTESESFERA)
	g++ $(LIBS) -fopenmp -larmadillo  $<  -o $@	
