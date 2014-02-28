all: helloMCI upscale_steadystate_analytic printForceInfo mirror_cpgrid

helloMCI:
	g++ -O2 -o helloMCI helloMCI.cpp -lopmcore

upscale_steadystate_analytic:
	g++ -O2 -o upscale_steadystate_analytic upscale_steadystate_analytic.cpp -lopmcore

printForceInfo:
	g++ -O2 -std=gnu++0x -o printForceInfo printForceInfo.cpp -lopmcore -lopmporsol -ldunecornerpoint -llapack -lblas

mirror_cpgrid:
	g++ -O2 -std=gnu++0x -o mirror_cpgrid mirror_cpgrid.cpp -lopmcore
