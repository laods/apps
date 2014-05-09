all: helloMCI upscale_steadystate_analytic printForceInfo mirror_cpgrid

helloMCI:
	g++ -O2 -o bin/helloMCI source/helloMCI.cpp -lopmcore

upscale_steadystate_analytic:
	g++ -O2 -o bin/upscale_steadystate_analytic source/upscale_steadystate_analytic.cpp -lopmcore

printForceInfo:
	g++ -O2 -std=gnu++0x -o bin/printForceInfo source/printForceInfo.cpp -lopmcore -lopmporsol -ldunecornerpoint -llapack -lblas

mirror_cpgrid:
	g++ -O2 -std=gnu++0x -o bin/mirror_cpgrid source/mirror_cpgrid.cpp -lopmcore

correlationlengths:
	g++ -O2 -std=gnu++0x -o bin/correlationlengths source/correlationlengths.cpp -lopmcore -lopmporsol -ldunecornerpoint -llapack -lblas
