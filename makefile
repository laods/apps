all: helloMCI upscale_steadystate_analytic printForceInfo mirror_cpgrid

helloMCI:
	g++ -O2 -o bin/helloMCI source/helloMCI.cpp -lopmcore

upscale_steadystate_analytic:
	g++ -O2 -o bin/upscale_steadystate_analytic source/upscale_steadystate_analytic.cpp -lopmcore

printForceInfo:
	g++ -O2 -std=gnu++0x -o bin/printForceInfo -I/private/laods/opm/RH6/april15/install/include -L/private/laods/opm/RH6/april15/install/lib64 -I/project/res/x86_64_RH_5/share/opm/dune-2.3/dune-geometry/include -I/project/res/x86_64_RH_5/share/opm/dune/dune-common/include -I/project/res/x86_64_RH_5/share/opm/dune/dune-grid/include -L/project/res/x86_64_RH_5/share/opm/dune/dune-grid/lib -I/project/res/x86_64_RH_5/share/ert/release/nightly/upstream/master/include source/printForceInfo.cpp -lopmcore -lopmporsol -ldunecornerpoint -llapack -lblas

mirror_cpgrid:
	g++ -O2 -std=gnu++0x -o bin/mirror_cpgrid source/mirror_cpgrid.cpp -lopmcore

correlationlengths:
	g++ -O2 -std=gnu++0x -o bin/correlationlengths source/correlationlengths.cpp -lopmcore -lopmporsol -ldunecornerpoint -llapack -lblas

debug:
	g++ -g3 -O0 -std=gnu++0x -o bin/correlationlengths source/correlationlengths.cpp -lopmcore -lopmporsol -ldunecornerpoint -llapack -lblas
