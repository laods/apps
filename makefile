all: helloMCI upscale_steadystate_analytic

helloMCI:
	g++ -O2 -o helloMCI helloMCI.cpp -lopmcore

upscale_steadystate_analytic:
	g++ -O2 -o upscale_steadystate_analytic upscale_steadystate_analytic.cpp -lopmcore