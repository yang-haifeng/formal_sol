all : fml_sol

fml_sol : *.cpp *.h
	g++ *.cpp -o fml_sol
