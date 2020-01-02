#include<iostream>
#include<fstream>
#include<cstdlib> //exit
#include<string>
#include<stdio.h>
#include<iomanip>
#include<ctime>
#include "Grid.h"

using namespace std;

int main()
{
	double H;
	double L;
	int nH;
	int nL;
	int num_nodes;
	int num_elements;
	double ro;
	double alfa;
	double t0; 
	double t_env;
	double conductivity;
	double specific_heat;
	double sym_time;
	double step_time;

	int i = 0;
	string::size_type st;
	string line;
	int numOf_line = 1;
	fstream file;

	file.open("dane.txt", ios::in);
	if (file.good() == false) {
		cout << "File doesnt exist";
		exit(0);
	}
	while (getline(file, line)){
		switch (numOf_line){
		case 1: H = stod(line, &st); break;
		case 2: L = stod(line, &st); break;
		case 3: nH = atoi(line.c_str()); break;
		case 4: nL = atoi(line.c_str()); break;
		case 5: ro = stod(line.c_str()); break;
		case 6: alfa = stod(line.c_str()); break;
		case 7: t0 = stod(line.c_str()); break;
		case 8: t_env = stod(line.c_str()); break;
		case 9: conductivity = stod(line.c_str()); break;
		case 10: specific_heat = stod(line.c_str()); break;
		case 11: sym_time = stod(line.c_str()); break;
		case 12: step_time = stod(line.c_str()); break;
		}
		numOf_line++;
	}
	file.close();

	Grid grid = Grid(H, L, nH, nL, ro, alfa, t0, t_env, conductivity, specific_heat, sym_time, step_time);

	system("pause");

}