#pragma once
#include "Node.h"
#include "Element.h"
#include "LocalLayout.h"

class Grid {

public:
	double H; //wysokosc siatki
	double L; //szerokosc siatki
	int nH; //ilosc wezlow na wysokosc
	int nL; //ilosc wezlow na szerokosc
	int num_nodes; //ilosc wszystkich wezlow
	int num_elements; //ilosc wszystkich elementow
	double ro;
	double alfa;
	double t0;
	double t_env;
	double conductivity;
	double specific_heat;
	double sym_time;
	double step_time;
	Node *arr_nodes; // nH * nL
	Element*arr_elem; // (nH-1) * (nL-1)
	LocalLayout*locally;
	double*det_j;

	double*global_vectorP; //(nH*nL)
	double**global_matrixH; // ((nH*nL)*(nH*nL))
	double**global_matrixC;
	double *temp;
	double**create_globalH();
	double**create_globalC();
	double**finally_H();
	double*create_globalP();
	double*finally_P(double*);
	double*gauss(double **, double *);

	Grid();
	Grid(double, double, int, int, double, double, double, double, double, double, double, double);
	void create_node();
	void create_elem();
	double**jacobian(int, int);
	double**create_matrixH_dV(int);
	double**create_matrixC(int);
	double** create_matrixH_dS(int);
	double**matrixC_dev_dtau();
	double** matrixH_full(int);
	double* create_vectorP(int);
	double detJ(int, int);
	double*size_length();

	
};