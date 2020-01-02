#include "Element.h"
#include <stdio.h>
#include <ctime>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

Element::Element() {
	this->id_element = 0;
	for (int i = 0; i < 4; i++) { //bo 4 elementy w id elementu
		this->id_node[i] = 0;
		this->id_walls[i] = 0;
		this->localvectorP[i] = 0;
	}

	this->localmatrixH = new double*[4];
	for (int i = 0; i < 4; i++) {
		this->localmatrixH[i] = new double[4];
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			this->localmatrixH[i][j] = 0;
		}
	}
}

Element::Element(int id, int*id_nodes, int*id_walls, double**localmatrixH, double*localvectorP) {
	this->id_element = id;
	for (int i = 0; i < 4; i++)
	{
		this->id_node[i] = id_nodes[i];
		this->id_walls[i] = id_walls[i];
		this->localvectorP[i] = localvectorP[i];
	}
	this->localmatrixH = new double*[4];
	for (int i = 0; i < 4; i++) {
		this->localmatrixH[i] = new double[4];
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			this->localmatrixH[i][j] = localmatrixH[i][j];	
		}
	}
}

Element::~Element(){}


