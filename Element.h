#pragma once
class Element {

public:
	int id_element;
	int id_node[4];
	int id_walls[4];
	double localvectorP[4];
	double**localmatrixH;


	Element();
	Element(int, int*, int*, double**, double*);
	~Element();
};