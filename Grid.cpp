#include "Grid.h"
#include <stdio.h>
#include <ctime>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
using namespace std;

Grid::Grid() {

	this->H = 0.0;
	this->L = 0.0;
	this->nH = 0;
	this->nL = 0;
	this->num_nodes = 0;
	this->num_elements = 0;
	this->ro = 0.0;
	this->alfa = 0.0;
	this->t0 = 0.0;
	this->t_env = 0.0;
	this->conductivity = 0.0;
	this->specific_heat = 0.0;
	this->sym_time = 0.0;
	this->step_time = 0.0;
	this->arr_elem = new Element[this->num_elements];
	this->arr_nodes = new Node[this->num_nodes];
	this->locally = new LocalLayout();
	this->det_j = new double[4];
	create_elem();
	create_node();
	
	for (int i = 0; i < 4; i++) {
		this->det_j[i] = 0;
	}
	this->global_matrixH = new double*[(nH*nL)];
	for (int i = 0; i < (nH*nL); i++) {
		this->global_matrixH[i] = new double[(nH*nL)];
	}
	for (int i = 0; i < (nH*nL); i++) {
		for (int j = 0; j < (nH*nL); j++) {
			this->global_matrixH[i][j] = 0;
		}
	}
	this->global_matrixC = new double*[(nH*nL)];
	for (int i = 0; i < (nH*nL); i++) {
		this->global_matrixC[i] = new double[(nH*nL)];
	}
	for (int i = 0; i < (nH*nL); i++) {
		for (int j = 0; j < (nH*nL); j++) {
			this->global_matrixC[i][j] = 0;
		}
	}
	this->global_vectorP = new double[(nH*nL)];
	for (int i = 0; i < (nH*nL); i++) {
		this->global_vectorP[i] = 0;
	}
	this->temp = new double[(nH*nL)];
	for (int i = 0; i < (nH*nL); i++) {
		this->temp[i] = 0;
	}
}///////////////////////////////////////////////////////////////// MAIN //////////////////////////////////////////////////////////////////////////////////////
Grid::Grid(double H, double L, int nH, int nL, double ro, double alfa, double t0, double t_env, double conductivity, double specific_heat, double sym_time, double step_time) {
	this->H = H;
	this->L = L;
	this->nH = nH;
	this->nL = nL;
	this->ro = ro;
	this->alfa = alfa;
	this->t0 = t0;
	this->t_env = t_env;
	this->conductivity = conductivity;
	this->specific_heat = specific_heat;
	this->sym_time = sym_time;
	this->step_time = step_time;
	this->num_nodes = nH*nL;
	this->num_elements = (nH-1)*(nL-1);
	this->arr_elem = new Element[num_elements];
	this->arr_nodes = new Node[num_nodes];
	this->locally = new LocalLayout();
	this->det_j = new double[4];
	create_elem();
	create_node();

	this->global_matrixH = new double*[(nH*nL)];
	for (int i = 0; i < (nH*nL); i++) {
		this->global_matrixH[i] = new double[(nH*nL)];
	}
	for (int i = 0; i < (nH*nL); i++) {
		for (int j = 0; j < (nH*nL); j++) {
			this->global_matrixH[i][j] = 0;
		}
	}
	this->global_matrixC = new double*[(nH*nL)];
	for (int i = 0; i < (nH*nL); i++) {
		this->global_matrixC[i] = new double[(nH*nL)];
	}
	for (int i = 0; i < (nH*nL); i++) {
		for (int j = 0; j < (nH*nL); j++) {
			this->global_matrixC[i][j] = 0;
		}
	}
	this->global_vectorP = new double[(nH*nL)];
	for (int i = 0; i < (nH*nL); i++) {
		this->global_vectorP[i] = 0;
	}/////////////////////////////
	double*T1 = new double[(nH*nL)];
	for (int i = 0; i < (nH*nL); i++) {
		T1[i] = t0;
	}////////////////////////////////
	double**Hglob = new double*[(nH*nL)];
	for (int i = 0; i < (nH*nL); i++) {
		Hglob[i] = new double[(nH*nL)];
	}
	for (int i = 0; i < (nH*nL); i++) {
		for (int j = 0; j < (nH*nL); j++) {
			Hglob[i][j] = 0;
		}
	}
	double*Pglob = new double[(nH*nL)];
	for (int i = 0; i < (nH*nL); i++) {
		Pglob[i] = 0;
	}
	double min, max;
	for (int t = step_time; t <= sym_time; t += step_time) {
		cout << "Time[s] : " << t << endl;
		//cout << endl<<"Matrix H:" << endl;
		Hglob=finally_H();

	/*	for (int i = 0; i < (nH*nL); i++) {
			for (int j = 0; j < (nH*nL); j++) {
				cout<<Hglob[i][j] << " ";
			}
			cout << endl;
		}*/
		//cout <<endl<< "Vector P:" << endl;
		Pglob = finally_P(T1);
		/*for (int i = 0; i < (nH*nL); i++){
			cout << Pglob[i] << " ";
		}*/
		double **HP = new double *[num_nodes];
		for (int j = 0; j < num_nodes; j++) {
			HP[j] = new double[num_nodes + 1];
		}
		for (int i = 0; i < num_nodes; i++) {
			for (int j = 0; j < num_nodes; j++) {
				HP[i][j] = Hglob[i][j];
			}
		}
		for (int i = 0; i < num_nodes; i++) {
			HP[i][(nH*nL)] = Pglob[i];
		}
		T1 = gauss(HP, T1);
		/*cout << "Temperatures: " << endl;
		for (int i = 0; i < (nH*nL); i++) {
			cout << T1[i] << " ";
		}*/
		for (int i = 0; i < (nH*nL); i++) {
			if (i == 0) {
				min = T1[i];
				max = T1[i];
			}
			else {
				if (T1[i] < min)	min = T1[i];
				if (T1[i] > max)	max = T1[i];
			}
		}
		cout<< "MIN: " << min << ", MAX: " << max << endl<<endl;
		min = 0;
		max = 0;
		for (int i = 0; i < (nH*nL); i++) {
			Pglob[i] = 0;
		}
	}
	
}//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double *Grid::size_length() {
	double*dimensions = new double[4];
	for (int i = 0; i < 4; i++) {
		if (i % 2 == 1) dimensions[i] = (this->H) / (this->nH - 1);
		else dimensions[i] = (this->L) / (this->nL - 1);
	}
	return dimensions;
}
void Grid::create_elem() {

	int iter = 1;

	for (int i = 0; i < this->num_elements; i++) {
		arr_elem[i].id_element = i + 1; //elementy numeruje od 1 

		arr_elem[i].id_node[0] = iter;
		arr_elem[i].id_node[1] = iter + nH;
		arr_elem[i].id_node[2] = iter + nH + 1;
		arr_elem[i].id_node[3] = iter + 1;
		iter++;
		if (iter % nH == 0) iter++;

		for (int j = 0; j < 4; j++) {
			arr_elem[i].id_walls[j] = j + 1;
		}
	}
}
void Grid::create_node() {
	int i = 0;
	for (int a = 0; a < this->nL; a++) {
		for (int b = 0; b < this->nH; b++) {
			arr_nodes[i].id = i + 1;
			arr_nodes[i].x = a * (L / (nL - 1));
			arr_nodes[i].y = b * (H / (nH - 1));
			if (a != 0 && a != (nL - 1) && b != 0 && b != (nH - 1)) {
				arr_nodes[i].edge = false;
			}
			else {
				arr_nodes[i].edge = true;
			}
			if (i < num_nodes) { 
				i++;
			}
		}
	}
}//w ukladzie lokalnym wspolrzedne punktow to 0,1,2,3, wiec podaje jeden z tych punktow
double** Grid::jacobian(int id_elem, int point) { //tak naprawde liczy macierz H po dV xD
	double**jacobian = new double*[2];
	for (int i = 0; i < 2; i++) {
		jacobian[i] = new double[2];
	}
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			jacobian[i][j] = 0;
		}
	}
	for (int i = 0; i < 4; i++) {
		det_j[i] = 0;
	}//[point][0] - 0 bo N1 tylko we wzorze, 1 bo N2 itd
	jacobian[0][0] = ((this->locally->dNdksi[point][0])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[0] - 1].x + (this->locally->dNdksi[point][1])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[1] - 1].x + (this->locally->dNdksi[point][2])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[2] - 1].x + (this->locally->dNdksi[point][3])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[3] - 1].x);
	jacobian[0][1] = ((this->locally->dNdksi[point][0])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[0] - 1].y + (this->locally->dNdksi[point][1])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[1] - 1].y + (this->locally->dNdksi[point][2])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[2] - 1].y + (this->locally->dNdksi[point][3])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[3] - 1].y);
	jacobian[1][0] = ((this->locally->dNdeta[point][0])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[0] - 1].x + (this->locally->dNdeta[point][1])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[1] - 1].x + (this->locally->dNdeta[point][2])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[2] - 1].x + (this->locally->dNdeta[point][3])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[3] - 1].x);
	jacobian[1][1] = ((this->locally->dNdeta[point][0])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[0] - 1].y + (this->locally->dNdeta[point][1])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[1] - 1].y + (this->locally->dNdeta[point][2])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[2] - 1].y + (this->locally->dNdeta[point][3])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[3] - 1].y);
	det_j[point] = (jacobian[0][0] * jacobian[1][1]) - (jacobian[0][1] * jacobian[1][0]);
	double**reverse_j = new double*[2];
	for (int i = 0; i < 2; i++) {
		reverse_j[i] = new double[2];
	}//odwracanie macierzy
	reverse_j[0][0] = (1 / det_j[point])*jacobian[1][1];
	reverse_j[0][1] = (1 / det_j[point])*(-jacobian[0][1]);
	reverse_j[1][0] = (1 / det_j[point])*(-jacobian[1][0]);
	reverse_j[1][1] = (1 / det_j[point])*jacobian[0][0];
	double**dNdx = new double*[2];
	for (int i = 0; i < 4; i++) {
		dNdx[i] = new double[2];
	}
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			dNdx[i][j] = 0;
		}
	}
	double**dNdy = new double*[2];
	for (int i = 0; i < 4; i++) {
		dNdy[i] = new double[2];
	}
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			dNdy[i][j] = 0;
		}
	}
	dNdx[0][0] = reverse_j[0][0] * (this->locally->dNdksi[point][0]) + reverse_j[0][1] * (this->locally->dNdeta[point][0]);
	dNdx[0][1] = reverse_j[0][0] * (this->locally->dNdksi[point][1]) + reverse_j[0][1] * (this->locally->dNdeta[point][1]);
	dNdx[1][0] = reverse_j[0][0] * (this->locally->dNdksi[point][2]) + reverse_j[0][1] * (this->locally->dNdeta[point][2]);
	dNdx[1][1] = reverse_j[0][0] * (this->locally->dNdksi[point][3]) + reverse_j[0][1] * (this->locally->dNdeta[point][3]);

	dNdy[0][0] = reverse_j[1][0] * (this->locally->dNdksi[point][0]) + reverse_j[1][1] * (this->locally->dNdeta[point][0]);
	dNdy[0][1] = reverse_j[1][0] * (this->locally->dNdksi[point][1]) + reverse_j[1][1] * (this->locally->dNdeta[point][1]);
	dNdy[1][0] = reverse_j[1][0] * (this->locally->dNdksi[point][2]) + reverse_j[1][1] * (this->locally->dNdeta[point][2]);
	dNdy[1][1] = reverse_j[1][0] * (this->locally->dNdksi[point][3]) + reverse_j[1][1] * (this->locally->dNdeta[point][3]);
	double **dNdxdNdxT = new double*[4];
	for (int i = 0; i < 4; i++){
		dNdxdNdxT[i] = new double[4];
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			dNdxdNdxT[i][j] = 0;
		}
	}
	double **dNdydNdyT = new double*[4];
	for (int i = 0; i < 4; i++){
		dNdydNdyT[i] = new double[4];
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			dNdydNdyT[i][j] = 0;
		}
	}//////////////////////// dNdx*{dNdx}T i dNdy*{dNdy}T
	double tab_dNdx[4];
	for (int i = 0; i < 4; i++) {
		tab_dNdx[i] = 0;
	}
	tab_dNdx[0] = dNdx[0][0];
	tab_dNdx[1] = dNdx[0][1];
	tab_dNdx[2] = dNdx[1][0];
	tab_dNdx[3] = dNdx[1][1];
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			dNdxdNdxT[i][j] = tab_dNdx[i] * tab_dNdx[j];	
		}
	}
	double tab_dNdy[4];
	for (int i = 0; i < 4; i++) {
		tab_dNdy[i] = 0;
	}
	tab_dNdy[0] = dNdy[0][0];
	tab_dNdy[1] = dNdy[0][1];
	tab_dNdy[2] = dNdy[1][0];
	tab_dNdy[3] = dNdy[1][1];
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			dNdydNdyT[i][j] = tab_dNdy[i] * tab_dNdy[j];
		}
	}////////// mnozenie razy wyznacznik //////////////////////
	double **dNdxdNdxT_det_j = new double*[4];
	for (int i = 0; i < 4; i++){
		dNdxdNdxT_det_j[i] = new double[4];
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			dNdxdNdxT_det_j[i][j] = 0;
		}
	}
	double **dNdydNdyT_det_j = new double*[4];
	for (int i = 0; i < 4; i++){
		dNdydNdyT_det_j[i] = new double[4];
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			dNdydNdyT_det_j[i][j] = 0;
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			dNdxdNdxT_det_j[i][j] = dNdxdNdxT[i][j] * det_j[point];
			dNdydNdyT_det_j[i][j] = dNdydNdyT[i][j] * det_j[point];
		}
	}
	double ** sum_multiply_k = new double*[4];
	for (int i = 0; i < 4; i++) {
		sum_multiply_k[i] = new double[4];
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			sum_multiply_k[i][j] = 0;
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			sum_multiply_k[i][j] = (dNdxdNdxT_det_j[i][j] + dNdydNdyT_det_j[i][j])*conductivity;
		}
	}
	return sum_multiply_k;
}
double** Grid::create_matrixH_dV(int id_elem) {
	double**jac1 = new double*[4];
	for (int i = 0; i < 4; i++) {
		jac1[i] = new double[4];
	}
	double**jac2 = new double*[4];
	for (int i = 0; i < 4; i++) {
		jac2[i] = new double[4];
	}
	double**jac3 = new double*[4];
	for (int i = 0; i < 4; i++) {
		jac3[i] = new double[4];
	}
	double**jac4 = new double*[4];
	for (int i = 0; i < 4; i++) {
		jac4[i] = new double[4];
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			jac1[i][j] = 0;
			jac2[i][j] = 0;
			jac3[i][j] = 0;
			jac4[i][j] = 0;
		}
	}
	jac1 = jacobian(id_elem, 0);
	jac2 = jacobian(id_elem, 1);
	jac3 = jacobian(id_elem, 2);
	jac4 = jacobian(id_elem, 3);
	double**matrixH = new double*[4];
	for (int i = 0; i < 4; i++) {
		matrixH[i] = new double[4];
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			matrixH[i][j] = 0;
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			matrixH[i][j] = jac1[i][j] + jac2[i][j] + jac3[i][j] + jac4[i][j];
		}
	}
	return matrixH;
}
double**Grid::create_matrixC(int id_elem) ///// [C] = c*ro*NNT*detJ /////////////////////
{
	double det0= (detJ(id_elem, 0));
	double det1 = (detJ(id_elem, 1));
	double det2 = (detJ(id_elem, 2));
	double det3 = (detJ(id_elem, 3));
	double **tab_NN1=new double*[4];
	double **tab_NN2 = new double*[4];
	double **tab_NN3 = new double*[4];
	double **tab_NN4 = new double*[4];
	for (int i = 0; i < 4; i++) {
		tab_NN1[i] = new double[4];
		tab_NN2[i] = new double[4];
		tab_NN3[i] = new double[4];
		tab_NN4[i] = new double[4];
	}
	for (int i = 0; i < 4; i++){
		for (int j = 0; j < 4; j++){
			tab_NN1[i][j] = (this->locally->NdV[0][i])*(this->locally->NdV[0][j])*specific_heat;
			tab_NN1[i][j] *= ro;
			tab_NN1[i][j] *= det0;
		}
		for (int j = 0; j < 4; j++){
			tab_NN2[i][j] = (this->locally->NdV[1][i])*(this->locally->NdV[1][j])*specific_heat;
			tab_NN2[i][j] *= ro;
			tab_NN2[i][j] *= det1;
		}
		for (int j = 0; j < 4; j++){
			tab_NN3[i][j] = (this->locally->NdV[2][i])*(this->locally->NdV[2][j])*specific_heat;
			tab_NN3[i][j] *= ro;
			tab_NN3[i][j] *= det2;
		}
		for (int j = 0; j < 4; j++){
			tab_NN4[i][j] = (this->locally->NdV[3][i])*(this->locally->NdV[3][j])*specific_heat;
			tab_NN4[i][j] *= ro;
			tab_NN4[i][j] *= det3;
		}
	}
	double**matrixC = new double*[4];
	for (int i = 0; i < 4; i++) {
		matrixC[i] = new double[4];
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			matrixC[i][j] = 0;
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			matrixC[i][j] = tab_NN1[i][j] + tab_NN2[i][j] + tab_NN3[i][j] + tab_NN4[i][j];
		}
	}
	return matrixC;
}
//// zwraca tylko wyznacznik (potrzebny np do usuniecia calki dV w macierzy C)
double Grid::detJ(int id_elem, int point) {
	double**jacobian = new double*[2];
	for (int i = 0; i < 2; i++) {
		jacobian[i] = new double[2];
	}
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			jacobian[i][j] = 0;
		}
	}
	for (int i = 0; i < 4; i++) {
		det_j[i] = 0;
	}
	jacobian[0][0] = ((this->locally->dNdksi[point][0])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[0] - 1].x + (this->locally->dNdksi[point][1])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[1] - 1].x + (this->locally->dNdksi[point][2])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[2] - 1].x + (this->locally->dNdksi[point][3])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[3] - 1].x);
	jacobian[0][1] = ((this->locally->dNdksi[point][0])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[0] - 1].y + (this->locally->dNdksi[point][1])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[1] - 1].y + (this->locally->dNdksi[point][2])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[2] - 1].y + (this->locally->dNdksi[point][3])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[3] - 1].y);
	jacobian[1][0] = ((this->locally->dNdeta[point][0])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[0] - 1].x + (this->locally->dNdeta[point][1])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[1] - 1].x + (this->locally->dNdeta[point][2])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[2] - 1].x + (this->locally->dNdeta[point][3])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[3] - 1].x);
	jacobian[1][1] = ((this->locally->dNdeta[point][0])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[0] - 1].y + (this->locally->dNdeta[point][1])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[1] - 1].y + (this->locally->dNdeta[point][2])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[2] - 1].y + (this->locally->dNdeta[point][3])*this->arr_nodes[this->arr_elem[id_elem - 1].id_node[3] - 1].y);
	det_j[point] = (jacobian[0][0] * jacobian[1][1]) - (jacobian[0][1] * jacobian[1][0]);
	return (det_j[point]);
}
double**Grid::create_matrixH_dS(int id_elem) { // [H]= alfa * NNT * wyznacznik (zeby znikla calka po dS)
	double** matrixH = new double*[4];
	double** NNT1 = new double*[4]; 
	double** NNT2 = new double*[4];
	double** NNT3 = new double*[4];
	double** NNT4 = new double*[4];
	double** NNT5 = new double*[4];
	double** NNT6 = new double*[4];
	double** NNT7 = new double*[4];
	double** NNT8 = new double*[4];
	double**suma1 = new double*[4];
	double**suma2 = new double*[4];
	double**suma3 = new double*[4];
	double**suma4 = new double*[4];
	for (int i = 0; i < 4; i++) {
		matrixH[i] = new double[4];
		NNT1[i] = new double[4];
		NNT2[i] = new double[4];
		NNT3[i] = new double[4];
		NNT4[i] = new double[4];
		NNT5[i] = new double[4];
		NNT6[i] = new double[4];
		NNT7[i] = new double[4];
		NNT8[i] = new double[4];
		suma1[i] = new double[4];
		suma2[i] = new double[4];
		suma3[i] = new double[4];
		suma4[i] = new double[4];
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			matrixH[i][j] = 0;
			NNT1[i][j] = 0;
			NNT2[i][j] = 0;
			NNT3[i][j] = 0;
			NNT4[i][j] = 0;
			NNT5[i][j] = 0;
			NNT6[i][j] = 0;
			NNT7[i][j] = 0;
			NNT8[i][j] = 0;
			suma1[i][j] = 0;
			suma2[i][j] = 0;
			suma3[i][j] = 0;
			suma4[i][j] = 0;
		}
	}
	double** N = new double*[8]; //funkcje ksztaltu
		for (int i = 0; i < 8; i++) {
			N[i] = new double[4];
		}
		for (int i = 0; i < 8; i++) {
			for (int j = 0; j < 4; j++) {
				N[i][j] = 0;
			}
		}
		double detJ = 0; //wyznacznik (bok/2)
	if (this->arr_nodes[this->arr_elem[id_elem - 1].id_node[0] - 1].edge == 1 && this->arr_nodes[this->arr_elem[id_elem - 1].id_node[1] - 1].edge == 1) { // cala krawedz jest boczna
		//warunki brzegowe sa
		//na kazdej scianie 2 punkty i dla kazdego ksi i eta 2x2x2 i tutaj liczymy funkcje ksztaltu z LocalLayout
		N[0][0] = this->locally->NdS[0][0];
		N[0][1] = this->locally->NdS[0][1];
		N[0][2] = this->locally->NdS[0][2];
		N[0][3] = this->locally->NdS[0][3];
		N[1][0] = this->locally->NdS[1][0];
		N[1][1] = this->locally->NdS[1][1];
		N[1][2] = this->locally->NdS[1][2];
		N[1][3] = this->locally->NdS[1][3];
	}
	if (this->arr_nodes[this->arr_elem[id_elem - 1].id_node[1] - 1].edge == 1 && this->arr_nodes[this->arr_elem[id_elem - 1].id_node[2] - 1].edge == 1) {	
		N[2][0] = this->locally->NdS[2][0];
		N[2][1] = this->locally->NdS[2][1];
		N[2][2] = this->locally->NdS[2][2];
		N[2][3] = this->locally->NdS[2][3];
		N[3][0] = this->locally->NdS[3][0];
		N[3][1] = this->locally->NdS[3][1];
		N[3][2] = this->locally->NdS[3][2];
		N[3][3] = this->locally->NdS[3][3];
	}
	if (this->arr_nodes[this->arr_elem[id_elem - 1].id_node[2] - 1].edge == 1 && this->arr_nodes[this->arr_elem[id_elem - 1].id_node[3] - 1].edge == 1) {
		N[4][0] = this->locally->NdS[4][0];
		N[4][1] = this->locally->NdS[4][1];
		N[4][2] = this->locally->NdS[4][2];
		N[4][3] = this->locally->NdS[4][3];
		N[5][0] = this->locally->NdS[5][0];
		N[5][1] = this->locally->NdS[5][1];
		N[5][2] = this->locally->NdS[5][2];
		N[5][3] = this->locally->NdS[5][3];
	}
	if (this->arr_nodes[this->arr_elem[id_elem - 1].id_node[3] - 1].edge == 1 && this->arr_nodes[this->arr_elem[id_elem - 1].id_node[0] - 1].edge == 1) {
		N[6][0] = this->locally->NdS[6][0];
		N[6][1] = this->locally->NdS[6][1];
		N[6][2] = this->locally->NdS[6][2];
		N[6][3] = this->locally->NdS[6][3];
		N[7][0] = this->locally->NdS[7][0];
		N[7][1] = this->locally->NdS[7][1];
		N[7][2] = this->locally->NdS[7][2];
		N[7][3] = this->locally->NdS[7][3];
	}//CTR + K; CTR + => C - dodaje lub U - usuwa
	double*size = new double[4];
	size = size_length();
	for (int i = 0; i < 4; i++) size[i] = size[i] / 2; // wyznacznik
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			NNT1[i][j] = N[0][i] * N[0][j];
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			NNT2[i][j] = N[1][i] * N[1][j];
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			NNT3[i][j] = N[2][i] * N[2][j];
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			NNT4[i][j] = N[3][i] * N[3][j];
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			NNT5[i][j] = N[4][i] * N[4][j];
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			NNT6[i][j] = N[5][i] * N[5][j];
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			NNT7[i][j] = N[6][i] * N[6][j];
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			NNT8[i][j] = N[7][i] * N[7][j];
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			suma1[i][j] = (NNT1[i][j] + NNT2[i][j]);
			suma1[i][j] *= size[0];
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			suma2[i][j] = (NNT3[i][j] + NNT4[i][j]);
			suma2[i][j] *= size[1];
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			suma3[i][j] = (NNT5[i][j] + NNT6[i][j]);
			suma3[i][j] *= size[2];
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			suma4[i][j] = (NNT7[i][j] + NNT8[i][j]);
			suma4[i][j] *= size[3];
		}
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			matrixH[i][j] = (suma1[i][j] + suma2[i][j] + suma3[i][j] + suma4[i][j])*alfa;
		}
	}
	return matrixH;
}
double** Grid::matrixH_full(int id_elem){
	double **dV = new double*[4];
	double **dS = new double*[4];
	double**H = new double*[4];
	for (int i = 0; i < 4; i++) {
		H[i] = new double[4];
		dV[i] = new double[4];
		dS[i] = new double[4];
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			H[i][j] = 0;
			dV[i][j] = 0;
			dS[i][j] = 0;
		}
	}
	dV = create_matrixH_dV(id_elem);
	dS = create_matrixH_dS(id_elem);
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			H[i][j] = dV[i][j] + dS[i][j];
		}
	}
	return H;
}
double**Grid::create_globalH() {
	double**localH = new double*[4];
	for (int i = 0; i < 4; i++){
		localH[i] = new double[4];
	}
	for (int i = 0; i < num_elements; i++) {
		localH = matrixH_full(i+1);
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				global_matrixH[arr_elem[i].id_node[j] - 1][arr_elem[i].id_node[k] - 1] += localH[j][k];
			}
		}
	}
	return global_matrixH;
}
double**Grid::create_globalC() {
	double**localC = new double*[4];
	for (int i = 0; i < 4; i++) {
		localC[i] = new double[4];
	} 
	for (int i = 0; i < num_elements; i++) {
		localC = create_matrixC(i + 1);
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				global_matrixC[arr_elem[i].id_node[j] - 1][arr_elem[i].id_node[k] - 1] += localC[j][k];
			}
		}
	}
	return global_matrixC;
}
double** Grid::matrixC_dev_dtau() {
	double **C = new double*[(nH*nL)];
	for (int i = 0; i < (nH*nL); i++) {
		C[i] = new double[(nH*nL)];
	}
	C = create_globalC();
	for (int i = 0; i < (nH*nL); i++) {
		for (int j = 0; j < (nH*nL); j++) {
			C[i][j] = (C[i][j]) / step_time;
		}
	}
	return C;
}
double**Grid::finally_H() {
	double**H = new double*[(nH*nL)];
	for (int i = 0; i < (nH*nL); i++) {
		H[i] = new double[(nH*nL)];
	}
	H = create_globalH();
	double**C = new double*[(nH*nL)];
	for (int i = 0; i < (nH*nL); i++) {
		C[i] = new double[(nH*nL)];
	}
	C = matrixC_dev_dtau();
	double**finallyH = new double*[(nH*nL)];
	for (int i = 0; i < (nH*nL); i++) {
		finallyH[i] = new double[(nH*nL)];
	}

	for (int i = 0; i < (nH*nL); i++) {
		for (int j = 0; j < (nH*nL); j++) {
			finallyH[i][j] = 0;
		}
	}

	for (int i = 0; i < (nH*nL); i++) {
		for (int j = 0; j < (nH*nL); j++) {
			finallyH[i][j] = H[i][j] + C[i][j];
		}
	}

	for (int i = 0; i < (nH*nL); i++) {
		for (int j = 0; j < (nH*nL); j++) {
			H[i][j] = 0;
		}
	}
	for (int i = 0; i < (nH*nL); i++) {
		for (int j = 0; j < (nH*nL); j++) {
			C[i][j] = 0;
		}
	}
	return finallyH;
}
double* Grid::create_vectorP(int id_elem) {
	double*vector_P = new double[4];
	for (int i = 0; i < 4; i++) {
		vector_P[i] = 0;
	}
	double** N = new double*[8];
	for (int i = 0; i < 8; i++) {
		N[i] = new double[4];
	}
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 4; j++) {
			N[i][j] = 0;
		}
	}
	if (this->arr_nodes[this->arr_elem[id_elem - 1].id_node[0] - 1].edge == 1 && this->arr_nodes[this->arr_elem[id_elem - 1].id_node[1] - 1].edge == 1) {	//warunek brzegowy spe³niony 
		N[0][0] = this->locally->NdS[0][0];
		N[0][1] = this->locally->NdS[0][1];
		N[0][2] = this->locally->NdS[0][2];
		N[0][3] = this->locally->NdS[0][3];
		N[1][0] = this->locally->NdS[1][0];
		N[1][1] = this->locally->NdS[1][1];
		N[1][2] = this->locally->NdS[1][2];
		N[1][3] = this->locally->NdS[1][3];
	}
	if (this->arr_nodes[this->arr_elem[id_elem - 1].id_node[1] - 1].edge == 1 && this->arr_nodes[this->arr_elem[id_elem - 1].id_node[2] - 1].edge == 1) {	//warunek brzegowy spe³niony 
		N[2][0] = this->locally->NdS[2][0];
		N[2][1] = this->locally->NdS[2][1];
		N[2][2] = this->locally->NdS[2][2];
		N[2][3] = this->locally->NdS[2][3];
		N[3][0] = this->locally->NdS[3][0];
		N[3][1] = this->locally->NdS[3][1];
		N[3][2] = this->locally->NdS[3][2];
		N[3][3] = this->locally->NdS[3][3];
	}
	if (this->arr_nodes[this->arr_elem[id_elem - 1].id_node[2] - 1].edge == 1 && this->arr_nodes[this->arr_elem[id_elem - 1].id_node[3] - 1].edge == 1) {	//warunek brzegowy spe³niony 
		N[4][0] = this->locally->NdS[4][0];
		N[4][1] = this->locally->NdS[4][1];
		N[4][2] = this->locally->NdS[4][2];
		N[4][3] = this->locally->NdS[4][3];
		N[5][0] = this->locally->NdS[5][0];
		N[5][1] = this->locally->NdS[5][1];
		N[5][2] = this->locally->NdS[5][2];
		N[5][3] = this->locally->NdS[5][3];
	}
	if (this->arr_nodes[this->arr_elem[id_elem - 1].id_node[3] - 1].edge == 1 && this->arr_nodes[this->arr_elem[id_elem - 1].id_node[0] - 1].edge == 1) {	//warunek brzegowy spe³niony 
		N[6][0] = this->locally->NdS[6][0];
		N[6][1] = this->locally->NdS[6][1];
		N[6][2] = this->locally->NdS[6][2];
		N[6][3] = this->locally->NdS[6][3];
		N[7][0] = this->locally->NdS[7][0];
		N[7][1] = this->locally->NdS[7][1];
		N[7][2] = this->locally->NdS[7][2];
		N[7][3] = this->locally->NdS[7][3];
	}
	double*suma1 = new double[4];
	double*suma2 = new double[4];
	double*suma3 = new double[4];
	double*suma4 = new double[4];
	for (int i = 0; i < 4; i++) {
		suma1[i] = 0;
		suma2[i] = 0;
		suma3[i] = 0;
		suma4[i] = 0;
	}
	double*size = new double[4];
	size = size_length();
	for (int i = 0; i < 4; i++) size[i] = size[i] / 2; // wyznacznik
	for (int i = 0; i < 4; i++) {
		suma1[i] = (N[0][i] + N[1][i])*size[0];
	}
	for (int i = 0; i < 4; i++) {
		suma2[i] = (N[2][i] + N[3][i])*size[1];
	}
	for (int i = 0; i < 4; i++) {
		suma3[i] = (N[4][i] + N[5][i])*size[2];
	}
	for (int i = 0; i < 4; i++) {
		suma4[i] = (N[6][i] + N[7][i])*size[3];
	}
	for (int i = 0; i < 4; i++) {
		vector_P[i] = (suma1[i] + suma2[i] + suma3[i] + suma4[i])*alfa;
		vector_P[i] *= t_env;
	}

	return vector_P;
}
double* Grid::create_globalP() {
	double*localP = new double[4];
	for (int i = 0; i < 4; i++) {
		localP[i] = 0;
	}
	for (int i = 0; i < num_elements; i++) {
		localP = create_vectorP(i + 1);
		for (int j = 0; j < 4; j++) {
			global_vectorP[arr_elem[i].id_node[j] - 1] += localP[j];
		}
	}
	return global_vectorP;
}
double*Grid::finally_P(double *temp) {
	double*P = new double[(nH*nL)];
	P = create_globalP();
	double**C = new double*[(nH*nL)];
	for (int i = 0; i < (nH*nL); i++) {
		C[i] = new double[(nH*nL)];
	}
	C = matrixC_dev_dtau();
	double*vectorC_mul_t0 = new double[(nH*nL)];
	for (int i = 0; i < (nH*nL); i++) {
		vectorC_mul_t0[i] = 0;
	}
	double*T = new double[(nH*nL)];
	for (int i = 0; i < (nH*nL); i++) {
		T[i] = temp[i];
	}

	for (int i = 0; i < (nH*nL); i++) {
		for (int j = 0; j < (nH*nL); j++) {
			vectorC_mul_t0[i] += (C[i][j] * T[j]);
		}
	}
	
	double*finallyP = new double[(nH*nL)];
	for (int i = 0; i < (nH*nL); i++) {
		finallyP[i] = 0;
	}
	for (int i = 0; i < (nH*nL); i++) {
		finallyP[i] = P[i] + vectorC_mul_t0[i];
	}
	for (int i = 0; i < (nH*nL); i++) {
		P[i] = 0;
	}
	for (int i = 0; i < (nH*nL); i++) {
		for (int j = 0; j < (nH*nL); j++) {
			C[i][j]=0 ;
		}
	}
	for (int i = 0; i < (nH*nL); i++) {
		T[i] = 0;
	}
	
	for (int j = 0; j < (nH*nL); j++) {
		vectorC_mul_t0[j] = 0;
	}
	return finallyP;
}
double* Grid::gauss( double ** HP, double * T){
	const double eps = 1e-12;
	int i, j, k;
	double m, s; // eliminacja wspó³czynników
	for (i = 0; i < num_nodes - 1; i++){
		for (j = i + 1; j < num_nodes; j++){
			if (fabs(HP[i][i]) < eps) return NULL;
			m = -HP[j][i] / HP[i][i];
			for (k = i + 1; k <= num_nodes; k++)
				HP[j][k] += m * HP[i][k];
		}
	}// wyliczanie niewiadomych
	for (i = num_nodes - 1; i >= 0; i--){
		s = HP[i][num_nodes];
		for (j = num_nodes - 1; j >= i + 1; j--)
			s -= HP[i][j] * T[j];
		if (fabs(HP[i][i]) < eps) return NULL;
		T[i] = s / HP[i][i];
	}
	return T;
}