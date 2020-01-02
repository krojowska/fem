#include "LocalLayout.h"
#include <math.h>
#include<iostream>
using namespace std;

#define pos_cord 1/sqrt(3);
#define neg_cord -1/sqrt(3);

LocalLayout::LocalLayout() {
	for (int i = 0; i < 4; i++) {
		this->ksi[i] = 0;
		this->eta[i] = 0;
	}

	for (int i = 0; i < 8; i++) {
		this->ksi_S[i] = 0;
		this->eta_S[i] = 0;
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			this->NdV[i][j] = 0;
			this->dNdksi[i][j] = 0;
			this->dNdeta[i][j] = 0;
		}
	}

	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 4; j++) {
			this->NdS[i][j] = 0;
		}
	}

	set_ksi();
	set_eta();

	calc_NdV();
	calc_NdS();

	calc_dNdksi();
	calc_dNdeta();
}

LocalLayout::LocalLayout(double *ksi, double *eta, double *ksi_S, double*eta_S, double**NdV, double**NdS, double **dNdksi, double**dNdeta) {
	for (int i = 0; i < 4; i++) {
		this->ksi[i] = 0;
		this->eta[i] = 0;
	}

	for (int i = 0; i < 8; i++) {
		this->ksi_S[i] = 0;
		this->eta_S[i] = 0;
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			this->NdV[i][j] = 0;
			this->dNdksi[i][j] = 0;
			this->dNdeta[i][j] = 0;
		}
	}
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 4; j++) {
			this->NdS[i][j] = 0;
		}
	}

	set_ksi();
	set_eta();

	calc_NdV();
	calc_NdS();

	calc_dNdksi();
	calc_dNdeta();
}

LocalLayout::~LocalLayout() {};

void LocalLayout::set_ksi() {
	this->ksi[0] = neg_cord;
	this->ksi[1] = pos_cord;
	this->ksi[2] = pos_cord;
	this->ksi[3] = neg_cord;
}

void LocalLayout::set_eta() {
	this->eta[0] = neg_cord;
	this->eta[1] = neg_cord;
	this->eta[2] = pos_cord;
	this->eta[3] = pos_cord;
}

void LocalLayout::calc_NdV() {
	//set_ksi();
	//set_eta();

	for (int i = 0; i < 4; i++) {
		this->NdV[i][0] = 0.25*(1 - this->ksi[i])*(1 - this->eta[i]);
		this->NdV[i][1] = 0.25*(1 + this->ksi[i])*(1 - this->eta[i]);
		this->NdV[i][2] = 0.25*(1 + this->ksi[i])*(1 + this->eta[i]);
		this->NdV[i][3] = 0.25*(1 - this->ksi[i])*(1 + this->eta[i]);
	}

	/*cout << "N: " << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << N[i][j] << " ";
		}
		cout << endl;
	}*/
}

void LocalLayout::calc_NdS() {
	this->ksi_S[0] = neg_cord;
	this->ksi_S[1] = pos_cord;
	this->ksi_S[2] = 1;
	this->ksi_S[3] = 1;
	this->ksi_S[4] = pos_cord;
	this->ksi_S[5] = neg_cord;
	this->ksi_S[6] = -1;
	this->ksi_S[7] = -1;

	this->eta_S[0] = -1;
	this->eta_S[1] = -1;
	this->eta_S[2] = neg_cord;
	this->eta_S[3] = pos_cord;
	this->eta_S[4] = 1;
	this->eta_S[5] = 1;
	this->eta_S[6] = pos_cord;
	this->eta_S[7] = neg_cord;

	for (int i = 0; i < 8; i++) {
		this->NdS[i][0] = 0.25*(1 - this->ksi_S[i])*(1 - this->eta_S[i]);
		this->NdS[i][1] = 0.25*(1 + this->ksi_S[i])*(1 - this->eta_S[i]);
		this->NdS[i][2] = 0.25*(1 + this->ksi_S[i])*(1 + this->eta_S[i]);
		this->NdS[i][3] = 0.25*(1 - this->ksi_S[i])*(1 + this->eta_S[i]);

	}
}

void LocalLayout::calc_dNdksi() {
	//setEta();
	for (int i = 0; i < 4; i++) {
		this->dNdksi[i][0] = -0.25*(1 - this->eta[i]);
		this->dNdksi[i][1] = 0.25*(1 - this->eta[i]);
		this->dNdksi[i][2] = 0.25*(1 + this->eta[i]);
		this->dNdksi[i][3] = -0.25*(1 + this->eta[i]);
	}

	/*cout << endl<< "dndksi: " << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << dNdksi[i][j] << " ";
		}
		cout << endl;
	}*/

}
void LocalLayout::calc_dNdeta() {
	//setKsi();
	for (int i = 0; i < 4; i++) {
		this->dNdeta[i][0] = -0.25*(1 - this->ksi[i]);
		this->dNdeta[i][1] = -0.25*(1 + this->ksi[i]);
		this->dNdeta[i][2] = 0.25*(1 + this->ksi[i]);
		this->dNdeta[i][3] = 0.25*(1 - this->ksi[i]);
	}

	/*cout << endl<<"dndeta: " << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << dNdeta[i][j] << " ";
		}
		cout << endl;
	}*/
}

