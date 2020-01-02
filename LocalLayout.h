#pragma once
class LocalLayout {

public:
	double ksi[4];
	double eta[4];

	double ksi_S[8];
	double eta_S[8];

	double NdV[4][4]; //funckje ksztaltu objetosc
	double NdS[8][4]; //funkcje ksztaltu po powierzchni
	double dNdksi[4][4];
	double dNdeta[4][4];

	LocalLayout();
	LocalLayout(double *, double *, double*, double*, double**, double**, double **, double**);
	~LocalLayout();

	void calc_NdV();
	void calc_NdS();
	void calc_dNdksi();
	void calc_dNdeta();
	void set_ksi();
	void set_eta();
};