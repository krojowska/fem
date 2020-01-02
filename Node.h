#pragma once
class Node {

public:
	int id;
	double x, y;
	bool edge;

	Node();
	Node(double, double);
	~Node();
	bool getEdge();
	void setEdge(bool);
};