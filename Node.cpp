#include "Node.h"

Node::Node() {
	this->x = 0;
	this->y = 0;
}

Node::Node(double x, double y) {
	this->x = x;
	this->y = y;
}
Node::~Node(){}

bool Node::getEdge() {
	return this->edge;
}


void Node::setEdge(bool edge) {
	this->edge = edge;
}