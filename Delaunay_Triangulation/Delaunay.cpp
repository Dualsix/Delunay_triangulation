#include <iostream>
#include <list>
#include "Delaunay.h"

std::vector<Edge> edges_list;

void triangulateMesh(Mesh &mesh){

	superTriangle(mesh);

	for (int i = 0; i < mesh.vertices.size() - 3; i++){
		addVertex(mesh, i);
		if (i % 100 == 0) std::cout << i << " of " << mesh.vertices.size() << " completed.\n";
	}

	for (int i = 0; i < mesh.faces.size(); i++){
		if ((mesh.faces[i].a == mesh.vertices.size() - 3) || (mesh.faces[i].a == mesh.vertices.size() - 2) || (mesh.faces[i].a == mesh.vertices.size() - 1) ||
			(mesh.faces[i].b == mesh.vertices.size() - 3) || (mesh.faces[i].b == mesh.vertices.size() - 2) || (mesh.faces[i].b == mesh.vertices.size() - 1) ||
			(mesh.faces[i].c == mesh.vertices.size() - 3) || (mesh.faces[i].c == mesh.vertices.size() - 2) || (mesh.faces[i].c == mesh.vertices.size() - 1))
		{
			mesh.faces.erase(mesh.faces.begin() + i);
			i = i - 1;
		}
	}
}

void superTriangle(Mesh &mesh){

	double
		xmin = mesh.vertices[0].x,
		xmax = xmin,
		ymin = mesh.vertices[0].y,
		ymax = ymin;

	for (int i = 1; i < mesh.vertices.size(); i++){
		if (mesh.vertices[i].x < xmin) xmin = mesh.vertices[i].x;
		if (mesh.vertices[i].x > xmax) xmax = mesh.vertices[i].x;
		if (mesh.vertices[i].y < ymin) ymin = mesh.vertices[i].y;
		if (mesh.vertices[i].y > ymax) ymax = mesh.vertices[i].y;
	}

	double
		dx = xmax - xmin,
		dy = ymax - ymin,
		dmax = (dx > dy) ? dx : dy,
		xmid = (xmax + xmin) / 2.0,
		ymid = (ymax + ymin) / 2.0;

	Vertex
		v1{ xmid - 20 * dmax, ymid - dmax, 0.0 },
		v2{ xmid + 20 * dmax, ymid - dmax, 0.0 },
		v3{ xmid, ymid + 20 * dmax, 0.0 };

	mesh.vertices.push_back(v1);
	mesh.vertices.push_back(v2);
	mesh.vertices.push_back(v3);
	mesh.faces.push_back(Face{ mesh.vertices.size() - 3, mesh.vertices.size() - 2, mesh.vertices.size() - 1 });

}

//Implement this function
void addVertex(Mesh &mesh, const int vertexIndex){
	for (int j = 0; j < mesh.faces.size(); j++) {
		if (isInsideCircumCircle(mesh.vertices[mesh.faces[j].a], mesh.vertices[mesh.faces[j].b], mesh.vertices[mesh.faces[j].c], mesh.vertices[vertexIndex])) {
			Edge aux;
			aux.v1 = mesh.faces[j].a;
			aux.v2 = mesh.faces[j].b;
			edges_list.push_back(aux);

			Edge aux1;
			aux1.v1 = mesh.faces[j].b;
			aux1.v2 = mesh.faces[j].c;
			edges_list.push_back(aux1);

			Edge aux2;
			aux2.v1 = mesh.faces[j].c;
			aux2.v2 = mesh.faces[j].a;
			edges_list.push_back(aux2);
			
			mesh.faces.erase(mesh.faces.begin() + j);
			j--;
		}
	}
	for (int  k = 0; k < edges_list.size(); k++) {
		Edge e1 = edges_list[k];
		for (int  j = k + 1; j < edges_list.size(); j++) {
			Edge e2 = edges_list[j];
			if (((e1.v1 == e2.v1) && (e1.v2 == e2.v2)) || ((e1.v2 == e2.v1) && (e1.v1 == e2.v2))){
				edges_list.erase(edges_list.begin() + j);
				edges_list.erase(edges_list.begin() + k);
				k--;
			}
		}
	}
	for (int j = 0; j < edges_list.size(); j++) {
		Face f;
		f.a = edges_list[j].v1;
		f.b = edges_list[j].v2;
		f.c = vertexIndex;
		mesh.faces.push_back(f);
	}
	edges_list.clear();
}


bool isInsideCircumCircle(const Vertex A, const Vertex B, const Vertex C, const Vertex &point){
	double m1, m2, mx1, mx2, my1, my2, xc, yc, r;
	double dx, dy, rsqr, drsqr;

	if (abs(A.y - B.y) < EPSILON && abs(B.y - C.y) < EPSILON)
		return(false);
	if (abs(B.y - A.y) < EPSILON){
		m2 = -(C.x - B.x) / (C.y - B.y);
		mx2 = (B.x + C.x) / 2.0;
		my2 = (B.y + C.y) / 2.0;
		xc = (B.x + A.x) / 2.0;
		yc = m2 * (xc - mx2) + my2;
	}
	else if (abs(C.y - B.y) < EPSILON){
		m1 = -(B.x - A.x) / (B.y - A.y);
		mx1 = (A.x + B.x) / 2.0;
		my1 = (A.y + B.y) / 2.0;
		xc = (C.x + B.x) / 2.0;
		yc = m1 * (xc - mx1) + my1;
	}
	else{
		m1 = -(B.x - A.x) / (B.y - A.y);
		m2 = -(C.x - B.x) / (C.y - B.y);
		mx1 = (A.x + B.x) / 2.0;
		mx2 = (B.x + C.x) / 2.0;
		my1 = (A.y + B.y) / 2.0;
		my2 = (B.y + C.y) / 2.0;
		xc = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2);
		yc = m1 * (xc - mx1) + my1;
	}
	dx = B.x - xc;
	dy = B.y - yc;
	rsqr = dx * dx + dy * dy;
	dx = point.x - xc;
	dy = point.y - yc;
	drsqr = dx * dx + dy * dy;
	return((drsqr <= rsqr) ? true : false);

}






