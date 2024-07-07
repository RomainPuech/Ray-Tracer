#include <cmath>
#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <iostream>
#include <algorithm>

#define PI 3.1415926535
#define GAMMA 2.2

#include <string>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <list>

#include <chrono>
using namespace std::chrono;


#include <random>
static std::default_random_engine engine(10);
static std::uniform_real_distribution<double> uniform(0,1);

void boxMuller(double stdev , double &x , double &y) {
	double r1 = uniform(engine) ;
	double r2 = uniform(engine) ;
	x = sqrt(-2 * log ( r1 ) ) *cos ( 2 * PI*r2 ) *stdev ;
	y = sqrt(-2 * log ( r1 ) ) *sin ( 2 * PI*r2 ) *stdev ;
}

class Vector {
public:
	explicit Vector(double x = 0, double y = 0, double z = 0) {
		data[0] = x;
		data[1] = y;
		data[2] = z;
	}
	double norm2() const {
		return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
	}
	double norm() const {
		return sqrt(norm2());
	}
	void normalize() {
		double n = norm();
		data[0] /= n;
		data[1] /= n;
		data[2] /= n;
	}
	Vector cross(const Vector& b) const {
		return Vector(data[1] * b[2] - data[2] * b[1],
			data[2] * b[0] - data[0] * b[2],
			data[0] * b[1] - data[1] * b[0]);
	}
	double operator[](int i) const { return data[i]; };
	double& operator[](int i) { return data[i]; };
	double data[3];
};

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator-(const Vector& a) {
	return Vector(-a[0], -a[1], -a[2]);
}
Vector operator*(const double a, const Vector& b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
	return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(const Vector& a, const Vector& b) {
	return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}
Vector operator/(const Vector& a, const double b) {
	return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}


// random cosine

Vector random_cos(const Vector &N){
	// find smallest element of N
	int min_index = 0;
	for(int i = 1; i<3; i++){
		if(std::abs(N[i])<std::abs(N[min_index])){
			min_index = i;
		}
	}
	Vector T1 = Vector(0.,0.,0.);
	if(min_index==0){
		T1 = Vector(0.,N[2],-N[1]);
		
	}else if(min_index==1){
		T1 = Vector(N[2],0.,-N[0]);
	}else{
		T1 = Vector(N[1],-N[0],0.);
	}
	T1.normalize();
	Vector T2 = N.cross(T1);
	T2.normalize();
	double r1 = (double)rand() / RAND_MAX;
	double r2 = (double)rand() / RAND_MAX;
	double theta = 2.*PI*r1;
	return cos(theta)*sqrt(1-r2)*T1 + sin(theta)*sqrt(1-r2)*T2 + sqrt(r2)*N;

}
/////////////////////

class Ray {
public:
	Vector u,O;
	explicit Ray(Vector u, Vector O):u(u),O(O){};
};


// geometry
class Geometry {
public:
	Vector rho;
	bool mirror;
	Geometry() {};
	Geometry(Vector rho):rho(rho),mirror(false){};
	Geometry(Vector rho, bool mirror):rho(rho),mirror(mirror){};
	virtual double intersect(const Ray& r, Vector &res, Vector &resN, Vector &resUV, double current_best_obj=0.) = 0;
};
class Sphere : public Geometry{
public:
	Vector C;
	double R;
	explicit Sphere(){};
	explicit Sphere(Vector C, double R, Vector rho):Geometry(rho,false),C(C),R(R){};
	explicit Sphere(Vector C, double R, Vector rho, bool mirror):Geometry(rho,mirror),C(C),R(R){};

	double intersect(const Ray& r, Vector &res, Vector &resN, Vector &resUV, double current_best = 0.){
		double delta = pow(dot(r.u, r.O - C),2) - ((r.O - C).norm2() - pow(R,2));
		if(delta<0){return 0;}
		double t0 = dot(r.u,C-r.O);
		double sqrtdelta = std::sqrt(delta);
		double t1 = t0 - sqrtdelta;
		double t2 = t0 + sqrtdelta;
		double t = 0;
		if(t1>0){
			t = t1;
		}else if(t2>0){
			t = t2;
		}
		// results
		res = r.O + t*r.u;
		resN = res - C;
		resUV = Vector(-1.,-1.,-1.);
		return t;
	}
};

// Mesh
class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
	};
	int vtxi, vtxj, vtxk; // indices for vertex coordinates array
	int uvi, uvj, uvk;  // uv coordinates array
	int ni, nj, nk;  // normals array
	int group;      
};

class BoundingBox {
public:
    Vector B_min;
    Vector B_max;
    explicit BoundingBox(Vector min = Vector(), Vector max = Vector()):B_min(min),B_max(max) {}
	bool intersect(const Ray& r, double &t_res){
		double txmin = (B_min[0] - r.O[0]) / r.u[0];
		double txmax = (B_max[0] - r.O[0]) / r.u[0];
		if (txmin > txmax) std::swap(txmin, txmax);
		double tymin = (B_min[1] - r.O[1]) / r.u[1];
		double tymax = (B_max[1] - r.O[1]) / r.u[1];
		if (tymin > tymax) std::swap(tymin, tymax);
		double tzmin = (B_min[2] - r.O[2]) / r.u[2];
		double tzmax = (B_max[2] - r.O[2]) / r.u[2];
		if (tzmin > tzmax) std::swap(tzmin, tzmax);

		if(std::min(txmax,std::min(tymax,tzmax))>std::max(txmin,std::max(tymin,tzmin))){
			// there is an intersection
			double res_t = std::max(txmin,std::max(tymin,tzmin));
			//res = r.O + res_t*r.u;
			t_res = res_t;
			return true;
			}
		return false;
	}
};

struct Node {
	BoundingBox bbox;
	Node* child_left = nullptr;
	Node* child_right = nullptr;
	int starting_triangle;
	int ending_triangle;
};

class TriangleMesh : public Geometry{
private:
	Node* root;
	unsigned char* texture;
	int texture_H;
	int texture_W;

	BoundingBox compute_bbox(int starting_triangle = 0, int ending_triangle = -1) {
		if (ending_triangle == -1) ending_triangle = indices.size();
		double infinity = std::numeric_limits<double>::max();
		double neg_infinity = -std::numeric_limits<double>::max();
		Vector min = Vector(infinity,infinity,infinity);
		Vector max = Vector(neg_infinity,neg_infinity,neg_infinity);
		for(int i = starting_triangle; i<ending_triangle; i++){
			Vector A = vertices[indices[i].vtxi];
			Vector B = vertices[indices[i].vtxj];
			Vector C = vertices[indices[i].vtxk];
			min[0] = std::min(min[0],std::min(A[0],std::min(B[0],C[0])));
			min[1] = std::min(min[1],std::min(A[1],std::min(B[1],C[1])));
			min[2] = std::min(min[2],std::min(A[2],std::min(B[2],C[2])));
			max[0] = std::max(max[0],std::max(A[0],std::max(B[0],C[0])));
			max[1] = std::max(max[1],std::max(A[1],std::max(B[1],C[1])));
			max[2] = std::max(max[2],std::max(A[2],std::max(B[2],C[2])));
		}
		BoundingBox res = BoundingBox(min,max);
		return res;

	}

	Vector compute_diag(const BoundingBox& bbox) {
		return bbox.B_max - bbox.B_min;
	}

	int get_longest(const Vector& diag) {
		if (std::abs(diag[0]) > std::abs(diag[1]) && std::abs(diag[0]) > std::abs(diag[2])) return 0;
		if (std::abs(diag[1]) > std::abs(diag[0]) && std::abs(diag[1]) > std::abs(diag[2])) return 1;
		return 2;
	}

	Vector compute_barycenter(const TriangleIndices& t) {
		Vector A = vertices[t.vtxi];
		Vector B = vertices[t.vtxj];
		Vector C = vertices[t.vtxk];
		return (A + B + C) / 3.;
	}

	void compute_BVH(Node* node, int starting_triangle = 0, int ending_triangle = -1) {
		if (ending_triangle == -1) ending_triangle = indices.size();
		node->bbox = compute_bbox(starting_triangle, ending_triangle); // BBox from starting triangle included to ending triangle excluded
		node->starting_triangle = starting_triangle;
		node->ending_triangle = ending_triangle;
		Vector diag = compute_diag(node->bbox);
		Vector middle_diag = node->bbox.B_min + diag * 0.5;
		int longest_axis = get_longest(diag);
		int pivot_index = starting_triangle;
		for (int i = starting_triangle; i < ending_triangle; i++) {
			Vector barycenter = compute_barycenter(indices[i]);

			
			if (barycenter[longest_axis] < middle_diag[longest_axis]) {
				std::swap(indices[i], indices[pivot_index]);
				pivot_index++;
			}
		}
		// stopping criterion
		if (pivot_index <= starting_triangle || pivot_index >= ending_triangle-1 || ending_triangle - starting_triangle < 5){
			return;
		}
		
		node->child_left = new Node();
		node->child_right = new Node();
		this->compute_BVH(node->child_left, starting_triangle, pivot_index);
		this->compute_BVH(node->child_right, pivot_index, ending_triangle);
	}

public:
	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
    ~TriangleMesh() {}
    TriangleMesh():Geometry(Vector(1.,1.,1.),false) {};
	TriangleMesh(Vector rho):Geometry(rho,false){};
	TriangleMesh(Vector rho, bool mirror):Geometry(rho,mirror){};

	double triangle_intersect(const Ray& r, Vector &res, Vector &resN, Vector &resUV, int starting_triangle = 0, int ending_triangle = -1){
		if (ending_triangle == -1) ending_triangle = indices.size();
		double t = 0.;
		for(int i = starting_triangle; i < ending_triangle; i++){
			Vector A = vertices[indices[i].vtxi];
			Vector B = vertices[indices[i].vtxj];
			Vector C = vertices[indices[i].vtxk];
			Vector e1 = B - A;
			Vector e2 = C - A;
			Vector N = cross(e1,e2); // not normalized!
			Vector u = r.u;
			Vector O = r.O;
			Vector AO = A - O;

			double beta = dot(e2,cross(AO,u))/dot(u,N);
			double gamma = -dot(e1,cross(AO,u))/dot(u,N);
			double alpha = 1 - beta - gamma;
			double tprime = dot(AO,N)/dot(u,N);
			if((alpha>=0 and beta>=0 and gamma>=0) and tprime>0 and (tprime < t or t==0)){
				t = tprime;
				res = r.O + t*r.u;
				resN = N;
				// new shading normal:
				resN =  (alpha*normals[indices[i].ni] + beta*normals[indices[i].nj] + gamma*normals[indices[i].nk]);
				resUV = alpha*uvs[indices[i].uvi] + beta*uvs[indices[i].uvj] + gamma*uvs[indices[i].uvk];
				resUV[0] = std::fmod(resUV[0],1.);
				resUV[1] = std::fmod(resUV[1],1.);
				
			}

		}
		return t;
	}
	double intersect(const Ray& r, Vector &res, Vector &resN, Vector &resUV, double current_best_object = 0.) {
		
		double inter_distance = 0;
		Vector res_tmp = Vector(0.,0.,0.);
		Vector resN_tmp = Vector(0.,0.,0.);
		Vector resUV_tmp = Vector(-1.,-1.,-1.);
		if(!root->bbox.intersect(r,inter_distance)){
			
			return 0;
		}
		/////////////////////////////////
		if(inter_distance>current_best_object and current_best_object!=0.) return 0;
		std::list<Node*> nodes_to_visit;
		nodes_to_visit.push_front(root);
		
		double best_inter_distance = std::numeric_limits<double>::max();
		while (!nodes_to_visit.empty()) {
			Node* curNode = nodes_to_visit.back();
			nodes_to_visit.pop_back();
			// if there is one child, then it is not a leaf, so test the bounding box
			if (curNode->child_left) {
				if (curNode->child_left->bbox.intersect(r, inter_distance)) {
					if (inter_distance < best_inter_distance) {
						nodes_to_visit.push_back(curNode->child_left);
					}
				}
			}
			if(curNode->child_right){
				if (curNode->child_right->bbox.intersect(r, inter_distance)) {
					if (inter_distance < best_inter_distance) {
						nodes_to_visit.push_back(curNode->child_right);
					}
				}
			}
			else {
				// if there is no child, then it is a leaf, so test all triangles 
				
				inter_distance = triangle_intersect(r, res_tmp, resN_tmp, resUV_tmp, curNode->starting_triangle, curNode->ending_triangle);
				if (inter_distance > 0 and inter_distance < best_inter_distance) {
					best_inter_distance = inter_distance;
					res = res_tmp;
					resN = resN_tmp;
					resUV = resUV_tmp;
				}
				
			}
		}
		// set the texture
		if(resUV[0]>=0 and resUV[1]>=0){
			
			int x = (int) ((resUV[0]) * (texture_W));
			int y = (int) ((resUV[1]) * (texture_H));
			y = texture_H - y - 1;
			
			int index = 3*(y*texture_W + x);
			
			resUV = Vector(pow((double)texture[index]/255.,GAMMA),pow((double)texture[index+1]/255.,GAMMA),pow((double)texture[index+2]/255.,GAMMA));
			
			
		}
		return best_inter_distance;
	}

	void readOBJ(const char* obj, double scale=1, Vector translation = Vector(0.,0.,0.)) {

		char matfile[255];
		char grp[255];

		FILE* f;
		f = fopen(obj, "r");
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				curGroup++;
			}

			if (line[0] == 'v' && line[1] == ' ') {
				Vector vec;

				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
					col[0] = std::min(1., std::max(0., col[0]));
					col[1] = std::min(1., std::max(0., col[1]));
					col[2] = std::min(1., std::max(0., col[2]));

					vertices.push_back(vec*scale + translation);
					vertexcolors.push_back(col);

				} else {
					sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
					vertices.push_back(vec*scale + translation);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector vec;
				sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;
				t.group = curGroup;

				char* consumedline = line + 1;
				int offset;

				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
					indices.push_back(t);
				} else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
						indices.push_back(t);
					} else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							indices.push_back(t);
						} else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}

				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.group = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					} else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						} else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;								
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							} else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								} else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}

			}

		}
		fclose(f);
		root = new Node();
		this->root = root;
		compute_BVH(root);

		// load the texture
		this->texture = stbi_load("cat_diff.png", &texture_W, &texture_H, NULL, 3);
	}
	
};

class Light{
	public:
	Vector L;
	double I;
	explicit Light(Vector L, double I):L(L),I(I){};

};

class Scene{
	std::vector<Geometry*> scene;
	public:
		Light l;
		explicit Scene(Light l):l(l){};
		void add(Geometry* s){
			scene.push_back(s);
		};
		bool intersect(Ray r, Vector &res, Vector &resN, Vector &resUV, Geometry**&obj){
			double t = 0;
			double tmpt=0;
			Vector tmpres = Vector();
			Vector tmpN = Vector();
			Vector tmpUV = Vector(-1.,-1.,-1.);
			for (auto it = begin(scene); it != end(scene); ++it) {
				tmpt = (*it)->intersect(r,tmpres,tmpN,tmpUV,t);
    			if(tmpt>0 and (t==0 or tmpt < t)){
					t = tmpt;
					*obj = *it;
					res[0] = tmpres[0];
					res[1] = tmpres[1];
					res[2] = tmpres[2];
					resN = tmpN;
					resUV = tmpUV;
				}
			}
			if(t==0){
				return false;
			}
			return true;
		}
		void del(){
			for(auto it = scene.begin(); it!=scene.end(); ++it){
				delete *it;
			}
		}
};

bool getShadow(const Vector& P, const Vector& N, Scene s){
			// N must be normalized 
			Vector P_epsilon = P + 0.1*N;
			Light l = s.l;
			Vector directions_light = l.L - P_epsilon;
			directions_light.normalize();
			Ray to_light = Ray(directions_light,P_epsilon);
			Vector P_prime = Vector(0.,0.,0.);
			Vector N_prime = Vector(0.,0.,0.);
			Vector UV_prime = Vector(-1.,-1.,-1.);

			Geometry* geomPtr = new Sphere();
    		Geometry** sphere_prime = &geomPtr;

			if(s.intersect(to_light,P_prime,N_prime,UV_prime,sphere_prime)){
				double dist1 = (P_prime-P_epsilon).norm2();
				double dist2 = (l.L-P_epsilon).norm2();
				if(dist1<=dist2){
					return true;
				}
			}
			return false;
		}

Vector getColor(Ray r, Scene s, int rec){
	if(rec<0){return Vector(0.,0.,0.);}
	Light l = s.l;
	Vector P = Vector(0,0,0);
	Vector N = Vector(0,0,0);
	Vector UV = Vector(-1.,-1.,-1.);
	Geometry* geomPtr = new Sphere();
    Geometry** objectptr = &geomPtr;
	if(!s.intersect(r,P,N,UV,objectptr)){
		return Vector(0.,0.,0.);
	}
	Geometry* object = (*objectptr); 
	N.normalize();
	if(object->mirror){
		Vector P_epsilon = P + 0.001*N;
		Vector iota = r.u;
		Vector r = iota-2*dot(iota,N)*N;
		r.normalize();
		Ray newray = Ray(r,P_epsilon);
		return getColor(newray,s,rec-1);
	}
	double distance_attenuation = l.I/(4*PI*(l.L - P).norm2());
	Vector material = object->rho;
	// texture
	if(UV[0]!=-1.){
		material = UV;

	}
	double solid_angle = std::max(0.,dot(N,(l.L-P)/(l.L-P).norm()));
	bool visibility = (double)!getShadow(P,N,s);
	Vector L0 = Vector(0.,0.,0.);
	L0 = distance_attenuation*solid_angle*material/PI*visibility;
	// indirect lighting
	Vector P_epsilon = P + 0.00001*N;
	Ray randomRay = Ray(random_cos(N),P_epsilon);
	Vector L1 = material*getColor(randomRay,s,rec-1);
	L0 = L0 + L1;
	return L0;
};

int main() {

	int W = 512;
	int H = 512;

	// camera center
	Vector Q = Vector(0,0,55);
	const double alpha = PI/3;

	//Define the scene
	// central ball
	Vector C = Vector(0,0,0);
	double R = 10;
	Vector rho = Vector(0.5,0.5,0.5);
	Sphere* central_sphere = new Sphere(C,R,rho);
	// mirror ball
	Vector mirror_C = Vector(0,0,0);
	Sphere* mirror_sphere = new Sphere(mirror_C,R,rho,true);
	//
	Vector CG = Vector(0,0,-1000);
	double RG = 940;
	Vector rhoG = Vector(0.,1.,0.);
	Sphere* sphereG = new Sphere(CG,RG,rhoG);
	//
	Vector CB = Vector(0,-1000,0);
	double RB = 990;
	Vector rhoB = Vector(0.,0.,1.);
	Sphere* sphereB = new Sphere(CB,RB,rhoB);
	//
	Vector CR = Vector(0,1000,0);
	double RR = 940;
	Vector rhoR = Vector(1.,0.,0.);
	Sphere* sphereR = new Sphere(CR,RR,rhoR);
	//
	Vector CP = Vector(0,0,1000);
	double RP = 940;
	Vector rhoP = Vector(1.,0.10,0.75);
	Sphere* sphereP = new Sphere(CP,RP,rhoP);
	//
	TriangleMesh* cat = new TriangleMesh();
	cat->readOBJ("cat.obj",0.6,Vector(0.,-10.,0));

	// define light

	Vector L = Vector(-10,20,40); 
	double I = 2e10;
	Light l = Light(L,I);
	Scene scene = Scene(l);
	//scene.add(central_sphere);
	scene.add(mirror_sphere);
	scene.add(sphereG);
	scene.add(sphereB);
	scene.add(sphereR);
	scene.add(sphereP);
	//scene.add(cat);

	// start clock now to compute rendering time only
	auto start = high_resolution_clock::now();

	double x,y;
	std::vector<unsigned char> image(W * H * 3, 0);
	#pragma omp parallel for
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			Vector color = Vector(0.,0.,0.);
			int num_rays = 1;
			double x,y;
			for(int raynb = 0; raynb<num_rays;++raynb){
				boxMuller(0.4, x, y);
                Vector pixel = Vector(Q[0]+(j+x)+0.5-W/2,
                                      Q[1]-(i+y)-0.5+H/2,
                                      Q[2]-W/(2*tan(alpha/2)));
				Vector u = (pixel - Q);
				u.normalize();
                Ray ray = Ray(u,Q);
				color = color + getColor(ray,scene,4);
			}
			color = color/num_rays;
			image[(i * W + j) * 3 + 0] = std::min(pow(color[0],1/GAMMA),255.);
			image[(i * W + j) * 3 + 1] = std::min(pow(color[1],1/GAMMA),255.);
			image[(i * W + j) * 3 + 2] = std::min(pow(color[2],1/GAMMA),255.);
		
			
		}
	}
	// rendering stopped
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	//std::cout << (int) duration.count()/1000 << " microseconds"<< std::endl;
	stbi_write_png("image.png", W, H, 3, &image[0], 0);
	scene.del();
	return 0;
}




