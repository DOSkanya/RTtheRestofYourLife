#pragma once
#include "hittable.h"
#include "matrix.h"
#include "vec3.h"

class triangle : public hittable {
public:
	triangle() {}
	triangle(point3 a, point3 b, point3 c, shared_ptr<material> m) : mat_ptr(m) {
		v[0] = a;
		v[1] = b;
		v[2] = c;
	}

	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const;
	virtual bool bounding_box(double time0, double tiem1, aabb& output_box) const;
	//virtual double pdf_value(const point3& o, const vec3& v) const;
	//virtual vec3 random(const vec3& o) const;
public:
	point3 v[3];
	shared_ptr<material> mat_ptr;
};

bool triangle::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
	matrix A;
	A.m[0][0] = v[0].x() - v[1].x();
	A.m[0][1] = v[0].x() - v[2].x();
	A.m[0][2] = r.direction().x();
	A.m[1][0] = v[0].y() - v[1].y();
	A.m[1][1] = v[0].y() - v[2].y();
	A.m[1][2] = r.direction().y();
	A.m[2][0] = v[0].z() - v[1].z();
	A.m[2][1] = v[0].z() - v[2].z();
	A.m[2][2] = r.direction().z();
	double A_d = A.determinant();

	matrix beta_m;
	beta_m.m[0][0] = v[0].x() - r.origin().x();
	beta_m.m[0][1] = v[0].x() - v[2].x();
	beta_m.m[0][2] = r.direction().x();
	beta_m.m[1][0] = v[0].y() - r.origin().y();
	beta_m.m[1][1] = v[0].y() - v[2].y();
	beta_m.m[1][2] = r.direction().y();
	beta_m.m[2][0] = v[0].z() - r.origin().z();
	beta_m.m[2][1] = v[0].z() - v[2].z();
	beta_m.m[2][2] = r.direction().z();
	double beta_d = beta_m.determinant();

	matrix gamma_m;
	gamma_m.m[0][0] = v[0].x() - v[1].x();
	gamma_m.m[0][1] = v[0].x() - r.origin().x();
	gamma_m.m[0][2] = r.direction().x();
	gamma_m.m[1][0] = v[0].y() - v[1].y();
	gamma_m.m[1][1] = v[0].y() - r.origin().y();
	gamma_m.m[1][2] = r.direction().y();
	gamma_m.m[2][0] = v[0].z() - v[1].z();
	gamma_m.m[2][1] = v[0].z() - r.origin().z();
	gamma_m.m[2][2] = r.direction().z();
	double gamma_d = gamma_m.determinant();

	double beta = beta_d / A_d;
	double gamma = gamma_d / A_d;
	double alpha = 1 - beta - gamma;
	if (alpha < 0 || beta < 0 || gamma < 0 || A_d == 0)
		return false;

	matrix t_m;
	t_m.m[0][0] = v[0].x() - v[1].x();
	t_m.m[0][1] = v[0].x() - v[2].x();
	t_m.m[0][2] = v[0].x() - r.origin().x();
	t_m.m[1][0] = v[0].y() - v[1].y();
	t_m.m[1][1] = v[0].y() - v[2].y();
	t_m.m[1][2] = v[0].y() - r.origin().y();
	t_m.m[2][0] = v[0].z() - v[1].z();
	t_m.m[2][1] = v[0].z() - v[2].z();
	t_m.m[2][2] = v[0].z() - r.origin().z();
	double t_d = t_m.determinant();
	double t = t_d / A_d;

	rec.p = r.at(t);
	rec.t = t;
	vec3 outward_normal = unit_vector(cross(v[1] - v[0], v[2] - v[0]));//Normalized
	rec.set_face_normal(r, outward_normal);
	rec.u = 0.0;
	rec.v = 0.0;
	rec.mat_ptr = mat_ptr;
	return true;
}

bool triangle::bounding_box(double time0, double tiem1, aabb& output_box) const {
	output_box = aabb(point3(std::min(v[0].x(), std::min(v[1].x(), v[2].x())),
		std::min(v[0].y(), std::min(v[1].y(), v[2].y())),
		std::min(v[0].z(), std::min(v[1].z(), v[2].z()))),
		point3(std::max(v[0].x(), std::max(v[1].x(), v[2].x())),
			std::max(v[0].y(), std::max(v[1].y(), v[2].y())),
			std::max(v[0].z(), std::max(v[1].z(), v[2].z()))));
	return true;
}
/*
double triangle::pdf_value(const point3& o, const vec3& v) const {
	hit_record rec;
	if (!this->hit(ray(o, v), 0.001, infinity, rec))
		return 0;

	vec3 a = this->v[1] - this->v[0];
	vec3 b = this->v[2] - this->v[0];
	auto area = cross(a, b).length() / 2.0;
	auto distance_squared = rec.t * rec.t * v.length_squared();
	auto cosine = fabs(dot(v, rec.normal) / v.length());

	return distance_squared / (cosine * area);
}

vec3 triangle::random(const vec3& o) const {
	auto alpha = random_double();
	auto beta = random_double(0.0, 1.0 - alpha);
	auto gamma = 1 - alpha - beta;

	auto random_point = alpha * v[0] + beta * v[1] + gamma * v[2];
	return random_point - o;
}*/