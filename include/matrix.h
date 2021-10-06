#pragma once
class matrix {
public:
	matrix() {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				m[i][j] = 0.0;
			}
		}
	}

	double determinant();
public:
	double m[3][3];
};

double matrix::determinant() {
	double determinant = m[0][0] * m[1][1] * m[2][2] + m[1][0] * m[2][1] * m[0][2] + m[0][1] * m[1][2] * m[2][0]
		- m[0][2] * m[1][1] * m[2][0] - m[1][0] * m[0][1] * m[2][2] - m[1][2] * m[2][1] * m[0][0];
	return determinant;
}