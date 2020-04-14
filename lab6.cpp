#include<iostream>
#include<math.h>
const int N = 13;
using namespace std;
//Сплайн интерполяция т.к для метода 3/8 нужна ещё одна точка
double *gauss(double **a, double *y, int n){
	double *x = new double[n];
	for(int i = 0; i < n; i++){
		y[i] /= a[i][i];
		for(int k = n - 1; k >= 0; k--){
			a[i][k] /= a[i][i];
		}
		for(int j = i + 1; j < n; j++){
			y[j] -= a[j][i] * y[i];
			for(int k = n - 1; k >= 0; k--){
				a[j][k] -= a[j][i] * a[i][k];
			}
		}
	}
	for(int i = n - 1; i >= 0; i--){
		for(int j = i - 1; j >= 0; j--){
			y[j] -= a[j][i] * y[i];
			a[j][i] = 0;
		}
	}
	for(int i = 0; i < n; i++)
		x[i] = y[i];
	return x;
}
double *progon(double **a, double *y, int n){
	double *x = new double[n];
	const int N = n - 1;
	double alp[N];
	double bet[N + 1];
 	alp[0] = -a[0][1] / a[0][0];
 	bet[0] = y[0] / a[0][0];
	for(int i = 1; i < N; i++){
		alp[i] = -a[i][i + 1] / (a[i][i] + a[i][i - 1] * alp[i - 1]);
		bet[i] = (y[i] - a[i][i - 1] * bet[i - 1]) / (a[i][i] + a[i][i - 1] * alp[i - 1]); 
	}
	bet[N] = (y[N] - a[N][N - 1] * bet[N - 1]) / (a[N][N] + a[N][N - 1] * alp[N - 1]);
	x[N] = bet[N];
	for(int i = N - 1; i >= 0; i--){
		x[i] = alp[i] * x[i + 1] + bet[i];
	}
	return x;
} 
double *find_k(double fx[2][9]){
	const int n = 9;
	double** matrix = new double*[n];
	for(int i = 0; i < n; i++)
		matrix[i] = new double[n];
	double y[n];
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
				matrix[i][j] = 0;
				if(i == j)
					matrix[i][j] = 4;
				if(i == j + 1 || i == j - 1)
					matrix[i][j] = 1;
		}
	}
	for(int i = 1; i < 8; i++)
		y[i] = 3 * (fx[1][i + 1] - fx[1][i - 1]);
	y[0] = 3 * (fx[1][1] - fx[1][0]);
	y[8] = 3 * (fx[1][8] - fx[1][7]);
	matrix[0][0] = 2;
	matrix[n - 1][n - 1] = 2;
	return progon(matrix, y, n);	
}
double *find_a(int i0, double *k, double fx[2][9]){
	int n = 4;
	double **a = new double*[n];
	for(int i = 0; i < n; i++)
		a[i] = new double[n];
	double a0[4][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {1, 1, 1, 1}, {0, 1, 2, 3}};
	for(int i = 0; i < n; i++)
		for(int j = 0; j < n; j++)			
			a[i][j] = a0[i][j];
	double *y = new double[n];
	double y0[4] = {fx[1][i0], k[i0], fx[1][i0 + 1], k[i0 + 1]};
	for(int i = 0; i < n; i++)
		y[i] = y0[i];
	double *c = gauss(a, y, n);
	delete[] a;
	delete[] y;
	return c; 
}
double cub_spline(double* a, double x, double x0){
	double f = 0;
	for(int i = 0; i <= 3; i++){
		f += a[i] * pow(x - x0, i);
	}
	return f;
}
double trap(double fx[2][N], double h0, int k){ // трапеций
	double trap = 0;
	double h = k * h0;
	for(int i = 0; i < N - 1; i += k){
		trap += (fx[1][i] + fx[1][i + k]);
	}
	return h * trap / 2;
}
double simp(double fx[2][N], double h0, int k){	// Симпсона
	double simp = 0;
	double h = k * h0;
	for(int i = 0; i < N - 1; i += 2 * k){
		simp += fx[1][i] + 4 * fx[1][i + k] + fx[1][i + 2 * k];
	}
	return h * simp / 3;
}
double m3d8(double fx[2][N], double h0, int k){ // 3 / 8
	double m3d8 = 0;
	double h = k * h0;
	for(int i = 0; i < N - 1; i += 3 * k){
		m3d8 += fx[1][i] + 3 * (fx[1][i + k] + fx[1][i + 2 * k]) + fx[1][i + 3 * k];
	}
	return h * m3d8 * 3 / 8;
}
double runge(double(*foo)(double fx[2][N], double h0, int k), double fx[2][N], double h0){
	return abs(foo(fx, h0, 1) - foo(fx, h0, 2));
}
int main(){
	double h0 = 0.125;
	double fx[2][9] = {{0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1}, {0, 0.02147, 0.29305, 0.494105, 0.541341, 0.516855, 0.468617, 0.416531, 0.367879}};
	// С помощью сплайн интерполяции поделим функцию на N частей, чтобы можно было применить все методы в одинаковых условиях.
	double fxN[2][N];
	double nanoh = 0.00001;
	double h1 = 1. / (N - 1);
	int k = 0;
	for(int i = 0; i < 8; i++){
		for(double j = i * h0; j < (i + 1) * h0; j += nanoh){
			if(fabs(j - k * h1) <= nanoh){
				fxN[0][k] = j;
				fxN[1][k] = cub_spline(find_a(i, find_k(fx), fx), j, i * h0);
				k++;
			}
		}
	}
	cout<<"trap: "<<trap(fxN, h1, 1)<<endl;
	cout<<"runge trap: "<<runge(trap, fxN, h1) / 3<<endl;
	cout<<"simp: "<<simp(fxN, h1, 1)<<endl;
	cout<<"runge simp: "<<runge(simp, fxN, h1) / 15<<endl;
	cout<<"3 / 8: "<<m3d8(fxN, h1, 1)<<endl;
	cout<<"runge 3/8: "<<runge(m3d8, fxN, h1) / 31<<endl;
	return 0;
}