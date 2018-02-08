#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>       /* sin */
#include <vector>
# include <iostream>
#include <fstream>
using std::vector;
typedef std::vector<double> double_vector;
int ind = 499;
// linspace equivalent of MATLAB, slightly different implementation from Python
vector<double> linspace(double a, double b, int n) {
	vector<double> array;
	double step = (b-a) / (n-1);
	while(a <= b) {
		array.push_back(a);
        a += step;           // could recode to better handle rounding errors
    }
    return array;
}

// scan protocol strucutre
struct scanProtocol {
		int Ns; // number of samples
		int NTheta; // number of detectors
		double mmPerSample;
	} mine, yours;

// phantom description
struct phantom {//struct('h', 0, 'k', 0,'a', ind, 'b', ind, 'alpha', 0);
		double h; // x-coordinate of center [cm]
		double k; // y-coordinate of center [cm]
		double a; // x semi-axis [cm]
		double b; // y semi-axis [cm]
		double alpha; // rotational angle of ellipse [radians]
	};

// parallelBeam projector function declaration	
	void parallelBeamProjection(scanProtocol scannerInfo, phantom phantomInfo, double **sinogram);
	void multiply(vector<double> x,vector<double> y,  double rotationMatrix[2][2],vector<double> & x_rot,vector<double> & y_rot, scanProtocol scannerInfo);

	void ray_phantom_intersection(vector<double> & x1, vector<double> & y1, vector<double> & x2, vector<double> & y2, phantom phantomInfo,vector<double> p,vector<double> q, scanProtocol scannerInfo) ;

	void find_x_y_intercept(vector<double> & p, vector<double> & q, vector<double> x1, vector<double> y1, vector<double> x2, vector<double> y2, scanProtocol scannerInfo);


	int main(){
		scanProtocol GE_scanner = {1000, 1000, 1};
		phantom circular = {0, 0, 100, 100, 0};

		double **sinogram;
		sinogram = new double*[GE_scanner.Ns]; // dynamic array (size ind) of pointers to int
		for (int i = 0; i < GE_scanner.Ns; ++i) {
			sinogram[i] = new double[GE_scanner.NTheta];
	  		// each i-th pointer is now pointing to dynamic array (size ind) of actual int values
		}

		parallelBeamProjection(GE_scanner, circular, sinogram);


		// std::ofstream myfile;
		// myfile.open ("sinogram.csv");
		// int j;
		// for (int i=0; i<GE_scanner.NTheta; i++){
		// 	for (j=0; j<GE_scanner.Ns ; j++){
		// 		myfile << sinogram[i][j] << std::endl;
		// 	}
		// }
		// myfile.close();
		return 0;  
	}	 

void multiply(vector<double> x,vector<double> y,  double rotationMatrix[2][2],vector<double> & x_rot,vector<double> & y_rot, scanProtocol scannerInfo){
	for (int i = 0; i<scannerInfo.Ns; i+=1){
		x_rot[i] = x[i]* rotationMatrix[0][0]+ y[i]*rotationMatrix[1][0] ;
		y_rot[i] = x[i]* (-rotationMatrix[1][0])+ y[i]*rotationMatrix[0][0] ;

		}
	return;
	};

void ray_phantom_intersection(vector<double> & x1, vector<double> & y1, vector<double> & x2, vector<double> & y2, phantom phantomInfo,vector<double> p,vector<double> q, scanProtocol scannerInfo) { // out: x1, x2, y1, y2
	// % calculate x and y intercept for these rotated array
	double M, N, O, determinant;
    // parameter for line
	std::vector<double> m(scannerInfo.Ns,0),c(scannerInfo.Ns,0); 
	for (int i=0; i<scannerInfo.Ns; i+=1){
		m[i] = -q[i]/p[i];
		c[i] = q[i];
	}
	// parameter for ellipse
	double A = pow(cos(phantomInfo.alpha) /phantomInfo.a,2)  + pow(sin(phantomInfo.alpha)/ phantomInfo.b,2) ;
	double B = - 2 * cos (phantomInfo.alpha) * sin(phantomInfo.alpha) * (pow(phantomInfo.a,-2) - pow(phantomInfo.b,-2));
	double C = pow(sin(phantomInfo.alpha)/phantomInfo.a,2) + pow(cos(phantomInfo.alpha) / phantomInfo.b,2);


	for (int i = 0; i < scannerInfo.Ns; i+=1){
		// handle case for inifinte slope .i.e. q = inf
		if (isinf(q[i])){
			M = C;
			N = B*p[i]-(2*C*phantomInfo.k+B*phantomInfo.h);
			O = A*pow(p[i],2) - p[i]*(2*A*phantomInfo.h + phantomInfo.k*B) + (A*pow(phantomInfo.h,2)+B*phantomInfo.h*phantomInfo.k+C*pow(phantomInfo.k,2)-1);
			determinant = pow(N,2) - 4* M * O;
			if (determinant<0){
				y1[i]  = 0;
				y2[i]  = 0;
				x1[i] = 0;
				x2[i] = 0;
			} else{
				y1[i]  = (-N + pow(determinant,0.5))/ (2*M);
				y2[i]  = (-N - pow(determinant,0.5))/ (2*M);
				x1[i] = p[i];
				x2[i] = p[i];
				}
		} else {
			M =  A + B*m[i] + C* pow(m[i],2);
			N = B*c[i] + 2*C*m[i]*c[i] - 2*A*phantomInfo.h - phantomInfo.k*B - m[i]* (2*C*phantomInfo.k+ B* phantomInfo.h);
			O = C*pow(c[i],2) - c[i]* (2*C*phantomInfo.k + B*phantomInfo.h) + A*pow(phantomInfo.h,2) + B*phantomInfo.h*phantomInfo.k + C*pow(phantomInfo.k,2) -1 ;
			determinant = pow(N, 2) - 4 * O * M;
			if (determinant<0){
				y1[i]  = 0;
				y2[i]  = 0;
				x1[i] = 0;
				x2[i] = 0;
			} else{
				x1[i] = (- N + pow(determinant,0.5)) / (2*M);
				x2[i] = (- N - pow(determinant,0.5)) / (2*M);;
				y1[i] = m[i]*x1[i] + c[i];
				y2[i] = m[i]*x2[i] + c[i];
			}
		}
	}
};


void find_x_y_intercept(vector<double> & p, vector<double> & q, vector<double> x1, vector<double> y1, vector<double> x2, vector<double> y2, scanProtocol scannerInfo){
	for (int i = 0; i<scannerInfo.Ns; i+=1){
		q[i] = (y1[i]/x1[i] - y2[i]/x2[i])/(1/x1[i]-1/x2[i]);
		p[i] =1/(1/x1[i] - y1[i]/(x1[i]*q[i]));
	}
};  // information will be out in p, q

void parallelBeamProjection(scanProtocol scannerInfo, phantom phantomInfo,double **sinogram){
	double rotationMatrix [2][2];
	double_vector x = linspace(-scannerInfo.Ns / 2* scannerInfo.mmPerSample, scannerInfo.Ns /2 * scannerInfo.mmPerSample,scannerInfo.Ns);    // make room for ind integers,
	std::vector<double> yDet_o(scannerInfo.Ns,-1500), x_detector(scannerInfo.Ns,0),y_detector(scannerInfo.Ns,0); 
	std::vector<double> ySource_o(scannerInfo.Ns,1500),x_source(scannerInfo.Ns,0), y_source(scannerInfo.Ns,0); 
	std::vector<double> p(scannerInfo.Ns,0),q(scannerInfo.Ns,0); 
	std::vector<double> x1(scannerInfo.Ns,0),y1(scannerInfo.Ns,0),x2(scannerInfo.Ns,0),y2(scannerInfo.Ns,0); 
	std::vector<double> test(2,0.5); 
	double_vector theta = linspace(M_PI/2+ 1e-5, M_PI/2+2*M_PI - 1e-5, scannerInfo.Ns);    // make room for ind integers,

	// loop through different view angles
	for (int i=0; i<scannerInfo.NTheta;i+=1){
		rotationMatrix[0][0] = cos(theta[i]);
		rotationMatrix[1][0] = sin(theta[i]);
		multiply(x, yDet_o, rotationMatrix, x_detector, y_detector, scannerInfo);  // information will be stored in x_detector, y_detector
		multiply(x, ySource_o, rotationMatrix, x_source, y_source, scannerInfo);   // information will be stored in x_source, y_source
		find_x_y_intercept(p, q, x_detector, y_detector, x_source, y_source, scannerInfo);  // information will be out in p, q
        ray_phantom_intersection(x1,y1, x2, y2, phantomInfo, p, q, scannerInfo);  // out: x1, y1, x2, y2
        // assign the distance between (x1, y1) and (x2, y2) to the sinogram row and repeat for each view angle
        for (int j = 0; j<scannerInfo.Ns; j+=1){  	
        	sinogram[i][j] = pow(pow(x1[j]-x2[j],2) + pow(y1[j]-y2[j],2),0.5);
        }
    }
    return;
}