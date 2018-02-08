	// This is my first C++ program

	//It prints a line of text
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>       /* sin */
#include <vector>
# include <iostream>
using std::vector;
typedef std::vector<double> double_vector;

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
struct phantom {//struct('h', 0, 'k', 0,'a', 10, 'b', 10, 'alpha', 0);
		double h; // x-coordinate of center [cm]
		double k; // y-coordinate of center [cm]
		double a; // x semi-axis [cm]
		double b; // y semi-axis [cm]
		double alpha; // rotational angle of ellipse [radians]
	};

// parallelBeam projector function declaration	
	void parallelBeamProjection(scanProtocol scannerInfo, phantom phantomInfo, double **sinogram);

	void multiply(vector<double> x,vector<double> yDet,  double rotationMatrix[2][2],vector<double> x_rot,vector<double> y_rot);

	int main()
	{
		scanProtocol GE_scanner = {1000, 1000, 1};
		phantom circular = {0, 0, 10, 10, 0};

	// scan
	// create sinogram ~2d array
	 // initialize 2d pointer to pointer (~array)
		int row = 10, col = 5;
		double **sinogram;
	sinogram = new double*[GE_scanner.Ns]; // dynamic array (size 10) of pointers to int
	for (int i = 0; i < GE_scanner.Ns; ++i) {
		sinogram[i] = new double[GE_scanner.NTheta];
	  // each i-th pointer is now pointing to dynamic array (size 10) of actual int values
	}

	parallelBeamProjection(GE_scanner, circular, sinogram);
	std:: cout << GE_scanner.Ns << std:: endl;

	//...
	size_t size = 10;

std::vector<double> array = linspace(0, 10, size);    // make room for 10 integers,
                                 // and initialize them to 0
// do something with them:
for(int i=0; i<size; ++i){
	array[i] = i;
	std::cout << array[i];
}
// no need to delete anything

return 0;  
} 

void parallelBeamProjection(scanProtocol scannerInfo, phantom phantomInfo,double **sinogram){
	double rotationMatrix [2][2];
	double_vector x = linspace(0, 10, scannerInfo.Ns);    // make room for 10 integers,
	std::vector<double> yDet(scannerInfo.Ns,-1500), x_detector(scannerInfo.Ns,0),y_detector(scannerInfo.Ns,0); 
	std::vector<double> ySource(scannerInfo.Ns,1500),x_source(scannerInfo.Ns,0), y_source(scannerInfo.Ns,0); 

	double_vector theta = linspace(1e-5, 2*M_PI - 1e-5, scannerInfo.Ns);    // make room for 10 integers,



		// loop through different view angles
	for (int i=0; i<scannerInfo.NTheta;i+=1){
   			// rotMatrix = [cos(theta(iTheta)) -sin(theta(iTheta));
        				// sin(theta(iTheta)) cos(theta(iTheta))];
		rotationMatrix[0][0] = cos(theta[i]);
		rotationMatrix[0][1] = -sin(theta[i]);
		rotationMatrix[1][0] = sin(theta[i]);
		rotationMatrix[1][1] = cos(theta[i]);

		multiply(x, yDet, rotationMatrix, x_detector, y_detector, scannerInfo);
		multiply(x, ySource, rotationMatrix, x_source, y_source, scannerInfo);


	}

	return;
}

	void multiply(vector<double> x,vector<double> yDet,  double rotationMatrix[2][2],vector<double> x_rot,vector<double> y_rot, scanProtocol scannerInfo){
		// 

		return;

	};
