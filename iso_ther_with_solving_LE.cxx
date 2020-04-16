//Extrapolation Part to calculate Isothermal Compressibility
// Using the formula S(0) = S(k_1) + (dk)^2 * alpha + (dk)^4 + beta
# include <iostream>
# include <fstream>
# include <string>
# include <vector>
# include <cmath>
# include <sstream>

using namespace std;

//struct LinearData{
//	double alpha,beta,cost;
//};

int main( int argc, char* argv [] ){
	cout<<"This is the testing line"<<endl;

	vector <double> kvector(3);
	vector <double> static_factor(3);
	double IssothermalCompressibility;
	double density = 0;
	
	string fname = "Isso_Compr_LE.txt";
	ofstream output(fname);
	
 	int nfiles = argc -1;
    	cout << "Calculating Issothermal Compressibility  " << nfiles <<'\t'<< " differnt density value " << endl;

	for( int j = 1; j < nfiles + 1; ++ j ) {
        	string fname = argv [ j ];
        	ifstream input (fname);
		
		for(int i=0; i<3; i++){
			kvector[i] = 0;
			static_factor[i] = 0; 
			input >> kvector[i] >> static_factor[i];
		}
	density += 0.1;
	// S(k_2) = S(k_1) + (dk)^2 * Alpha + (dk)^4 * Beta
	// S(k_3) = S(k_2) + (dk)^2 * Alpha + (dk)^4 * Beta
	// Averaging the factor "(dk)^2 * Alpha + (dk)^4 * Beta"
	// (dk)^2 * Alpha + (dk)^4 * Beta = ( S(K_3) - S(k_1) )/2

	// S( k = 0) = S(k_1) - ( (dk)^2 * Alpha + (dk)^4 * Beta )
/*

	double static_factor_0 = static_factor[0] - ( static_factor[2] - static_factor[0] ) / (2*density);

*/
	double dk = kvector[1] - kvector[0];
	double alpha = (static_factor[2] - 2*static_factor[1] + static_factor[0]);
	double static_factor_0 = static_factor[0] - alpha;	

	output << density <<'\t'<< static_factor_0 << endl;
	
	};

	return 0;
};














