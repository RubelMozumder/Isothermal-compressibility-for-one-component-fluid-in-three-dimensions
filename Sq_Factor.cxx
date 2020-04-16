//This part for calculating Static Structure Factor
//Using Kz = Ky = Kz = Nkx * 2*PI /L; Nkx = 30 
# include <iostream>
# include <fstream>
# include <string>
# include <vector>
# include <cmath>
# include <sstream>

using namespace std;

// Struct for K and its corresponding Sk
struct staticFactor{
	double	abs_k;
	double	Sk;
	staticFactor () : abs_k (0), Sk (0){}
	staticFactor ( staticFactor const& S )
		: abs_k ( S.abs_k )
		, Sk (S.Sk) {} 
};

//Struct three dimensional Vector
struct coordinate_type{
	double x ;
	double y ;
	double z ;
	coordinate_type () : x (0), y (0), z (0) {}
	// copy - constructor
	coordinate_type ( coordinate_type const & c )
		: x ( c.x )
		, y ( c.y ) 
		, z ( c.z ) {}
};

void readPositionFile ( string fname, int N_Part, vector < coordinate_type >& position );
void arrangeK_vecComp ( vector < coordinate_type >& K_vec, const int Nkx, const double del_k );
void staticFactorCal ( vector < coordinate_type > const& position, vector < staticFactor >& Sta_Fac_i_file, vector <coordinate_type > K_vec, const double del_k );
void Write_staticFactor_Data( string fname, vector < staticFactor > const& Sta_Fac );


int main( int argc , char* argv [] ){
 
    int N_Part = 12 * 12 * 12;
    const double density = 0.2; // in sigma cube unit
    const double L = pow ( N_Part / density, 1./3. ); 
    
    const int Nkx = 30 ; // Nkx = Nky = Nkz for K-vector componet
    
    const double del_k = double( 2 * M_PI / L ); //del_k = Kmin
    const double Kmax = pow (3, 1./2.) * Nkx * del_k ;
    //No of bin for K-vector
    const int N_bin =  int( ( Kmax - del_k ) / del_k ) ; 
    
    //Vector for averaged Static Factor over few time Steps
    vector < staticFactor > Sta_Fac_avg ( N_bin ) ;
    
    for ( int i = 0; i < N_bin; i++ ){
            
            Sta_Fac_avg [ i ]. abs_k = ( i + 1 ) * del_k ;
    };
    
    //  Vector For all combination of nx,ny and nz and ignoring the combination ( 0 0 0 )
    vector < coordinate_type > K_vec ( Nkx * Nkx * Nkx - 1 ); 
    //  Filling up the K Vector Component combination
    arrangeK_vecComp ( K_vec, Nkx, del_k );

    
    int nfiles = argc -1 ;
    cout << " Caluculating S(K) from " << nfiles << " files " << endl ;
    for ( int j = 1; j < nfiles +1; ++ j ) {
       
        vector <coordinate_type> position ;
        
        readPositionFile( argv [ j ], N_Part, position ) ;
        //Introducing Vector for Static Factor for each inputed file
        vector <staticFactor> Sta_Fac_i_file ( N_bin ) ;
        staticFactorCal ( position, Sta_Fac_i_file, K_vec, del_k );
             
        // Averaging for of Static Factor
        for ( int i = 0; i < N_bin; i++ ){
            
            Sta_Fac_avg [ i ]. Sk += Sta_Fac_i_file [ i ] . Sk / nfiles;
        }
        
        
        
    };
    
    stringstream fname ;
    fname << "StaticFactor_density" << density;
    
    Write_staticFactor_Data( fname.str(), Sta_Fac_avg );
    
};

void readPositionFile ( string fname, int N_Part, vector < coordinate_type >& position )
{
	double vx, vy, vz ;
	ifstream input( fname ) ;
	for( int i = 0; i < N_Part; i++){
		coordinate_type p ;
		input >> p.x >> p.y >> p.z >> vx >> vy >> vz ;
		position.push_back( p ) ;
	};
};

void arrangeK_vecComp ( vector < coordinate_type >& K_vec, const int Nkx, const double del_k )
{
    //for qubic system Nx=Ny=Nz
    for( int i = 0; i < Nkx; i++){
        for( int j = 0; j < Nkx; j++){
            for( int l = 0; l < Nkx; l++){
                if( i == 0 && j == 0 && l == 0)// test & & &
                    continue;
                // [here index - 1] as we skipt 0,0,0 combination
                K_vec [ Nkx*Nkx*i + Nkx*j + l - 1].x = l * del_k ;
                K_vec [ Nkx*Nkx*i + Nkx*j + l - 1].y = j * del_k ;
                K_vec [ Nkx*Nkx*i + Nkx*j + l - 1].z = i * del_k ;
            };
        };
    };
};

void staticFactorCal ( vector < coordinate_type > const& position, vector < staticFactor >& Sta_Fac_i_file, vector <coordinate_type > K_vec, const double del_k )
{
    const int N_K_vec = K_vec.size () ; //Total combination of Nkx,Nky, Nkz
    const int N_Part = position.size () ;
    const int N_bin = Sta_Fac_i_file.size () ;
    
     // For Counting the number of event in the K bin to average Sk
    vector < int > eventNo( N_bin ) ;
    for( int i = 0; i < N_bin; i++){
        
        eventNo [i] = 0 ;
    };
    
    // for loop for N_k_vec 
    for( int i = 0; i < N_K_vec; i++ ){
        
        double Sk_cos = 0 ;
        double Sk_sin = 0 ;
        
        coordinate_type const& K = K_vec [i] ;
        for( int j = 0; j < N_Part; j++ ) {
        
            coordinate_type const& P = position [ j ] ;
            Sk_cos += cos( K.x * P.x + K.y * P.y + K.z * P.z) ; // Test cos and sin value
            Sk_sin += sin( K.x * P.x + K.y * P.y + K.z * P.z) ;
        }
        
        double SkValue = ( Sk_cos * Sk_cos + Sk_sin * Sk_sin ) / N_Part ; //test double/int
        double abs_k_value = sqrt( K.x * K.x + K.y * K.y + K.z * K.z ) ; //test sqrt
        
        //minimum abs_k_value will del_k not 0 So K_bin_index is not < 1;
        int K_bin_index = int ( abs_k_value / del_k ) ;
        
        Sta_Fac_i_file [ K_bin_index - 1 ]. abs_k += abs_k_value ;
        Sta_Fac_i_file [ K_bin_index - 1 ]. Sk += SkValue ;
       
        eventNo [ K_bin_index - 1 ] += 1; 
        
    };
    
    for( int i = 0; i < N_bin; i++ ){
    
        if ( eventNo [ i ] == 0 )
            continue;
        
        Sta_Fac_i_file [ i ] . abs_k /= eventNo [ i ];
        Sta_Fac_i_file [ i ] . Sk /= eventNo [ i ];
        
    };
    
}

void Write_staticFactor_Data( string fname, vector < staticFactor > const & Sta_Fac )
{
    ofstream output( fname ) ;
    
    // the following lines set the output precision very high
    output. precision (15) ;
    output. setf ( ios :: scientific , ios :: floatfield ) ;
    
    const int N = Sta_Fac. size () ;
   
    for ( int i = 0; i < N ; ++ i ) {
        
        // i for indication line number in the output file
        output << Sta_Fac [ i ]. abs_k <<"\t"<< Sta_Fac [ i ]. Sk<<endl;// <<"\ti\n";
    
    };
    
};




















