// Molecular Dynamic Part To Create equilibrium Phase Space
# include <fstream>
# include <string>
# include <vector>
# include <cmath>
# include <iostream>

using namespace std ;

double logistic_map (){//Used for creating initial velocity
	static double x = 0.5; // start value
	double const mu = 1.99;
	x = mu * x *(1 - x ) ;
	return x ;
};
struct coordinate_type{
	double x ;
	double y ;
	double z ;
	//Constructor
	coordinate_type () : x (0) , y (0) , z (0) {}
	// copy - constructor
	coordinate_type ( coordinate_type const & c ): x(c.x) , y (c.y) , z (c.z) {}

};
struct simulation_data{
	int N ;// number of particles
	double dt ;// delta t
	int step ; // current simulation step
	double potential_energy ; // per particle
	double kinetic_energy ;// per particle
	double temperature ;
	double pressure ;
	double virial ;
	bool with_thermostat ; // whether the termostat is enabled or not
	double thermostat_temperature ;
	double thermostat_coupling ; // coupling constant of the thermostat (known as Q )
	double zeta ; // friction coefficient
	coordinate_type total_momentum ;
	coordinate_type box_length ;
	vector < coordinate_type > position ;
	vector < coordinate_type > position_unfolded ;
	vector < coordinate_type > position0 ;
	vector < coordinate_type > velocity ;
	vector < coordinate_type > force ;
	struct potential_data {
	double rr_c ;// square of cutoff of the potential cutoff
	double shift_epot ; // shift of the potential
	} potential ;
	vector <coordinate_type> obstacle_pos ;
};

void initialize ( simulation_data & sim , int N ,const double density );
void set_initial_phase_space ( simulation_data & sim );
void write_phase_space ( simulation_data const & sim, const double density );
void pair_potential ( double rr , double & fval , double & epot , simulation_data const& sim );
void calculate_force ( simulation_data & sim );
//void read_phase_space ( simulation_data & sim , string fname );
void md_step ( simulation_data & sim );
void calc_thermodynamic_variables ( simulation_data & sim );


int main (){
	double Tmax ;
	simulation_data sim;
	const double density = 0.2;// density in sigma
	int N = 12*12*12;
	initialize (sim, N, density);
	write_phase_space ( sim, density ) ;
	cout << " dt : " ;
	cin >> sim . dt ;
	cout << " Tmax : " ;
	cin >> Tmax ;
	int steps = Tmax / sim . dt ;
	cout << " Simulating " << steps << " steps .. " ;
	cout << " with canonical ( NVT ) measurement " << endl ;
	calculate_force ( sim ) ;
	ofstream thermodynamics ( " thermodynamics " ) ;
	int print_interval = max ( steps /5000 , 1) ;
	sim.thermostat_coupling = 1.0 / sim . N ;
	sim.thermostat_temperature = 1.0;
	sim.with_thermostat = true ;
	for ( int i = 0; i < steps ; ++ i ) {
		if (( i % print_interval ) == 0 || i == 0) {
			calc_thermodynamic_variables ( sim ) ;
			double etot = sim.kinetic_energy + sim.potential_energy ;
			thermodynamics << i * sim.dt
			<< " " << sim.kinetic_energy
			<< " " << sim.potential_energy
			<< " " << etot
			<< " " << sim.pressure	
			<< " " << sim.total_momentum.x
			<< " " << sim.total_momentum.y
			<< " " << sim.total_momentum.z
			<< endl ;
		}

		sim . step ++ ;
		md_step ( sim ) ;
        if ( sim . step > ( steps - 5 ) ){
            write_phase_space ( sim, density ) ;
        
        }
            
	}
	// write out the final configuration
	
	return 0;
}

void initialize ( simulation_data & sim , const int N, const double density ){
	sim.step = 0;
	sim.N = N ;
	sim.temperature = 1.5;// temperature used for initialization
	sim.with_thermostat = false ; // disable thermostat at the beginning
	sim.thermostat_temperature = 1.0;
	sim.zeta = 0; // initial value for the friction
	// set up arrays
	sim.position.resize (N) ;
	sim.velocity.resize (N) ;
	sim.position0.resize (N) ;
	sim.position_unfolded.resize (N) ;
	sim.force.resize (N) ;
	const double edge_length = pow(N/density, 1./3.);
	sim.box_length.x = edge_length;
	sim.box_length.y = edge_length;
	sim.box_length.z = edge_length;
	// intitialize the potential
	sim.potential.shift_epot = 0;
	// calculate the potential energy shift
	double r_c = pow(2, 1./6.);
	sim.potential.rr_c = r_c * r_c ;
	double fval , epot_shift ;
	pair_potential ( r_c * r_c , fval , epot_shift , sim ) ;
	sim.potential.shift_epot = epot_shift ;
	// set the starting positions
	set_initial_phase_space (sim) ;
}

void set_initial_phase_space ( simulation_data& sim ){
	// lattice constant
	const double a = pow( sim.box_length.x*sim.box_length.y*sim.box_length.z/sim.N , 1./3. ) ;
	// place particles on lattice
	int Nx = 12;
	int Ny = 12;
	int Nz = 12;
	for ( int i = 0; i < Nx ; ++ i ) {
		for ( int j = 0; j < Ny ; ++ j ) {
			for (int k = 0; k < Nz; ++ k ){ 
				sim.position [ i + Nx*j + Nx*Ny*k ]. x = a *( i +0.5) ;
				sim.position [ i + Nx*j + Nx*Ny*k ]. y = a *( j +0.5) ;
				sim.position [ i + Nx*j + Nx*Ny*k ]. z = a *( k +0.5) ;
			}
		}
	}
	sim.total_momentum = coordinate_type() ; // reset to zero
	for ( int i = 0; i < sim.N ; ++i ) {
		coordinate_type& v = sim.velocity [i];
		v.x = logistic_map () - 0.5;
		v.y = logistic_map () - 0.5;
		v.z = logistic_map () - 0.5;
		sim.total_momentum.x += v.x ;
		sim.total_momentum.y += v.y ;
		sim.total_momentum.z += v.z ;
	}
	// shift total momentum to zero
	coordinate_type shift ;
	shift.x = sim.total_momentum.x / sim.N ;
	shift.y = sim.total_momentum.y / sim.N ;
	shift.z = sim.total_momentum.z / sim.N ;
	double vv = 0;
	for (int i = 0; i < sim.N ; ++i) {
		coordinate_type& v = sim.velocity [i];
		v.x -= shift.x ;
		v.y -= shift.y ;
		v.z -= shift.z ;
		// calculate new v ^2
		vv += v.x * v.x + v.y * v.y + v.z * v.z ;
	}
	// scale the velocities using the scaling factor, d = 3
	double scale_factor = sqrt (3* sim.N * sim.temperature / vv ) ;
	for ( int i = 0; i < sim . N ; ++ i ) {
		sim.velocity [i].x *= scale_factor ;
		sim.velocity [i].y *= scale_factor ;
		sim.velocity [i].z *= scale_factor ;
	}
};

void write_phase_space ( simulation_data const& sim, const double density ){
	ofstream output ( " Phase_space_ " + to_string ( sim.step ) + "_Density" + to_string ( density) ) ;
	// the following lines set the output precision very high
	output . precision (15) ;
	output . setf ( ios :: scientific , ios :: floatfield ) ;
	// write the phase space data
	for ( int i = 0; i < sim.N ; ++ i ) {
		output <<sim.position [i].x <<'\t'<< sim.position [i].y <<'\t'<< sim.position [i].z <<'\t'
		<< sim.velocity [i].x <<'\t'<< sim.velocity [i].y<<'\t'<< sim.velocity [i].z <<'\n';
	}
}

void pair_potential ( double rr , double & fval , double & epot , simulation_data const& sim ){
	// try to calculate the exponents of r in an effective manner :
	// reuse the expontentials of 1/ r ^2 for both 1/ r ^{12} and 1/ r ^{6}
	double const epsilon = 1;
	double rri = 1/rr ;// 1/ r ^2
	double r6i = rri * rri * rri ; // 1/ r ^6
	double eps_r6i = epsilon * r6i ;
	fval = 48 * rri * eps_r6i * ( r6i - 0.5) ;
	epot = 4 * eps_r6i * ( r6i - 1) - sim.potential.shift_epot ;
}

void calculate_force ( simulation_data & sim ){
	// reset the potential energy , as it will be recalcuated in this step
	double epot = 0;
	double virial = 0;
	// reset the force array
	for ( int i = 0; i < sim . N ; ++ i ) {
		sim.force[i]. x = 0;
		sim.force[i]. y = 0;
		sim.force[i]. z = 0;
	}
	for ( int i = 0; i < sim.N ; ++i ) {
		coordinate_type const& x1 = sim . position [ i ];
		// calculate the interaction beween two particles only once
		// ( ie . use newtons third law )
		for ( int j = i +1; j < sim.N ; ++j ) {
			double pot , fval ;
			coordinate_type const& x2 = sim . position [ j ];
			// distance between two particles
			coordinate_type dx ;
			dx . x = x1 . x - x2 . x ;
			dx . y = x1 . y - x2 . y ;
			dx . z = x1 . z - x2 . z ;
			// apply minimum image convention
			if ( dx . x > sim . box_length . x /2) {
				dx . x -= sim . box_length . x ;
			}
			if ( dx . x < - sim . box_length . x /2) {
				dx . x += sim . box_length . x ;
			}
			if ( dx . y > sim . box_length . y /2) {
				dx . y -= sim . box_length . y ;
			}
			if ( dx . y < - sim . box_length . y /2) {
				dx . y += sim . box_length . y ;
			}
			if ( dx . z > sim . box_length . z /2) {
				dx . z -= sim . box_length . z ;
			}
			if ( dx . z < - sim . box_length . z /2) {
				dx . z += sim . box_length . z ;
			}
			// calculate r ^2
			double rr = dx.x * dx.x + dx.y * dx.y + dx.z * dx.z ;
			// if we are over the cutoff , there is no force contribution
			if ( rr > sim.potential.rr_c )
				continue ;
			// calculate pair interaction
			pair_potential ( rr , fval , pot , sim ) ;
			epot += pot ;
			// convert F / r into vectorial value
			coordinate_type f ;
			f.x = fval * dx.x ;
			f.y = fval * dx.y ;
			f.z = fval * dx.z ;
			// set the force value using Newton ’s third law
			sim.force [ i ]. x += f.x ;
			sim.force [ i ]. y += f.y ;
			sim.force [ i ]. z += f.z ;

			sim.force [ j ]. x += -f.x ;
			sim.force [ j ]. y += -f.y ;
			sim.force [ j ]. z += -f.z ;
			virial += f.x * dx.x + f.y * dx.y + f.z * dx.z ;
			}
	}
	// we want the potential energy per particle
	sim . potential_energy = epot / sim . N ;
	sim . virial = virial ;
}

void md_step ( simulation_data & sim ){
	// this function assumes that the force has been calculated
	// the step before
	double const dt = sim.dt ;
	double const Q = sim.thermostat_coupling ;
	double inst_temp = 0; // instantaneous temperature
	double zeta_app = 0;
	double vel_factor = 1;
	if ( sim . with_thermostat ) {
		vel_factor = 1 - 0.5* dt * zeta_app ;
		zeta_app = sim . zeta + dt / Q *( sim.temperature - sim.thermostat_temperature );
	}
	// first half - step of the verlet integration
	// update the positions and the velocities
	for ( int i = 0; i < sim . N ; ++ i ) {
		coordinate_type & pos = sim . position [ i ];
		coordinate_type & pos_unfolded = sim . position_unfolded [ i ];
		coordinate_type & vel = sim . velocity [ i ];
		coordinate_type const& force = sim . force [ i ];
		// dfx and dfy is the total force -- including friction --
		// on the particle for dt/2
		double dfx = force . x - sim . zeta * vel . x ;
		double dfy = force . y - sim . zeta * vel . y ;
		double dfz = force . z - sim . zeta * vel . z ;
		double dx = vel . x * dt + 0.5* dfx * dt * dt ;
		double dy = vel . y * dt + 0.5* dfy * dt * dt ;
		double dz = vel . z * dt + 0.5* dfz * dt * dt ;
		pos . x += dx ;
		pos . y += dy ;
		pos . z += dz ;
		pos_unfolded . x += dx ;
		pos_unfolded . y += dy ;
		pos_unfolded . z += dz ;
		// velocity includes approximated friction coefficient
		vel . x += 0.5* dt *( dfx - zeta_app * vel . x ) * vel_factor ;
		vel . y += 0.5* dt *( dfy - zeta_app * vel . y ) * vel_factor ;
		vel . z += 0.5* dt *( dfz - zeta_app * vel . z ) * vel_factor ;
		// apply periodic boundary conditions
		if ( pos . x > sim . box_length . x ) {
		pos . x -= sim . box_length . x ;
		}
		if ( pos . x < 0) {
		pos . x += sim . box_length . x ;
		}
		if ( pos . y > sim . box_length . y ) {
		pos . y -= sim . box_length . y ;
		}
		if ( pos . y < 0) {
		pos . y += sim . box_length . y ;
		}
		if ( pos . z > sim . box_length . z ) {
		pos . z -= sim . box_length . z ;
		}
		if ( pos . z < 0) {
		pos . z += sim . box_length . z ;
		}
	}
	// update the forces as the positions have changed
	calculate_force ( sim ) ;
	// second half - step of the verlet integration
	// updates the velocity
	for ( int i = 0; i < sim . N ; ++ i ) {
		coordinate_type & vel = sim . velocity [ i ];
		coordinate_type const & force = sim . force [ i ];
		vel . x += 0.5* dt * force . x * vel_factor ;
		vel . y += 0.5* dt * force . y * vel_factor ;
		vel . z += 0.5* dt * force . z * vel_factor ;
		inst_temp += vel . x * vel . x + vel . y * vel . y + vel . z * vel . z ;
	}
	inst_temp /= 3.0* sim . N ;
	// update zeta
	if ( sim.with_thermostat ) {
		sim.zeta += dt /(2*Q ) *( sim.temperature + inst_temp - 2* sim.thermostat_temperature ) ;
	}
	else {
		sim . zeta = 0;
	}
	sim . temperature = inst_temp ;
}

void calc_thermodynamic_variables ( simulation_data& sim ){
	double ekin = 0;
	sim.total_momentum = coordinate_type () ;
	for ( int i = 0; i < sim . N ; ++ i ) {
		coordinate_type const& vel = sim.velocity [i];
		// assume unit mass
		ekin += vel.x * vel.x + vel.y * vel.y + vel.z * vel.z ;
		sim.total_momentum. x += vel.x ;
		sim.total_momentum. y += vel.y ;
		sim.total_momentum. z += vel.z ;
	}
	// missing factor of two , and divide by the number of particles
	ekin /= 2.* sim . N ;
	sim.kinetic_energy = ekin ;
	sim.temperature = ekin ;
	double volume = sim . box_length . x * sim . box_length . y * sim . box_length . z ;
	sim.pressure = ( sim.temperature * sim.N +  sim.virial / 3. ) / volume ;
};
