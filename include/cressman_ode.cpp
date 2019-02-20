#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <boost/numeric/odeint.hpp>


typedef std::vector<double> state_type;


class Cressman {
    public:

        Cressman( const double Cm, const double GNa, const double GK, const double GAHP,
            const double GKL, const double GNaL, const double GClL, const double GCa, const double Gglia,
            const double Koinf, const double gamma1, const double tau, const double control)
            : par_Cm(Cm), par_GNa(GNa), par_GK(GK), par_GAHP(GAHP), par_GKL(GKL), par_GNaL(GNaL),
              par_GClL(GClL), par_GCa(GCa), par_Gglia(Gglia), par_Koinf(Koinf),
              par_gamma1(gamma1), par_tau(tau), par_control(control) { }

        template<class State>
        void operator() (const State &x, State &dxdt, const double /* t */ ) {
            const double Ipump = rho*(1.0/(1.0 + exp((25.0 - x[6])/3.0)))*(1.0/(1.0 + exp(5.5 - x[5])));
            const double IGlia = par_Gglia/(1.0 + exp((18.0 - x[5])/2.5));
            const double Idiff = eps0*(x[5] - par_Koinf);

            const double Ki = 140.0 + (18.0 - x[6]);
            const double Nao = 144.0 - beta0*(x[6] - 18.0);
            const double ENa = 26.64*log(Nao/x[6]);
            const double EK = 26.64*log(par_control*x[5]/Ki);

            const double a_m = (3.0 + 0.1*x[0])/(1.0 - exp(-3.0 - 1.0/10*x[0]));
            const double b_m = 4*exp(-55.0/18 - 1.0/18.0*x[0]);
            const double ah = 0.07*exp(-11.0/5.0 - 1.0/20.0*x[0]);
            const double bh = 1.0/(1.0 + exp(-7.0/5.0 - 1.0/10.0*x[0]));
            const double an = (0.34 + 0.01*x[0])/(1.0 - exp(-17.0/5.0 - 1.0/10.0*x[0]));
            const double bn = 0.125*exp(-11.0/20.0 - 1.0/80.0*x[0]);
            const double minf = a_m/(a_m + b_m);
            const double taum = 1.0/(a_m + b_m);
            const double h_inf = ah/(bh + ah);
            const double tauh = 1.0/(bh + ah);
            const double ninf = an/(an + bn);
            const double taun = 1.0/(an + bn);

            const double INa = par_GNa*std::pow(x[1], 3)*x[3]*(x[0] - ENa) + par_GNaL*(x[0] - ENa);
            const double IK = (par_GK*std::pow(x[2], 4) + par_GAHP*x[4]/(1.0 + x[4]) + par_GKL)*(x[0] - EK);
            const double ICl = par_GClL*(x[0]-ECl);

            dxdt[0] = -(INa + IK + ICl)/par_Cm;
            dxdt[1] = phi*(minf - x[1])/taum;
            dxdt[2] = phi*(ninf - x[2])/taun;
            dxdt[3] = phi*(h_inf - x[3])/tauh;
            dxdt[4] = -x[4]/80.0 - 0.002*par_GCa*(x[0] - ECa)/(1.0 + exp(-(x[0] + 25.0)/2.5));
            dxdt[5] = (par_gamma1*beta0*IK - 2*beta0*Ipump - IGlia - Idiff)/par_tau;
            dxdt[6] = -(par_gamma1*INa + 3*Ipump)/par_tau;
        }
    private:
        const double rho = 1.25;
        const double eps0 = 1.2;
        const double beta0 = 7;
        const double ECa = 120;
        const double Cli = 6;
        const double Clo = 130;
        const double ECl = 26.64*log(Cli/Clo);
        const double phi = 3;

        const double par_Cm;
        const double par_GNa;
        const double par_GK;
        const double par_GAHP;
        const double par_GKL;
        const double par_GNaL;
        const double par_GClL;
        const double par_GCa;
        const double par_Gglia;
        const double par_Koinf;
        const double par_gamma1;
        const double par_tau;
        const double par_control;
};

jac = [
    (7*par_gamma1*(par_GKL + x[2]^4*par_GK + (x[3]*par_GAHP)/(x[4] + 1)))/par_tau , 0 , (28*x[2]^3*par_GK*par_gamma1*(x[0] - (666*log(-(x[5]*par_control)/(x[6] - 158)))/25))/par_tau , 0 , (7*par_gamma1*(par_GAHP/(x[4] + 1) - (x[4]*par_GAHP)/(x[4] + 1)^2)*(x[0] - (666*log(-(x[5]*par_control)/(x[6] - 158)))/25))/par_tau , -((4662*par_gamma1*(par_GKL + x[2]^4*par_GK + (x[4]*par_GAHP)/(x[4] + 1)))/(25*x[5]) + (35*exp(11/2 - x[5]))/(2*(exp(11/2 - x[5]) + 1)^2*(exp(25/3 - x[6]/3) + 1)) + (2*par_Gglia*exp(36/5 - (2*x[5])/5))/(5*(exp(36/5 - (2*x[5])/5) + 1)^2) + 6/5)/par_tau , -((35*exp(25/3 - x[6]/3))/(6*(exp(11/2 - x[5]) + 1)*(exp(25/3 - x[6]/3) + 1)^2) - (4662*par_gamma1*(par_GKL + x[2]^4*par_GK + (x[4]*par_GAHP)/(x[4] + 1)))/(25*(x[6] - 158)))/par_tau ;
    -(par_gamma1*(par_GNaL + x[3]*x[1]^3*par_GNa))/par_tau , -(3*x[3]*x[1]^2*par_GNa*par_gamma1*(x[0] - (666*log(-(7*x[6] - 270)/x[6]))/25))/par_tau , 0 , -(x[1]^3*par_GNa*par_gamma1*(x[0] - (666*log(-(7*x[6] - 270)/x[6]))/25))/par_tau , 0 , -(15*exp(11/2 - x[5]))/(4*par_tau*(exp(11/2 - x[5]) + 1)^2*(exp(25/3 - x[6]/3) + 1)) , -(par_gamma1*((666*x[6]*par_GNaL*((7*x[6] - 270)/x[6]^2 - 7/x[6]))/(25*(7*x[6] - 270)) + (666*x[3]*x[1]^3*x[6]*par_GNa*((7*x[6] - 270)/x[6]^2 - 7/x[6]))/(25*(7*x[6] - 270))) + (5*exp(25/3 - x[6]/3))/(4*(exp(11/2 - x[5]) + 1)*(exp(25/3 - x[6]/3) + 1)^2))/par_tau ];

class Cressman_jacobi {
    public:
        Cressman_jacobi(){}

        template<class State, class Matrix>
        void operator()(const State &x, Matrix &J, const double &t, State &dfdt){
           J(0, 0) = -(par_GKL + par_GClL + par_GNaL + std::pow(x[2], 4)*par_GK + (x[4]*par_GAHP)/(x[4] + 1) + x[3]*std::pow(x[1], 3)*par_GNa)/par_Cm;
           J(0, 1) = -(3*x[3]*std::pow(x[1], 2)*par_GNa*(x[0] - (666*log(-(7*x[6] - 270)/x[6]))/25))/par_Cm ;
           J(0, 2) = -(4*std::pow(x[2], 3)*par_GK*(x[0] - (666*log(-(x[5]*par_control)/(x[6] - 158)))/25))/par_Cm ;
           J(0, 3) = -(std::pow(x[1], 3)*par_GNa*(x[0] - (666*log(-(7*x[6] - 270)/x[6]))/25))/par_Cm ;
           J(0, 4) = -((par_GAHP/(x[4] + 1) - (x[4]*par_GAHP)/std::pow(x[4] + 1, 2))*(x[0] - (666*log(-(x[5]*par_control)/(x[6] - 158)))/25))/par_Cm ;
           J(0, 5) = (666*(par_GKL + std::pow(x[2], 4)*par_GK + (x[4]*par_GAHP)/(x[4] + 1)))/(25*x[5]*par_Cm) ;
           J(0, 6) = -((666*(par_GKL + std::pow(x[3], 4)*par_GK + (x[4]*par_GAHP)/(x[4] + 1)))/(25*(x[6] - 158)) + (666*x[6]*par_GNaL*((7*x[6] - 270)/std::pow(x[7], 2) - 7/x[6]))/(25*(7*x[6] - 270)) + (666*x[3]*std::pow(x[1], 3)*x[6]*par_GNa*((7*x[6] - 270)/x[6]^2 - 7/x[6]))/(25*(7*x[7] - 270)))/par_Cm;

           J(1, 0) = (3*x[1] + (3*(x[0]/10 + 3))/((4*exp(- x[0]/18 - 55/18) - (x[0]/10 + 3)/(exp(- x[0]/10 - 3) - 1))*(exp(- x[0]/10 - 3) - 1)))*((2*exp(- x[1]/18 - 55/18))/9 + 1/(10*(exp(- x[1]/10 - 3) - 1)) + (exp(- x[0]/10 - 3)*(x[0]/10 + 3))/(10*(exp(- x[0]/10 - 3) - 1)^2)) - (4*exp(- x[0]/18 - 55/18) - (x[0]/10 + 3)/(exp(- x[0]/10 - 3) - 1))*(3/(10*(4*exp(- x[0]/18 - 55/18) - (x[0]/10 + 3)/(exp(- x[0]/10 - 3) - 1))*(exp(- x[0]/10 - 3) - 1)) + (3*(x[0]/10 + 3)*((2*exp(- x[0]/18 - 55/18))/9 + 1/(10*(exp(- x[0]/10 - 3) - 1)) + (exp(- x[0]/10 - 3)*(x[0]/10 + 3))/(10*std::pow(exp(- x[0]/10 - 3) - 1, 2))))/(std::pow(4*exp(- x[0]/18 - 55/18) - (x[0]/10 + 3)/(exp(- x[0]/10 - 3) - 1), 2)*(exp(- x[0]/10 - 3) - 1)) + (3*exp(- x[0]/10 - 3)*(x[0]/10 + 3))/(10*(4*exp(- x[0]/18 - 55/18) - (x[0]/10 + 3)/(exp(- x[0]/10 - 3) - 1))*std::pow(exp(- x[0]/10 - 3) - 1, 2)));
           J(1, 1) = (3*(x[0]/10 + 3))/(exp(- x[0]/10 - 3) - 1) - 12*exp(- x[0]/18 - 55/18);
           J(1, 2) = 0.0;
           J(1, 3) = 0.0;
           J(1, 4) = 0.0;
           J(1, 5) = 0.0;
           J(1, 6) = 0.0;

           J(2, 0) = (3*x[2] + (3*(x[0]/100 + 17/50))/((exp(- x[0]/80 - 11/20)/8 - (x[0]/100 + 17/50)/(exp(- x[0]/10 - 17/5) - 1))*(exp(- x[0]/10 - 17/5) - 1)))*(exp(- x[0]/80 - 11/20)/640 + 1/(100*(exp(- x[0]/10 - 17/5) - 1)) + (exp(- x[0]/10 - 17/5)*(x[0]/100 + 17/50))/(10*(exp(- x[0]/10 - 17/5) - 1)^2)) - (exp(- x[0]/80 - 11/20)/8 - (x[0]/100 + 17/50)/(exp(- x[0]/10 - 17/5) - 1))*(3/(100*(exp(- x[0]/80 - 11/20)/8 - (x[0]/100 + 17/50)/(exp(- x[0]/10 - 17/5) - 1))*(exp(- x[0]/10 - 17/5) - 1)) + (3*(x[0]/100 + 17/50)*(exp(- x[0]/80 - 11/20)/640 + 1/(100*(exp(- x[0]/10 - 17/5) - 1)) + (exp(- x[0]/10 - 17/5)*(x[0]/100 + 17/50))/(10*(exp(- x[0]/10 - 17/5) - 1)^2)))/((exp(- x[0]/80 - 11/20)/8 - (x[0]/100 + 17/50)/(exp(- x[0]/10 - 17/5) - 1))^2*(exp(- x[0]/10 - 17/5) - 1)) + (3*exp(- x[0]/10 - 17/5)*(x[0]/100 + 17/50))/(10*(exp(- x[0]/80 - 11/20)/8 - (x[0]/100 + 17/50)/(exp(- x[0]/10 - 17/5) - 1))*std::pow(exp(- x[0]/10 - 17/5) - 1, 2)));
           J(2, 1) = 0.0;
           J(2, 2) = (3*(x[0]/100 + 17/50))/(exp(- x[0]/10 - 17/5) - 1) - (3*exp(- x[0]/80 - 11/20))/8
           J(2, 3) = 0.0;
           J(2, 4) = 0.0;
           J(2, 5) = 0.0;
           J(2, 6) = 0.0;

           J(3, 0) = ((7*exp(- x[0]/20 - 11/5))/2000 - exp(- x[0]/10 - 7/5)/(10*(exp(- x[0]/10 - 7/5) + 1)^2))*(3*x[3] - (21*exp(- x[0]/20 - 11/5))/(100*((7*exp(- x[0]/20 - 11/5))/100 + 1/(exp(- x[0]/10 - 7/5) + 1)))) - ((21*exp(- x[0]/20 - 11/5))/(2000*((7*exp(- x[0]/20 - 11/5))/100 + 1/(exp(- x[0]/10 - 7/5) + 1))) - (21*exp(- x[0]/20 - 11/5)*((7*exp(- x[0]/20 - 11/5))/2000 - exp(- x[0]/10 - 7/5)/(10*(exp(- x[0]/10 - 7/5) + 1)^2)))/(100*((7*exp(- x[0]/20 - 11/5))/100 + 1/(exp(- x[0]/10 - 7/5) + 1))^2))*((7*exp(- x[0]/20 - 11/5))/100 + 1/(exp(- x[0]/10 - 7/5) + 1));
           J(3, 1) = 0.0;
           J(3, 2) = 0.0;
           J(3, 3) = -(21*exp(- x[0]/20 - 11/5))/100 - 3/(exp(- x[0]/10 - 7/5) + 1);
           J(3, 4) = 0.0;
           J(3, 5) = 0.0;
           J(3, 6) = 0.0;

           J(4, 0) = -par_GCa/(500*(exp(- (2*x[0])/5 - 10) + 1)) - (par_GCa*exp(- (2*x[0])/5 - 10)*(x[0] - 120))/(1250*(exp(- (2*x[0])/5 - 10) + 1)^2) ;
           J(4, 1) = 0.0;
           J(4, 2) = 0.0;
           J(4, 3) = -1.0/80.0;
           J(4, 4) = 0.0;
           J(4, 5) = 0.0;
           J(4, 6) = 0.0;
        }
};


class push_back_state_and_time {
 ;   public:
        push_back_state_and_time(std::vector<state_type> &states , std::vector<double> &times)
            : m_states(states) , m_times(times) { }

        void operator()(const state_type &x , double t) {
            m_states.push_back(x);
            m_times.push_back(t);
        }

    private:
        std::vector<state_type>& m_states;
        std::vector<double>& m_times;
};


void take1() {
    using namespace std;
    using namespace boost::numeric::odeint;

    // State initialisation
    state_type x(7);
    x[0] = -50;
    x[1] = 0.0936;
    x[2] = 0.96859;
    x[3] = 0.08553;
    x[4] = 0.0;
    x[5] = 7.8;
    x[6] = 15.5;

    // integrate_observ
    vector<state_type> x_vec;
    vector<double> times;

    // RHS class
    Cressman rhs(1, 100, 40, 0.01, 0.05, 0.0175, 0.05, 0.1, 66, 8, 0.0445, 1000, 1);
    const double dt = 1e-6;
    const double T = 100000.0;
    size_t steps = integrate(rhs, x, 0.0, T, dt, push_back_state_and_time(x_vec, times));

    ofstream myfile;
    myfile.open("v_solution.txt");

    for(size_t i = 0; i <= steps; ++i) {
        myfile << times[i] << ", " << x_vec[i][0] << "\n";
    }
    myfile.close();
    std::cout << "Success!" << std::endl;
}


int main(int argc, char** argv) {
    take1();
}
