#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

// Timing
#include <ctime>

// BOOST
#include <boost/numeric/odeint.hpp>

#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>


typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;



class Cressman
{
    public:
        Cressman(const double Cm, const double GNa, const double GK, const double GAHP,
            const double GKL, const double GNaL, const double GClL, const double GCa, const double Gglia,
            const double Koinf, const double gamma1, const double tau, const double control)
            : Cm(Cm), GNa(GNa), GK(GK), GAHP(GAHP), GKL(GKL), GNaL(GNaL),
              GClL(GClL), GCa(GCa), Gglia(Gglia), Koinf(Koinf),
              gamma1(gamma1), tau(tau), control(control) { }

        template< class State >
        void operator() (const State &x, State &dxdt, const double /* t */ ) {
            const double Ipump = rho*(1.0/(1.0 + exp((25.0 - x[6])/3.0)))*(1.0/(1.0 + exp(5.5 - x[5])));
            const double IGlia = Gglia/(1.0 + exp((18.0 - x[5])/2.5));
            const double Idiff = eps0*(x[5] - Koinf);

            const double Ki = 140.0 + (18.0 - x[6]);
            const double Nao = 144.0 - beta0*(x[6] - 18.0);
            const double ENa = 26.64*log(Nao/x[6]);
            const double EK = 26.64*log(control*x[5]/Ki);

            const double a_m = (3.0 + 0.1*x[0])/(1.0 - exp(-3.0 - 0.1*x[0]));
            const double b_m = 4*exp(-55.0/18 - 1.0/18.0*x[0]);
            const double ah = 0.07*exp(-2.2 - 0.05*x[0]);
            const double bh = 1.0/(1.0 + exp(-1.4 - 0.1*x[0]));
            const double an = (0.34 + 0.01*x[0])/(1.0 - exp(-3.4 - 0-1*x[0]));
            const double bn = 0.125*exp(-0.55 - 0.0125*x[0]);

            const double minf = a_m/(a_m + b_m);
            const double taum = 1.0/(a_m + b_m);
            const double h_inf = ah/(bh + ah);
            const double tauh = 1.0/(bh + ah);
            const double ninf = an/(an + bn);
            const double taun = 1.0/(an + bn);

            const double INa = GNa*pow(x[1], 3)*x[2]*(x[0] - ENa) + GNaL*(x[0] - ENa);
            const double IK = (GK*pow(x[3], 4) + GAHP*x[4]/(1.0 + x[4]) + GKL)*(x[0] - EK);
            const double ICl = GClL*(x[0] - ECl);

            dxdt[0] = -(INa + IK + ICl)/Cm;
            dxdt[1] = phi*(minf - x[1])/taum;
            dxdt[2] = phi*(h_inf - x[2])/tauh;
            dxdt[3] = phi*(ninf - x[3])/taun;
            dxdt[4] = -x[4]/80.0 - 0.002*GCa*(x[0] - ECa)/(1.0 + exp(-(x[0] + 25.0)/2.5));
            dxdt[5] = (gamma1*beta0*IK - 2*beta0*Ipump - IGlia - Idiff)/tau;
            dxdt[6] = -(gamma1*INa + 3*Ipump)/tau;
        }

    private:
        const double rho = 1.25;
        const double eps0 = 1.2;
        const double beta0 = 7.0;
        const double ECa = 120.0;
        const double Cli = 6.0;
        const double Clo = 130.0;
        const double ECl = 26.64*log(Cli/Clo);
        const double phi = 3.0;

        const double Cm;
        const double GNa;
        const double GK;
        const double GAHP;
        const double GKL;
        const double GNaL;
        const double GClL;
        const double GCa;
        const double Gglia;
        const double Koinf;
        const double gamma1;
        const double tau;
        const double control;
};


class Cressman_jacobi
{
    public:
        Cressman_jacobi(const double Cm, const double GNa, const double GK, const double GAHP,
            const double GKL, const double GNaL, const double GClL, const double GCa, const double Gglia,
            const double Koinf, const double gamma1, const double tau, const double control)
            : Cm(Cm), GNa(GNa), GK(GK), GAHP(GAHP), GKL(GKL), GNaL(GNaL),
              GClL(GClL), GCa(GCa), Gglia(Gglia), Koinf(Koinf),
              gamma1(gamma1), tau(tau), control(control) { }

        void operator()(const vector_type &x, matrix_type &J, const double /* &t */, vector_type &dfdt)
        {
            J(0, 0) = (-GAHP*x[4]/(x[4] + 1.0) - GClL - GK*pow(x[3], 4) - GKL - GNa*pow(x[1], 3)*x[2] - GNaL)/Cm;
            J(0, 1) = -3*GNa*pow(x[1], 2)*x[2]*(-26.64*log((-beta0*(x[6] - 18.0) + 144.0)/x[6]) + x[0])/Cm;
            J(0, 2) = -GNa*pow(x[1], 3)*(-26.64*log((-beta0*(x[6] - 18.0) + 144.0)/x[6]) + x[0])/Cm;
            J(0, 3) = 4*GK*pow(x[3], 3)*(26.64*log(control*x[5]/(-x[6] + 158.0)) - x[0])/Cm;
            J(0, 4) = (26.64*log(control*x[5]/(-x[6] + 158.0)) - x[0])*(-GAHP*x[4]/pow(x[4] + 1.0, 2) + GAHP/(x[4] + 1.0))/Cm;
            J(0, 5) = -26.64*(-GAHP*x[4]/(x[4] + 1.0) - GK*pow(x[3], 4) - GKL)/(Cm*x[5]);
            J(0, 6) = (26.64*GNa*pow(x[1], 3)*x[2]*x[6]*(-beta0/x[6] - (-beta0*(x[6] - 18.0) + 144.0)/pow(x[6], 2))/(-beta0*(x[6] - 18.0) + 144.0) + 26.64*GNaL*x[6]*(-beta0/x[6] - (-beta0*(x[6] - 18.0) + 144.0)/pow(x[6], 2))/(-beta0*(x[6] - 18.0) + 144.0) - 26.64*(-GAHP*x[4]/(x[4] + 1.0) - GK*pow(x[3], 4) - GKL)/(-x[6] + 158.0))/Cm;

            J(1, 0) = phi*(-x[1] + (0.1*x[0] + 3.0)/(((0.1*x[0] + 3.0)/(-exp(-0.1*x[0] - 3.0) + 1.0) + 4*exp(-1.0/18.0*x[0] - 55.0/18.0))*(-exp(-0.1*x[0] - 3.0) + 1.0)))*(-1.0*0.1*(0.1*x[0] + 3.0)*exp(-0.1*x[0] - 3.0)/pow(-exp(-0.1*x[0] - 3.0) + 1.0, 2) + 1.0*0.1/(-exp(-0.1*x[0] - 3.0) + 1.0) - 4.0*1.0/18.0*exp(-1.0/18.0*x[0] - 55.0/18.0)) + phi*(1.0*(0.1*x[0] + 3.0)/(-exp(-0.1*x[0] - 3.0) + 1.0) + 4.0*exp(-1.0/18.0*x[0] - 55.0/18.0))*(-0.1*(0.1*x[0] + 3.0)*exp(-0.1*x[0] - 3.0)/(((0.1*x[0] + 3.0)/(-exp(-0.1*x[0] - 3.0) + 1.0) + 4*exp(-1.0/18.0*x[0] - 55.0/18.0))*pow(-exp(-0.1*x[0] - 3.0) + 1.0, 2)) + 0.1/(((0.1*x[0] + 3.0)/(-exp(-0.1*x[0] - 3.0) + 1.0) + 4*exp(-1.0/18.0*x[0] - 55.0/18.0))*(-exp(-0.1*x[0] - 3.0) + 1.0)) + (0.1*x[0] + 3.0)*(0.1*(0.1*x[0] + 3.0)*exp(-0.1*x[0] - 3.0)/pow(-exp(-0.1*x[0] - 3.0) + 1.0, 2) - 0.1/(-exp(-0.1*x[0] - 3.0) + 1.0) + 4*1.0/18.0*exp(-1.0/18.0*x[0] - 55.0/18.0))/(pow((0.1*x[0] + 3.0)/(-exp(-0.1*x[0] - 3.0) + 1.0) + 4*exp(-1.0/18.0*x[0] - 55.0/18.0), 2)*(-exp(-0.1*x[0] - 3.0) + 1.0)));
            J(1, 1) = -phi*(1.0*(0.1*x[0] + 3.0)/(-exp(-0.1*x[0] - 3.0) + 1.0) + 4.0*exp(-1.0/18.0*x[0] - 55.0/18.0));
            J(1, 3) = 0.0;
            J(1, 4) = 0.0;
            J(1, 5) = 0.0;
            J(1, 6) = 0.0;

            J(2, 0) = phi*(1.0*0.07*exp(-0.05*x[0] - 2.2) + 1.0/(exp(-0.1*x[0] - 1.4) + 1.0))*(-0.05*0.07*exp(-0.05*x[0] - 2.2)/(0.07*exp(-0.05*x[0] - 2.2) + 1.0/(exp(-0.1*x[0] - 1.4) + 1.0)) + 0.07*(0.05*0.07*exp(-0.05*x[0] - 2.2) - 1.0*0.1*exp(-0.1*x[0] - 1.4)/pow(exp(-0.1*x[0] - 1.4) + 1.0, 2))*exp(-0.05*x[0] - 2.2)/pow(0.07*exp(-0.05*x[0] - 2.2) + 1.0/(exp(-0.1*x[0] - 1.4) + 1.0), 2)) + phi*(-1.0*0.05*0.07*exp(-0.05*x[0] - 2.2) + 1.0*0.1*exp(-0.1*x[0] - 1.4)/pow(exp(-0.1*x[0] - 1.4) + 1.0, 2))*(0.07*exp(-0.05*x[0] - 2.2)/(0.07*exp(-0.05*x[0] - 2.2) + 1.0/(exp(-0.1*x[0] - 1.4) + 1.0)) - x[2]);
            J(2, 1) = 0.0;
            J(2, 2) = -phi*(1.0*0.07*exp(-0.05*x[0] - 2.2) + 1.0/(exp(-0.1*x[0] - 1.4) + 1.0));
            J(2, 3) = 0.0;
            J(2, 4) = 0.0;
            J(2, 5) = 0.0;
            J(2, 6) = 0.0;

            J(3, 0) = phi*(-x[3] + (0.01*x[0] + 0.34)/((0.125*exp(-0.0125*x[0] - 0.55) + (0.01*x[0] + 0.34)/(-exp(-0.1*x[0] - 3.4) + 1.0))*(-exp(-0.1*x[0] - 3.4) + 1.0)))*(1.0*0.01/(-exp(-0.1*x[0] - 3.4) + 1.0) - 1.0*0.0125*0.125*exp(-0.0125*x[0] - 0.55) - 1.0*0.1*(0.01*x[0] + 0.34)*exp(-0.1*x[0] - 3.4)/pow(-exp(-0.1*x[0] - 3.4) + 1.0, 2)) + phi*(1.0*0.125*exp(-0.0125*x[0] - 0.55) + 1.0*(0.01*x[0] + 0.34)/(-exp(-0.1*x[0] - 3.4) + 1.0))*(0.01/((0.125*exp(-0.0125*x[0] - 0.55) + (0.01*x[0] + 0.34)/(-exp(-0.1*x[0] - 3.4) + 1.0))*(-exp(-0.1*x[0] - 3.4) + 1.0)) - 0.1*(0.01*x[0] + 0.34)*exp(-0.1*x[0] - 3.4)/((0.125*exp(-0.0125*x[0] - 0.55) + (0.01*x[0] + 0.34)/(-exp(-0.1*x[0] - 3.4) + 1.0))*pow(-exp(-0.1*x[0] - 3.4) + 1.0, 2)) + (0.01*x[0] + 0.34)*(-0.01/(-exp(-0.1*x[0] - 3.4) + 1.0) + 0.0125*0.125*exp(-0.0125*x[0] - 0.55) + 0.1*(0.01*x[0] + 0.34)*exp(-0.1*x[0] - 3.4)/pow(-exp(-0.1*x[0] - 3.4) + 1.0, 2))/(pow(0.125*exp(-0.0125*x[0] - 0.55) + (0.01*x[0] + 0.34)/(-exp(-0.1*x[0] - 3.4) + 1.0), 2)*(-exp(-0.1*x[0] - 3.4) + 1.0)));
            J(3, 1) = 0.0;
            J(3, 2) = 0.0;
            J(3, 3) = -phi*(1.0*0.125*exp(-0.0125*x[0] - 0.55) + 1.0*(0.01*x[0] + 0.34)/(-exp(-0.1*x[0] - 3.4) + 1.0));
            J(3, 4) = 0.0;
            J(3, 5) = 0.0;
            J(3, 6) = 0.0;

            J(4, 0) = -0.002*GCa/(exp(10 - x[0]/2.5) + 1.0) - 0.002*GCa*(-ECa + x[0])*exp(10 - x[0]/2.5)/(2.5*pow(exp(10 - x[0]/2.5) + 1.0, 2));
            J(4, 1) = 0.0;
            J(4, 2) = 0.0;
            J(4, 3) = 0.0;
            J(4, 4) = -0.0125;
            J(4, 5) = 0.0;
            J(4, 6) = 0.0;

            J(5, 0) = beta0*gamma1*(GAHP*x[4]/(x[4] + 1.0) + GK*pow(x[3], 4) + GKL)/tau;
            J(5, 1) = 0.0;
            J(5, 2) = 0.0;
            J(5, 3) = 4*GK*beta0*gamma1*pow(x[3], 3)*(-26.64*log(control*x[5]/(-x[6] + 158.0)) + x[0])/tau;
            J(5, 4) = beta0*gamma1*(-26.64*log(control*x[5]/(-x[6] + 158.0)) + x[0])*(-GAHP*x[4]/pow(x[4] + 1.0, 2) + GAHP/(x[4] + 1.0))/tau;
            J(5, 5) = (-1.0/2.5*Gglia*exp(-1.0/2.5*x[5] + 7.2)/pow(exp(-1.0/2.5*x[5] + 7.2) + 1.0, 2) - 26.64*beta0*gamma1*(GAHP*x[4]/(x[4] + 1.0) + GK*pow(x[3], 4) + GKL)/x[5] - 2.0*beta0*rho*exp(5.5 - x[5])/(pow(exp(5.5 - x[5]) + 1.0, 2)*(exp(-1.0/3.0*x[6] + 25.0/3.0) + 1.0)) - eps0)/tau;
            J(5, 6) = (-2.0*1.0/3.0*beta0*rho*exp(-1.0/3.0*x[6] + 25.0/3.0)/((exp(5.5 - x[5]) + 1.0)*pow(exp(-1.0/3.0*x[6] + 25.0/3.0) + 1.0, 2)) - 26.64*beta0*gamma1*(GAHP*x[4]/(x[4] + 1.0) + GK*pow(x[3], 4) + GKL)/(-x[6] + 158.0))/tau;

            J(6, 0) = -gamma1*(GNa*pow(x[1], 3)*x[2] + GNaL)/tau;
            J(6, 1) = -3*GNa*gamma1*pow(x[1], 2)*x[2]*(-26.64*log((-beta0*(x[6] - 18.0) + 144.0)/x[6]) + x[0])/tau;
            J(6, 2) = -GNa*gamma1*pow(x[1], 3)*(-26.64*log((-beta0*(x[6] - 18.0) + 144.0)/x[6]) + x[0])/tau;
            J(6, 3) = 0.0;
            J(6, 4) = 0.0;
            J(6, 5) = -3.0*rho*exp(5.5 - x[5])/(tau*pow(exp(5.5 - x[5]) + 1.0, 2)*(exp(-1.0/3.0*x[6] + 25.0/3.0) + 1.0));
            J(6, 6) = (-3.0*1.0/3.0*rho*exp(-1.0/3.0*x[6] + 25.0/3.0)/((exp(5.5 - x[5]) + 1.0)*pow(exp(-1.0/3.0*x[6] + 25.0/3.0) + 1.0, 2)) - gamma1*(-26.64*GNa*pow(x[1], 3)*x[2]*x[6]*(-beta0/x[6] - (-beta0*(x[6] - 18.0) + 144.0)/pow(x[6], 2))/(-beta0*(x[6] - 18.0) + 144.0) - 26.64*GNaL*x[6]*(-beta0/x[6] - (-beta0*(x[6] - 18.0) + 144.0)/pow(x[6], 2))/(-beta0*(x[6] - 18.0) + 144.0)))/tau;

            dfdt[0] = 0.0;
            dfdt[1] = 0.0;
            dfdt[2] = 0.0;
            dfdt[3] = 0.0;
            dfdt[4] = 0.0;
            dfdt[5] = 0.0;
            dfdt[6] = 0.0;
        }

    private:
        const double rho = 1.25;
        const double eps0 = 1.2;
        const double beta0 = 7.0;
        const double ECa = 120.0;
        const double Cli = 6.0;
        const double Clo = 130.0;
        const double ECl = 26.64*log(Cli/Clo);
        const double phi = 3.0;

        const double Cm;
        const double GNa;
        const double GK;
        const double GAHP;
        const double GKL;
        const double GNaL;
        const double GClL;
        const double GCa;
        const double Gglia;
        const double Koinf;
        const double gamma1;
        const double tau;
        const double control;
};


template< class State_type >
class push_back_state_and_time
{
    public:
        push_back_state_and_time(std::vector< State_type > &states , std::vector< double > &times)
            : m_states(states), m_times(times) { }

        void operator()(const State_type &x , double t)
        {
            m_states.push_back(x);
            m_times.push_back(t);
        }

    private:
        std::vector< State_type > &m_states;
        std::vector< double > &m_times;
};


void take1()
{
    using namespace std;
    using namespace boost::numeric::odeint;

    namespace phoenix = boost::phoenix;

    // State initialisation
    vector_type x(7);
    x[0] = -50;
    x[1] = 0.0936;
    x[2] = 0.96859;
    x[3] = 0.08553;
    x[4] = 0.0;
    x[5] = 7.8;
    x[6] = 15.5;

    // integrate_observ
    vector< vector_type > x_vec;
    vector<double> times;

    // RHS class
    Cressman rhs(1, 100, 40, 0.01, 0.05, 0.0175, 0.05, 0.1, 66, 8, 0.0445, 1000, 1);
    Cressman_jacobi jacobi(1, 100, 40, 0.01, 0.05, 0.0175, 0.05, 0.1, 66, 8, 0.0445, 1000, 1);

    const double dt = 1e-1;
    const double T = 1e1;

    // Main solver call
    auto tick = clock();
    size_t num_of_steps = integrate_const(make_dense_output< rosenbrock4< double > >(1e-6, 1e-6),
            make_pair(rhs, jacobi), x, 0.0, T, dt,
            push_back_state_and_time< vector_type >(x_vec, times));
    auto tock = clock();
    cout << num_of_steps << " " << tock - tick << endl;


    auto x2 = x;
    tick = clock();
    size_t num_of_steps2 = integrate_const(make_dense_output< runge_kutta_dopri5< vector_type > >(1.0e-6, 1.0e-6),
            rhs, x2, 0.0, T, dt,
            push_back_state_and_time< vector_type >(x_vec, times));
    tock = clock();
    cout << num_of_steps2 << " " << tock - tick << endl;


    assert(times.size() == x_vec.size());
    ofstream myfile;
    myfile.open("v_solution.txt");

    for (size_t i = 0; i < times.size(); ++i)
    {
        myfile << times[i] << ", " << x_vec[i][0] << "\n";
    }
    myfile.close();
    std::cout << "Success!" << std::endl;
}


int main(int argc, char** argv)
{
    take1();
}
