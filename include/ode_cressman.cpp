#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

// Timing
#include <ctime>

// ODEINT
#include <boost/numeric/odeint.hpp>

#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>

typedef double float_type;
typedef boost::numeric::ublas::vector< float_type > vector_type;
typedef boost::numeric::ublas::matrix< float_type > matrix_type;
/* typedef std::vector< double > vector_type; */



template< typename FT >
class Cressman
{
    public:
        Cressman(const FT Cm, const FT GNa, const FT GK, const FT GAHP,
            const FT GKL, const FT GNaL, const FT GClL, const FT GCa, const FT Gglia,
            const FT Koinf, const FT gamma1, const FT tau, const FT control)
            : Cm(Cm), GNa(GNa), GK(GK), GAHP(GAHP), GKL(GKL), GNaL(GNaL),
              GClL(GClL), GCa(GCa), Gglia(Gglia), Koinf(Koinf),
              gamma1(gamma1), tau(tau), control(control) { }

        template< class State >
        void operator() (const State &x, State &dxdt, const FT /* t */ ) {
            using namespace std;
            const FT Ipump = rho*(1/(1 + exp((25 - x[6])/3.)))*(1/(1 + exp(5.5- x[5])));
            const FT IGlia = Gglia/(1 + exp((18.- x[5])/2.5));
            const FT Idiff = eps0*(x[5] - Koinf);

            const FT Ki = 140 + (18 - x[6]);
            const FT Nao = 144 - beta0*(x[6] - 18);
            const FT ENa = 26.64*log(Nao/x[6]);
            const FT EK = 26.64*log(control*x[5]/Ki);
            const FT ECl = 26.64*log(Cli/Clo);

            const FT am = (0.1*x[0] + 3.0)/(-exp(-0.1*x[0] - 3) + 1);
            const FT bm = 4*exp(-1./18.*x[0] - 55./18.);
            const FT ah = 0.07*exp(-0.05*x[0] - 2.2);
            const FT bh = 1.0/(exp(-0.1*x[0] - 1.4) + 1);
            const FT an = (0.01*x[0] + 0.34)/(-exp(-0.1*x[0] - 3.4) + 1);
            const FT bn = 0.125*exp(-0.0125*x[0] - 0.55);

            const FT taum = 1./(am + bm);
            const FT minf = 1./(am + bm)*am;
            const FT taun = 1./(an + bn);
            const FT ninf = 1./(an + bn)*an;
            const FT tauh = 1./(bh + ah);
            const FT hinf = 1./(bh + ah)*ah;

            const FT INa = GNa*pow(x[1], 3)*x[3]*(x[0] - ENa) + GNaL*(x[0] - ENa);
            const FT IK = (GK*pow(x[2], 4) + GAHP*x[4]/(1 + x[4]) + GKL)*(x[0] - EK);
            const FT ICl = GClL*(x[0] - ECl);

            dxdt[0] = -(INa + IK + ICl)/Cm;
            dxdt[1] = phi*(minf - x[1])/taum;
            dxdt[2] = phi*(ninf - x[2])/taun;
            dxdt[3] = phi*(hinf - x[3])/tauh;
            dxdt[4] = x[4]/80. - 0.002*GCa*(x[0] - ECa)/(1 + exp(-(x[0] + 25.)/2.5));
            dxdt[5] = (gamma1*beta0*IK - 2*beta0*Ipump - IGlia - Idiff)/tau;
            dxdt[6] = -(gamma1*INa + 3*Ipump)/tau;
        }

    private:
        const FT rho = 1.25;
        const FT eps0 = 1.2;
        const FT beta0 = 7.0;
        const FT ECa = 120.0;
        const FT Cli = 6.0;
        const FT Clo = 130.0;
        const FT ECl = 26.64*log(Cli/Clo);
        const FT phi = 3.0;

        const FT Cm;
        const FT GNa;
        const FT GK;
        const FT GAHP;
        const FT GKL;
        const FT GNaL;
        const FT GClL;
        const FT GCa;
        const FT Gglia;
        const FT Koinf;
        const FT gamma1;
        const FT tau;
        const FT control;
};


template< typename FT >
class Cressman_jacobi
{
    public:
        Cressman_jacobi(const FT Cm, const FT GNa, const FT GK, const FT GAHP,
            const FT GKL, const FT GNaL, const FT GClL, const FT GCa, const FT Gglia,
            const FT Koinf, const FT gamma1, const FT tau, const FT control)
            : Cm(Cm), GNa(GNa), GK(GK), GAHP(GAHP), GKL(GKL), GNaL(GNaL),
              GClL(GClL), GCa(GCa), Gglia(Gglia), Koinf(Koinf),
              gamma1(gamma1), tau(tau), control(control) { }

        void operator()(const vector_type &x, matrix_type &J, const FT /* &t */, vector_type &dfdt)
        {
            J(0, 0) = (-GAHP*x[4]/(x[4] + 1) - GClL - GK*pow(x[2], 4) - GKL - GNa*pow(x[1], 3)*x[3] - GNaL)/Cm;
            J(0, 1) = -3*GNa*pow(x[1], 2)*x[3]*(-26.64*log((-beta0*(x[6] - 18) + 144)/x[6]) + x[0])/Cm;
            J(0, 2) = 4*GK*pow(x[2], 3)*(26.64*log(control*x[5]/(-x[6] + 158)) - x[0])/Cm;
            J(0, 3) = -GNa*pow(x[1], 3)*(-26.64*log((-beta0*(x[6] - 18) + 144)/x[6]) + x[0])/Cm;
            J(0, 4) = (26.64*log(control*x[5]/(-x[6] + 158)) - x[0])*(-GAHP*x[4]/pow(x[4] + 1, 2) + GAHP/(x[4] + 1))/Cm;
            J(0, 5) = -26.64*(-GAHP*x[4]/(x[4] + 1) - GK*pow(x[2], 4) - GKL)/(Cm*x[5]);
            J(0, 6) = (26.64*GNa*pow(x[1], 3)*x[3]*x[6]*(-beta0/x[6] - (-beta0*(x[6] - 18) + 144)/pow(x[6], 2))/(-beta0*(x[6] - 18) + 144) + 26.64*GNaL*x[6]*(-beta0/x[6] - (-beta0*(x[6] - 18) + 144)/pow(x[6], 2))/(-beta0*(x[6] - 18) + 144) - 26.64*(-GAHP*x[4]/(x[4] + 1) - GK*pow(x[2], 4) - GKL)/(-x[6] + 158))/Cm;



            J(1, 0) = phi*(-x[1] + 1.0*(0.1*x[0] + 3.0)/(((0.1*x[0] + 3.0)/(-exp(-0.1*x[0] - 3) + 1) + 4*exp(-1./18.*x[0] - 55./18.))*(-exp(-0.1*x[0] - 3) + 1)))*(-1.0*0.1*(0.1*x[0] + 3.0)*exp(-0.1*x[0] - 3)/pow(-exp(-0.1*x[0] - 3) + 1, 2) + 1.0*0.1/(-exp(-0.1*x[0] - 3) + 1) - 4.0*1./18.*exp(-1./18.*x[0] - 55./18.)) + phi*(1.0*(0.1*x[0] + 3.0)/(-exp(-0.1*x[0] - 3) + 1) + 4.0*exp(-1./18.*x[0] - 55./18.))*(-1.0*0.1*(0.1*x[0] + 3.0)*exp(-0.1*x[0] - 3)/(((0.1*x[0] + 3.0)/(-exp(-0.1*x[0] - 3) + 1) + 4*exp(-1./18.*x[0] - 55./18.))*pow(-exp(-0.1*x[0] - 3) + 1, 2)) + 1.0*0.1/(((0.1*x[0] + 3.0)/(-exp(-0.1*x[0] - 3) + 1) + 4*exp(-1./18.*x[0] - 55./18.))*(-exp(-0.1*x[0] - 3) + 1)) + 1.0*(0.1*x[0] + 3.0)*(0.1*(0.1*x[0] + 3.0)*exp(-0.1*x[0] - 3)/pow(-exp(-0.1*x[0] - 3) + 1, 2) - 0.1/(-exp(-0.1*x[0] - 3) + 1) + 4*1./18.*exp(-1./18.*x[0] - 55./18.))/(pow((0.1*x[0] + 3.0)/(-exp(-0.1*x[0] - 3) + 1) + 4*exp(-1./18.*x[0] - 55./18.), 2)*(-exp(-0.1*x[0] - 3) + 1)));
            J(1, 1) = -phi*(1.0*(0.1*x[0] + 3.0)/(-exp(-0.1*x[0] - 3) + 1) + 4.0*exp(-1./18.*x[0] - 55./18.));
            J(1, 3) = 0.0;
            J(1, 4) = 0.0;
            J(1, 5) = 0.0;
            J(1, 6) = 0.0;

            J(2, 0) = phi*(-x[2] + 1.0*(0.01*x[0] + 0.34)/((0.125*exp(-0.0125*x[0] - 0.55) + (0.01*x[0] + 0.34)/(-exp(-0.1*x[0] - 3.4) + 1))*(-exp(-0.1*x[0] - 3.4) + 1)))*(1.0*0.01/(-exp(-0.1*x[0] - 3.4) + 1) - 1.0*0.0125*0.125*exp(-0.0125*x[0] - 0.55) - 1.0*0.1*(0.01*x[0] + 0.34)*exp(-0.1*x[0] - 3.4)/pow(-exp(-0.1*x[0] - 3.4) + 1, 2)) + phi*(1.0*0.125*exp(-0.0125*x[0] - 0.55) + 1.0*(0.01*x[0] + 0.34)/(-exp(-0.1*x[0] - 3.4) + 1))*(1.0*0.01/((0.125*exp(-0.0125*x[0] - 0.55) + (0.01*x[0] + 0.34)/(-exp(-0.1*x[0] - 3.4) + 1))*(-exp(-0.1*x[0] - 3.4) + 1)) - 1.0*0.1*(0.01*x[0] + 0.34)*exp(-0.1*x[0] - 3.4)/((0.125*exp(-0.0125*x[0] - 0.55) + (0.01*x[0] + 0.34)/(-exp(-0.1*x[0] - 3.4) + 1))*pow(-exp(-0.1*x[0] - 3.4) + 1, 2)) + 1.0*(0.01*x[0] + 0.34)*(-0.01/(-exp(-0.1*x[0] - 3.4) + 1) + 0.0125*0.125*exp(-0.0125*x[0] - 0.55) + 0.1*(0.01*x[0] + 0.34)*exp(-0.1*x[0] - 3.4)/pow(-exp(-0.1*x[0] - 3.4) + 1, 2))/(pow(0.125*exp(-0.0125*x[0] - 0.55) + (0.01*x[0] + 0.34)/(-exp(-0.1*x[0] - 3.4) + 1), 2)*(-exp(-0.1*x[0] - 3.4) + 1)));
            J(2, 1) = 0.0;
            J(2, 2) = -phi*(1.0*0.125*exp(-0.0125*x[0] - 0.55) + 1.0*(0.01*x[0] + 0.34)/(-exp(-0.1*x[0] - 3.4) + 1));
            J(2, 3) = 0.0;
            J(2, 4) = 0.0;
            J(2, 5) = 0.0;
            J(2, 6) = 0.0;

            J(3, 1) = phi*(1.0*0.07*exp(-0.05*x[0] - 2.2) + 1.0/(exp(-0.1*x[0] - 1.4) + 1))*(-1.0*0.05*0.07*exp(-0.05*x[0] - 2.2)/(0.07*exp(-0.05*x[0] - 2.2) + 1.0/(exp(-0.1*x[0] - 1.4) + 1)) + 1.0*0.07*(0.05*0.07*exp(-0.05*x[0] - 2.2) - 0.1*exp(-0.1*x[0] - 1.4)/pow(exp(-0.1*x[0] - 1.4) + 1, 2))*exp(-0.05*x[0] - 2.2)/pow(0.07*exp(-0.05*x[0] - 2.2) + 1.0/(exp(-0.1*x[0] - 1.4) + 1), 2)) + phi*(-1.0*0.05*0.07*exp(-0.05*x[0] - 2.2) + 1.0*0.1*exp(-0.1*x[0] - 1.4)/pow(exp(-0.1*x[0] - 1.4) + 1, 2))*(1.0*0.07*exp(-0.05*x[0] - 2.2)/(0.07*exp(-0.05*x[0] - 2.2) + 1.0/(exp(-0.1*x[0] - 1.4) + 1)) - x[3]);
            J(3, 0) = 0.0 ;
            J(3, 2) = 0.0;
            J(3, 3) = -phi*(1.0*0.07*exp(-0.05*x[0] - 2.2) + 1.0/(exp(-0.1*x[0] - 1.4) + 1));
            J(3, 4) = 0.0;
            J(3, 5) = 0.0;
            J(3, 6) = 0.0;

            J(4, 0) = -0.002*GCa/(exp((-25. - x[0])/2.5) + 1) - 0.002*GCa*(-ECa + x[0])*exp((-25. - x[0])/2.5)/(2.5*pow(exp((-25. - x[0])/2.5) + 1, 2));
            J(4, 1) = 0.0;
            J(4, 2) = 0.0;
            J(4, 3) = 0.0;
            J(4, 4) = 1.0/80.0;
            J(4, 5) = 0.0;
            J(4, 6) = 0.0;

            J(5, 0) = beta0*gamma1*(GAHP*x[4]/(x[4] + 1) + GK*pow(x[2], 4) + GKL)/tau;
            J(5, 1) = 0.0;
            J(5, 2) = 4*GK*beta0*gamma1*pow(x[2], 3)*(-26.64*log(control*x[5]/(-x[6] + 158)) + x[0])/tau;
            J(5, 3) = 0.0;
            J(5, 4) = beta0*gamma1*(-26.64*log(control*x[5]/(-x[6] + 158)) + x[0])*(-GAHP*x[4]/pow(x[4] + 1, 2) + GAHP/(x[4] + 1))/tau;
            J(5, 5) = (-26.64*beta0*gamma1*(GAHP*x[4]/(x[4] + 1) + GK*pow(x[2], 4) + GKL)/x[5] - 2*beta0*rho*exp(5.5 - x[5])/(pow(exp(5.5 - x[5]) + 1, 2)*(exp(25./3. - x[6]/3.) + 1)) - eps0 - Gglia*exp(18./2.5 - x[5]/2.5)/(2.5*pow(exp(18./2.5 - x[5]/2.5) + 1, 2)))/tau;
            J(5, 6) = (-26.64*beta0*gamma1*(GAHP*x[4]/(x[4] + 1) + GK*pow(x[2], 4) + GKL)/(-x[6] + 158) - 2*beta0*rho*exp(25./3. - x[6]/3.)/(3.*(exp(5.5 - x[5]) + 1)*pow(exp(25./3. - x[6]/3.) + 1, 2)))/tau;

            J(6, 0) = -gamma1*(GNa*pow(x[1], 3)*x[3] + GNaL)/tau;
            J(6, 1) = -3*GNa*gamma1*pow(x[1], 2)*x[3]*(-26.64*log((-beta0*(x[6] - 18) + 144)/x[6]) + x[0])/tau;
            J(6, 2) = 0.0;
            J(6, 3) = -GNa*gamma1*pow(x[1], 3)*(-26.64*log((-beta0*(x[6] - 18) + 144)/x[6]) + x[0])/tau;
            J(6, 4) = 0.0;
            J(6, 5) = -3*rho*exp(5.5 - x[5])/(tau*pow(exp(5.5 - x[5]) + 1, 2)*(exp(25./3. - x[6]/3.) + 1));
            J(6, 6) = (-gamma1*(-26.64*GNa*pow(x[1], 3)*x[3]*x[6]*(-beta0/x[6] - (-beta0*(x[6] - 18) + 144)/pow(x[6], 2))/(-beta0*(x[6] - 18) + 144) - 26.64*GNaL*x[6]*(-beta0/x[6] - (-beta0*(x[6] - 18) + 144)/pow(x[6], 2))/(-beta0*(x[6] - 18) + 144)) - 3*rho*exp(25./3. - x[6]/3.)/(3.*(exp(5.5 - x[5]) + 1)*pow(exp(25./3. - x[6]/3.) + 1, 2)))/tau;

            dfdt[0] = 0.0;
            dfdt[1] = 0.0;
            dfdt[2] = 0.0;
            dfdt[3] = 0.0;
            dfdt[4] = 0.0;
            dfdt[5] = 0.0;
            dfdt[6] = 0.0;
        }

    private:
        const FT rho = 1.25;
        const FT eps0 = 1.2;
        const FT beta0 = 7.0;
        const FT ECa = 120.0;
        const FT Cli = 6.0;
        const FT Clo = 130.0;
        const FT ECl = 26.64*log(Cli/Clo);
        const FT phi = 3.0;

        const FT Cm;
        const FT GNa;
        const FT GK;
        const FT GAHP;
        const FT GKL;
        const FT GNaL;
        const FT GClL;
        const FT GCa;
        const FT Gglia;
        const FT Koinf;
        const FT gamma1;
        const FT tau;
        const FT control;
};


template< class State_type, typename FT >
class push_back_state_and_time
{
    public:
        push_back_state_and_time(std::vector< State_type > &states , std::vector< double > &times)
            : m_states(states), m_times(times) { }

        void operator()(const State_type &x , FT t)
        {
            m_states.push_back(x);
            m_times.push_back(t);
        }

    private:
        std::vector< State_type > &m_states;
        std::vector< FT > &m_times;
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
    Cressman< float_type > rhs(1., 100., 40., 0.01, 0.05, 0.0175, 0.05, 0.1, 66, 8., 0.0445, 1000, 1.);
    Cressman_jacobi< float_type > jacobi(1., 100., 40., 0.01, 0.05, 0.0175, 0.05, 0.1, 66., 8., 0.0445, 1000, 1.);

    const double dt = 1e-1;
    const double T = 1e3;
    const double abs_tol = 1e-3, rel_tol=1e-2;

    typedef rosenbrock4< float_type > error_stepper_type;       // Should be FT
    typedef rosenbrock4_controller < error_stepper_type  > controller_stepper_type;
    controller_stepper_type controller_stepper(abs_tol, rel_tol);

    /* typedef rosenbrock4< vector_type > error_stepper_type; */
    /* typedef rosenbrock4_controller < error_stepper_type  > controller_stepper_type; */
    /* typedef rosenbrock4_dense_output < controller_stepper_type > dense_stepper_type; */

    /* controller_stepper_type controlled_stepper(abs_tol, rel_tol, error_stepper_type()); */
    /* dense_stepper_type dense_stepper(controlled_stepper); */
    /* dense_stepper.initialize(x, 1.0e-6, 1.0e-6); */

    auto tick = clock();
    size_t num_of_steps = integrate_adaptive(controller_stepper, make_pair(rhs, jacobi),
            x, 0.0, T, dt,
            push_back_state_and_time< vector_type, float_type >(x_vec, times));
    double tock = static_cast< double >(clock() - tick);
    cout << num_of_steps << " " << tock/CLOCKS_PER_SEC << endl;

    /* auto tick = clock(); */
    /* size_t num_of_steps = integrate_adaptive(controlled_stepper, */
    /*         make_pair(rhs, jacobi), x, 0.0, T, dt, */ 
    /*         push_back_state_and_time< vector_type >(x_vec, times)); */

    /* double tock = static_cast< double >(clock() - tick); */
    /* cout << num_of_steps << " " << tock/CLOCKS_PER_SEC << endl; */


    /* typedef runge_kutta_cash_karp54< vector_type > error_stepper_type; */
    /* typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type; */
    /* controlled_stepper_type controlled_stepper(abs_tol, rel_tol, error_stepper_type()); */
    /* auto tick = clock(); */
    /* size_t num_of_steps = integrate_adaptive(controlled_stepper, rhs, x, 0.0, T, dt); */
    /* double tock = static_cast< double >(clock() - tick); */
    /* cout << num_of_steps << " " << tock/CLOCKS_PER_SEC << endl; */


    /* typedef runge_kutta_dopri5< vector_type > error_stepper_type; */
    /* typedef runge_kutta_fehlberg78< vector_type > error_stepper_type; */
    /* typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type; */
    /* typedef dense_output_runge_kutta <controlled_stepper_type > dense_stepper_type; */

    /* typedef bulirsch_stoer< vector_type > controlled_stepper_type; */

    /* controlled_stepper_type controlled_stepper(abs_tol, rel_tol); */
    /* auto tick = clock(); */
    /* size_t num_of_steps = integrate_adaptive(controlled_stepper, rhs, x, 0.0, T, dt, */
    /*         push_back_state_and_time< vector_type >(x_vec, times)); */
    /* double tock = static_cast< double >(clock() - tick); */
    /* cout << num_of_steps << " " << tock/CLOCKS_PER_SEC << endl; */
    /* cout << times.back() << endl; */

    /* vector_type dxdt(7); */
    /* rhs(x, dxdt, 0.0); */

    // ***********************************************************************************
    // Main solver call
    // ***********************************************************************************
    /* auto tick = clock(); */
    /* size_t num_of_steps = integrate_const(make_dense_output< rosenbrock4< float_type > >(1e-6, 1e-6), */
    /*         make_pair(rhs, jacobi), x, 0.0, T, dt, */
    /*         push_back_state_and_time< vector_type >(x_vec, times)); */
    /* auto tock = clock(); */
    /* cout << num_of_steps << " " << tock - tick << endl; */

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
