library(BH)
library(RcppArmadillo)
src <- "
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <boost/numeric/odeint.hpp>
#include <vector>
namespace boost { namespace numeric { namespace odeint {
      template <>
      struct is_resizeable<arma::vec>
      {
	typedef boost::true_type type;
	const static bool value = type::value;
      };
      template <>
      struct same_size_impl<arma::vec, arma::vec>
      {
	static bool same_size(const arma::vec& x, const arma::vec& y)
	{
	  return x.size() == y.size();
	}
      };
      template<>
      struct resize_impl<arma::vec, arma::vec>
      {
	static void resize(arma::vec& v1, const arma::vec& v2)
	{
	  v1.resize(v2.size());
	}
      };
    } } } // namespace boost::numeric::odeint
typedef arma::vec state_type;
typedef std::vector<state_type> vector_state_type;
typedef std::vector<double> vector_times;
struct push_back_state_and_time
{
    std::vector< state_type >& m_states;
    std::vector< double >& m_times;
    push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times )
    : m_states( states ) , m_times( times ) { }
    void operator()( const state_type &x , double t )
    {
        m_states.push_back( x );
        m_times.push_back( t );
    }
};
class SEIR {
public:
  double p, epsilon, theta, t_today, N, p_lower_inf, eta, p_symp, gammaD, gamma_pos, t_start, t_end;
  SEIR(double p, double epsilon, double theta, double t_today, double N, double p_lower_inf,
       double eta, double p_symp, double gammaD, double gamma_pos, double t_start, double t_end) : p(p), epsilon(epsilon), theta(theta),
						   t_today(t_today), N(N), p_lower_inf(p_lower_inf),
						   eta(eta), p_symp(p_symp), gammaD(gammaD), gamma_pos(gamma_pos), t_start(t_start), t_end(t_end) { }
  double b(double t) {
    return ((1-p)/(1+exp(-epsilon*(t-t_today))) +p)*theta;
  }
  void operator() ( const state_type &x , state_type &dxdt , const double t )
  {
    double S =  x(0);
    double E =  x(1);
    double I_symp = x(2);
    double I_asymp = x(3);
    double R1 = x(4);
    // double R2 = x(5);
    dxdt(0) = -b(t) * S * I_symp/N - p_lower_inf*b(t) * S * I_asymp/N;
    dxdt(1) = b(t) * S * I_symp/N + p_lower_inf*b(t) * S * I_asymp/N - eta*E;
    dxdt(2) = p_symp * eta * E      - gammaD * I_symp;
    dxdt(3) = (1 - p_symp)* eta * E - gammaD * I_asymp;
    dxdt(4)   = gammaD * (I_symp + I_asymp) - gamma_pos * R1;
    dxdt(5)   =  gamma_pos * R1;
  }
  arma::mat run(arma::vec init, double from, double to, double step) {
    vector_state_type states;
    vector_times times;
    BOOST_STATIC_ASSERT( boost::numeric::odeint::is_resizeable<state_type>::value == true );
    boost::numeric::odeint::runge_kutta4< state_type > stepper;
    size_t steps = boost::numeric::odeint::integrate_const(stepper,
                                                           *this, init,
							   from, to, step,
							   push_back_state_and_time(states, times));
    size_t nx = states[0].size(), nTimes = times.size();
    // combine the results
    arma::mat combined(nTimes,nx+1);
    combined.col(0) = arma::vec(&times[0], times.size());
    for (size_t i=0; i<nTimes; i++)
      combined(arma::span(i,i),arma::span(1,nx)) = states[i].t();
    return(combined);
  }
};
// [[Rcpp::export]]
arma::mat SEIRmodel2(Rcpp::NumericVector state, Rcpp::NumericVector parms) {
  double p = parms[\"p\"];
  double epsilon = parms[\"epsilon\"];
  double theta = parms[\"theta\"];
  double t_today = parms[\"t_today\"];
  double N = parms[\"N\"];
  double p_lower_inf = parms[\"p_lower_inf\"];
  double eta = parms[\"eta\"];
  double p_symp = parms[\"p_symp\"];
  double gammaD = parms[\"gammaD\"];
  double gamma_pos = parms[\"gamma_pos\"];
  double t_start = parms[\"t_start\"];
  double t_end = parms[\"t_end\"];
  arma::vec init(&state[0], 6);
  SEIR model(p, epsilon, theta, t_today, N, p_lower_inf, eta, p_symp, gammaD, gamma_pos, t_start, t_end);
  return model.run(init, t_start, t_end, 1.0); 
}"
  cat(src,file="SEIRmodel2.cpp")
  sourceCpp("SEIRmodel2.cpp")
  