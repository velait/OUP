functions {

  vector integrand(vector x) {
    return(exp((0.5*2*square(x) - 0.25*square(square(x)))/1));
  }

  real[] integrand_ode(real r, real[] f, real[] theta, real[] x_r, int[] x_i) {
    real df_dx[1];
    real x = logit(r);
    df_dx[1] = exp((theta[1]*x + 0.5*theta[2]*square(x) - 0.25*square(square(x)))/theta[4]) * 1/(r * (1-r));
    return(df_dx);
  }

  real ode_integrate(real[] theta) {
    
    int x_i[0];
    return(
      integrate_ode_rk45(
        integrand_ode,
        rep_array(0.0, 1),
        1E-5,
        rep_array(1.0-1E-5, 1),
        theta,
        rep_array(0.0, 0),
        x_i)[1,1]
      );
      
  }

}
data {
}
model {}
