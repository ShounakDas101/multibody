#include <stdio.h>
#include <math.h>

/*Function representing the 2nd order differential equation (modify this according to your equation)
  My Equation is:
  */
double f(double t, double y, double z) {
// Example: dz/dt = d^2y/dt^2 =  = -A*Sin(2y) +B*Sin(wt)*Sin(y) -Cz + D
double A=1, B=2, C=3, D=4, w=5;
return -A*sin(2*y) + B*sin(w*t)*sin(y) -C*z + D;
}

// Euler method implementation
void eulerMethod(double t0, double y0, double z0, double dt, double t_end) {
double t = t0;
double y = y0;
double z = z0;

while (t <= t_end) {
printf("t = %lf, y = %lf\n", t, y);

double y_next = y + dt * z;
double z_next = z + dt * f(t, y, z);

t += dt;
y = y_next;
z = z_next;
}
}

int main() {
// Initial conditions and step size
double t0 = 0.0;      // Initial value of x
double y0 = 1.0;      // Initial value of y
double z0 = 0.0;      // Initial value of dy/dx (z)
double dt = 0.1;       // Step size
double t_end = 10.0; // Target value of x

eulerMethod(t0, y0, z0, dt, t_end);

return 0;
}