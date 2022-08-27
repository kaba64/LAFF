#include <stdio.h>
#include <cblas-openblas.h>

int main(int argc, char **argv) {
  double x[3] = {1.0, 2.0, 3.0};
  double alpha ;

  printf("Pre dscal\t: %g %g %g\n", x[0], x[1], x[2]);

  alpha = cblas_ddot(3, x, 1, x, 1);

  printf("%lf\n",alpha);
  //cblas_dscal(3, alpha, x, 1);

  //printf("Post dscal\t: %g %g %g\n", x[0], x[1], x[2]);
  return 0;
}
