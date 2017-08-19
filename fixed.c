
#include <math.h>     // for fabs
#include <stdlib.h>   // for strtod
#include <stdio.h>    // for fprintf, printf, stderr

typedef struct {
  double x;
  double p;
  double y;
  double q;
  double z;
  double r;
  double a;
  double b;
} Model;

double fa(Model* m) {
  return (m->r * m->z) / (m->r + m->p)
    + ((m->p) / (m->r + m->p)) * (m->x * m->b) /
      (m->x * m->b + (1 - m->x) * (1 - m->b));
} 

double fb(Model* m) {
  return (m->q * m->y) / (m->q + m->p)
    + ((m->p) / (m->q + m->p)) * ((1 - m->x) * m->a) /
      ((1 - m->x) * m->a + m->x * (1 - m->a));
} 

int main(int argc, char* argv[])
{
  if (argc != 9) {
    fprintf(stderr, "fixed x p y q z r a b\n");
    return 1;
  }

  Model m;
  m.x = strtod(argv[1], NULL);
  m.p = strtod(argv[2], NULL);
  m.y = strtod(argv[3], NULL);
  m.q = strtod(argv[4], NULL);
  m.z = strtod(argv[5], NULL);
  m.r = strtod(argv[6], NULL);
  m.a = strtod(argv[7], NULL);
  m.b = strtod(argv[8], NULL);

  double a;
  double b;
  do {
    a = m.a;
    b = m.b;
    printf("%1.5f  %1.5f\n", a, b);
    m.a = fa(&m);
    m.b = fb(&m);
  } while (fabs(m.a - a) + fabs(m.b - b) > 0.000001);
  return 0;
}

