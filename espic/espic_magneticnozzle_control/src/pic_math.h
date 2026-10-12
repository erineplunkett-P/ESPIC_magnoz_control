#ifndef ESPIC_MATH_H
#define ESPIC_MATH_H

#define SQR(a) ((a) * (a))

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

#define VEC_DOT(a, b) ((a[0]) * (b[0]) + (a[1]) * (b[1]) + (a[2]) * (b[2]));

#endif
