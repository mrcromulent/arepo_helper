#include "const.h"

const double msol = 1.989e33; // solar mass
const double G    = 6.6738e-8; // gravitational constant
const double electron_charge_esu = 4.80320427e-10;

#define X(a, b) b,
const char* const f [] =
{
        NAMES
};
#undef X