#include <iostream>
#include "src/CIsing.h"

int main() {
    CIsing<float> Ising(5,5);
    // Ising.update_J(1);
    int Niter = 500;
    Ising.evolve(Niter);

    Ising.printl();
    Ising.prints();
}
