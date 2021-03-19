#include <iostream>
#include "src/CIsing.h"
#include "src/visual/Visual.h"
#include "src/visual/Mesh.h"

int main() {
    CIsing<float> Ising(10,10);
//    Ising.evolve(500);
    Ising.printl();
    Ising.stream(800,800);


}
