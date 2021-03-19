#include <iostream>
#include "src/CIsing.h"
#include "src/visual/Visual.h"
#include "src/visual/Mesh.h"

int main() {
    CIsing<float> Ising(300,300);
    Ising.update_J(-1);
    Ising.update_t(1);
    Ising.printl();
    Ising.stream(800,800,10000);


}
