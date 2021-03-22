#include <iostream>
#include "src/CIsing.h"
#include "src/visual/Visual.h"
#include "src/visual/Mesh.h"
#include "src/visual/Display.h"

int main() {
    CIsing<float> Ising(200,200);
    Ising.update_J(-1);
    Ising.update_t(1);
    Ising.printl();
    Ising.stream(600,600,4000);

    std::cout << .3*8 << std::endl;

}
