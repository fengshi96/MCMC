#include <iostream>
#include "src/CIsing.h"
#include "src/visual/Visual.h"
#include "src/visual/Mesh.h"
#include "src/visual/Display.h"

int main() {
    CIsing<float> Ising(130,130);
    Ising.update_J(1);
    Ising.update_t(1);
    Ising.printl();
    Ising.stream(960,960,8000);



}
