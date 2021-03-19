//
// Created by shifeng on 3/13/21.
//

#ifndef MCMC_CISING_H
#define MCMC_CISING_H
#include "Matrix.h"
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include "visual/Visual.h"

template<class T>
class CIsing {
public:
    CIsing() : J(-1), t(1), e(0), f(0), m(0), h(2), w(2), Lat(2,2), spin(1) {
        config_init();
    }

    CIsing(int size1, int size2) : J(-1), t(1), e(0), f(0), m(0), h(size1), w(size2), Lat(size1,size2), spin(1) {
        config_init();
    }

    inline T temp() {return t;}
    inline T spin_type() {return spin;}
    inline int height() {return h;}
    inline int width() {return w;}
    inline int size() {
        std::cout << "(" << h << ", " << w << ") Square Lattice" << std::endl;
        return h * w;
    }
    void printl() { Lat.print();}  // print lattice
    void prints() {
        // print attribute status
        std::cout << "Lattice shape = (" << w << "," << h << ")" << std::endl;
        std::cout << "Temperature = " << t << std::endl;
        std::cout << "J = " << J << std::endl;
        std::cout << "Field = " << f << std::endl;
        std::cout << "Energy = " << e << std::endl;
        std::cout << "Magnitization = " << m << std::endl;
    }

    void config_init();

    T& operator()(int i, int j) {
        // assert(i < h && j < w);  // OBC
        if(i<0) i += h;      // PBC
        if(i>=h) i -= h;
        if(j<0) j += w;
        if(j>=w) j -= w;
        return Lat(i, j);
    }  // returns a reference of that element

    bool tryflip();  // flip or not (for M-H method)
    bool tryclusterflip();
    void evolve(const int nstep);  // number of iterations

    T energy();  // calculate energy
    T mag() {return m;}  // magnitization

    void update_t(T temp) {t = temp;}
    void update_e(T energy) {e = energy;}
    void update_m(T mag) {m = mag;}
    void update_f(T field) {f = field;}
    void update_J(T j) {J = j;}
    void stream(const int width, const int height, const int evolution_per_frame);
    // void update_all();  // update all attributes after change of configuration


private:
    T J;  // interaction type: FM = -1; AFM = 1
    T t, e, f, m;  // temperature, energy, field, magnitization
    int h, w;  // N = h * w = height * width
    T spin;  // type of spin
    Matrix<T> Lat;
};

// member function definitions
template<class T>
void CIsing<T>::stream(const int width, const int height, const int evolution_per_frame) {

    Display display(width, height, "ShowQuads");
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    Visual<T> visual(Lat);
    visual.Frame(width,height);

    while (!display.IsClosed()) {
        glClear(GL_COLOR_BUFFER_BIT);
        evolve(evolution_per_frame);
        // printl();
        visual.Load();
        visual.Plot();
        display.Update();
    }
}


// no annealing for now
template<class T>
void CIsing<T>::evolve(const int nstep) {
    for (int i = 0; i < nstep; ++i)
        tryflip();
}

// try flip one spin
template<class T>
bool CIsing<T>::tryflip() {
    int iw = static_cast<int>(drand48() * w);
    int ih = static_cast<int>(drand48() * h);

    // suppose I flip the spin at (indw, indh)
    T delta_mag = -2 * (*this)(iw, ih);
    T delta_energy = - 2 * J * (*this)(iw, ih) * ((*this)(iw+1, ih) + (*this)(iw, ih+1)
                          + (*this)(iw-1, ih) + (*this)(iw, ih-1)) - f * delta_mag;

    // accept or reject
    if (delta_energy <=0 || drand48() < exp(-delta_energy/t)) {
        (*this)(iw, ih) *= -1;
        m += delta_mag;
        e += delta_energy;
        return true;
    }
    return false;
}


template<class T>
void CIsing<T>::config_init() {
    Lat.fillRand();
    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            // std::cout << "i = " << i << ", j = " << j << std::endl;
            if (Lat(i,j) < 0.5) {
                Lat(i,j) = 1.0;
                m += spin;
            }
            else{
                Lat(i,j) = -1.0;
                m -= spin;
            }
        }
    }
    T energy_init = energy();
    update_e(energy_init);
}

template<class T>
T CIsing<T>::energy() {
    /* calculate energy of current configuration  */
    int Nalig = 0;  // number of aligned n.n pairs
    int Nanti = 0;  // number of anti-alighed n.n paris

    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            // check two bonds for each site
            if ((*this)(i, j) == (*this)(i + 1, j))
                Nalig++;
            else
                Nanti++;
            if ((*this)(i, j) == (*this)(i, j + 1))
                Nalig++;
            else
                Nanti++;
        }
    }
    T E = J * Nalig * spin * spin- J * Nanti * spin * spin + f * m;
    return E;
}

#endif //MCMC_CISING_H
