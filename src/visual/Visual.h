//
// Created by shifeng on 3/18/21.
//

#ifndef MCMC_VISUAL_H
#define MCMC_VISUAL_H

#include "../Matrix.h"
#include "Mesh.h"
#include "Display.h"
#include <glm/glm.hpp>
#include <vector>
#include <omp.h>

template<class T>
class Visual {
public:
    Visual(Matrix<T>& m)  : width(800),
                            height(800),
                            pixel_h(0),
                            pixel_w(0),
                            indx(0),
                            indy(0),
                            num2draw(0) {
        config = &m;
    }

    virtual ~Visual() {}

    void Frame(int width, int height) {
        this -> width = width;
        this -> height = height;
        pixel_w = 2.0 / config -> cols();
        pixel_h = 2.0 / config -> rows();
    }

    void Load();  // load vertex array from config matrix
    void Show();  // load vertex array and draw a static config
    void Plot();  // draw a single frame using current vertex array in gl buffer

private:
    std::vector<int> indx;
    std::vector<int> indy;
    int width;
    int height;
    int num2draw;
    T pixel_w;
    T pixel_h;
    Matrix<T> *config;
    std::vector<Vertex> vertices;
};


// member function definitions
template <class T>
void Visual<T>::Load() {
    indx.clear();
    indy.clear();
    vertices.clear();

    int row = config->rows();
    int col = config->cols();

    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            if ((*config)(i,j) > 0.0) {  // change later
                indx.push_back(i);
                indy.push_back(j);
            }
        }
    }
    num2draw = indx.size();

    vertices.reserve(4*num2draw);
//    vertices.resize(4*num2draw);

//        #pragma omp parallel for num_threads(4) default(none)
    for (int i = 0; i < num2draw; ++i) {
//            // thread safe
//            vertices[i*4] = (Vertex(glm::vec3(-1 + indx[i] * pixel_w, -1 + indy[i] * pixel_h, 0)));
//            vertices[1+i*4] = (Vertex(glm::vec3(-1 + (indx[i] + 1) * pixel_w, -1 + indy[i] * pixel_h, 0)));
//            vertices[2+i*4] = (Vertex(glm::vec3(-1 + (indx[i] + 1) * pixel_w, -1 + (indy[i] + 1) * pixel_h, 0)));
//            vertices[3+i*4] = (Vertex(glm::vec3(-1 + indx[i] * pixel_w, -1 + (indy[i] + 1) * pixel_h, 0)));

        // thread-unsafe
        vertices.emplace_back(Vertex(glm::vec3(-1 + indx[i] * pixel_w, -1 + indy[i] * pixel_h, 0)));
        vertices.emplace_back(Vertex(glm::vec3(-1 + (indx[i] + 1) * pixel_w, -1 + indy[i] * pixel_h, 0)));
        vertices.emplace_back(Vertex(glm::vec3(-1 + (indx[i] + 1) * pixel_w, -1 + (indy[i] + 1) * pixel_h, 0)));
        vertices.emplace_back(Vertex(glm::vec3(-1 + indx[i] * pixel_w, -1 + (indy[i] + 1) * pixel_h, 0)));

//            std::cout << "( " << -1 + indx[i] * pixel_w << "," << -1 + indy[i] * pixel_h << ")" << std::endl;
//            std::cout << "( " << -1 + (indx[i] + 1) * pixel_w << "," <<  -1 + indy[i] * pixel_h << ")" << std::endl;
//            std::cout << "( " << -1 + (indx[i] + 1) * pixel_w << "," <<  -1 + (indy[i] + 1) * pixel_h << ")" << std::endl;
//            std::cout << "( " << -1 + indx[i] * pixel_w << "," <<  -1 + (indy[i] + 1) * pixel_h << ")" << std::endl << std::endl;
    }

}

template <class T>
void Visual<T>::Plot()  {
    Mesh mesh(&(vertices[0]), 4*num2draw);
    glClear(GL_COLOR_BUFFER_BIT);
    mesh.Draw_quads();
}

template <class T>
void Visual<T>::Show() {
    Load();
    Display display(width, height, "ShowQuads");
    Mesh mesh(&(vertices[0]), vertices.size());
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    while (!display.IsClosed()) {
        glClear(GL_COLOR_BUFFER_BIT);
        mesh.Draw_quads();
        display.Update();
    }
}

#endif //MCMC_VISUAL_H
