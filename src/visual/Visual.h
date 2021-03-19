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

template<class T>
class Visual {
public:
    Visual(Matrix<T>& m)  : width(800), height(800), pixel_h(0), pixel_w(0),
    indx(0), indy(0), num2draw(0) { config = &m;  }

    virtual ~Visual() {}

    void Frame(int width, int height) {
        this -> width = width;
        this -> height = height;
        pixel_w = 2.0 / config -> cols();
        pixel_h = 2.0 / config -> rows();
    }

    void Load() {
        indx.resize(0);
        indy.resize(0);
        vertices.resize(0);

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
//        std::cout << "num2draw = " <<  num2draw << std::endl;

        for (int i = 0; i < num2draw; ++i) {
            vertices.emplace_back(Vertex(glm::vec3(-1 + indx[i] * pixel_w, -1 + indy[i] * pixel_h, 0)));
            vertices.emplace_back(Vertex(glm::vec3(-1 + (indx[i] + 1) * pixel_w, -1 + indy[i] * pixel_h, 0)));
            vertices.emplace_back(Vertex(glm::vec3(-1 + (indx[i] + 1) * pixel_w, -1 + (indy[i] + 1) * pixel_h, 0)));
            vertices.emplace_back(Vertex(glm::vec3(-1 + indx[i] * pixel_w, -1 + (indy[i] + 1) * pixel_h, 0)));

//            std::cout << "( " << -1 + indx[i] * pixel_w << "," << -1 + indy[i] * pixel_h << ")" << std::endl;
//            std::cout << "( " << -1 + (indx[i] + 1) * pixel_w << "," <<  -1 + indy[i] * pixel_h << ")" << std::endl;
//            std::cout << "( " << -1 + (indx[i] + 1) * pixel_w << "," <<  -1 + (indy[i] + 1) * pixel_h << ")" << std::endl;
//            std::cout << "( " << -1 + indx[i] * pixel_w << "," <<  -1 + (indy[i] + 1) * pixel_h << ")" << std::endl;
//            std::cout << std::endl;
//            std::cout << std::endl;

        }

    }

    void Show() {
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

    void Plot() {
        Mesh mesh(&(vertices[0]), 4*num2draw);
        glClear(GL_COLOR_BUFFER_BIT);
        mesh.Draw_quads();
    }

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

#endif //MCMC_VISUAL_H
