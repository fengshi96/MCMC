//
// Created by shifeng on 3/18/21.
//

#ifndef MCMC_MESH_H
#define MCMC_MESH_H


#include <glm/glm.hpp>
#include <GL/glew.h>


// class representing data of vertices
class Vertex {
public:
    Vertex() : pos(0,0,0) {};

    Vertex(const glm::vec3& pos) {
        this -> pos = pos;
    }

    Vertex(Vertex&&) = default;  // enforcing move constructor

    Vertex(Vertex& vertex){
        pos = *vertex.GetPos();
    }

    inline glm::vec3* GetPos() { return &pos; }
private:
    glm::vec3 pos;
};


// class of mesh
class Mesh {
public:
    Mesh(Vertex* vertices, unsigned int numVertices);
    void Draw();  // take the mesh and draw with GPU (drop into pipline and make into image)
    void Draw_quads();
    virtual ~Mesh();

private:
    Mesh(const Mesh& other);
    void operator=(const Mesh& other);

    // user-defined data type that consists of integral constants
    enum {
        POSITION_VB,
        NUM_BUFFERS
    };

    GLuint m_vertexArrayObject;
    GLuint m_vertexArrayBuffers[NUM_BUFFERS];
    unsigned int m_drawCount;  // how many things to draw
};



#endif //MCMC_MESH_H
