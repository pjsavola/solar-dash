#ifndef SD_GL_H
#define SD_GL_H

#define GLM_FORCE_RADIANS

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/gtx/transform.hpp>
#include <vector>
#include <string>

struct GLObjectData {
    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> colors;
};

// Base class for visible objects
class GLObject {
protected:
    GLObject(const GLObjectData &data);
    ~GLObject();
    void DrawAt(GLuint id, const glm::mat4 &model, const glm::mat4 &view) const;

private:
    unsigned int bufferSize;
    GLuint vertexbuffer;
    GLuint colorbuffer;
};

class Program {
public:
    Program(int width, int height, const std::string &title);
    ~Program();
    int Run() const;

private:
    GLFWwindow *window;
    GLuint VertexArrayID;
    GLuint programID;
    GLuint MatrixID;
    bool initialized;
};

#endif // SD_GL_H
