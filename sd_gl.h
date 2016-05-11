#ifndef SD_GL_H
#define SD_GL_H

#define GLM_FORCE_RADIANS

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/gtx/transform.hpp>

#include <ft2build.h>
#include FT_FREETYPE_H

#include <vector>
#include <string>
#include <map>

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

class FontRenderer {
private:
    struct Character {
        GLuint     TextureID;  // ID handle of the glyph texture
        glm::ivec2 Size;       // Size of glyph
        glm::ivec2 Bearing;    // Offset from baseline to left/top of glyph
        GLuint     Advance;    // Offset to advance to next glyph
    };

    std::map<GLchar, Character> Characters;

    FT_Library ft;
    FT_Face face;
    GLuint VBO;
    GLuint VAO;
    GLuint program;

public:
    ~FontRenderer();
    void Init();
    void RenderText(std::string text, GLfloat x, GLfloat y, GLfloat scale, glm::vec3 color) const;
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
    FontRenderer fontRenderer;
};

#endif // SD_GL_H
