#ifndef SD_GL_H
#define SD_GL_H

#define GLM_FORCE_RADIANS

#include <GL/glew.h>
#include <glfw3.h>
#include <glm/gtx/transform.hpp>

#include <ft2build.h>
#include FT_FREETYPE_H

#include <vector>
#include <string>
#include <map>

struct GLObjectData {
    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> colors;
    glm::vec3 color;
};

// Base class for visible objects
class GLObject {
protected:
    GLObject(const GLObjectData &data);
    ~GLObject();
    void DrawAt(GLuint id, const glm::mat4 &model, const glm::mat4 &view) const;
public:
    const glm::vec3& GetColor() const {
        return color;
    }
	void AdjustBrightness(float newBrightness);
private:
    unsigned int bufferSize;
    GLuint vertexbuffer;
    GLuint colorbuffer;
    glm::vec3 color;
	const std::vector<glm::vec3> colorData;
};

class Shader {
public:
    Shader(bool blend) : blend(blend) { }
    void Load(const char *vs_path, const char *fs_path);
    void Unload() {
        glDeleteProgram(id);
    }
    void Use() const {
        glUseProgram(id);
        if (blend) {
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        } else {
            glDisable(GL_BLEND);
        }
    }
    GLuint Get(const char *name) const {
        return glGetUniformLocation(id, name);
    }
private:
    GLuint id;
    bool blend;
};

class TextRenderer {
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
    Shader shader;

public:
    TextRenderer() : shader(true) { }
    ~TextRenderer();
    void Init();
    void RenderText(std::string text, GLfloat x, GLfloat y, GLfloat scale, glm::vec3 color) const;
};

class Program {
public:
    Program(int width, int height, const std::string &title);
    ~Program();
    int Run(const char *input) const;

private:
    GLFWwindow *window;
    GLuint VertexArrayID;
    GLuint MatrixID;
    bool initialized;
    TextRenderer textRenderer;
    Shader shader;
};

#endif // SD_GL_H
