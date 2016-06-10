#include <string>
#include <vector>
#include <deque>
#include <map>
#include <fstream>
#include <sstream>

#define PI 3.14159265f

#define SECTOR_SIZE 0.1f

#define ZERO glm::vec3(0.0f)
#define UNIT_T glm::vec3(0.0f, 1.0f, 0.0f)
#define UNIT_B glm::vec3(0.0f, -1.0f, 0.0f)
#define UNIT_L glm::vec3(-1.0f, 0.0f, 0.0f)
#define UNIT_R glm::vec3(1.0f, 0.0f, 0.0f)
#define UNIT_TL glm::vec3(-1.0f, 1.0f, 0.0f)
#define UNIT_TR glm::vec3(1.0f, 1.0f, 0.0f)
#define UNIT_BL glm::vec3(-1.0f, -1.0f, 0.0f)
#define UNIT_BR glm::vec3(1.0f, -1.0f, 0.0f)

struct ObjectData {
    float radius;
    float force;
    float mass;
    float maxVelocity;
    float elasticity;
    int health;
    int id;
    GLObjectData objectData;
};

enum SectorType {
    EMPTY,
    SQUARE,
    TRIANGLE_TL,
    TRIANGLE_TR,
    TRIANGLE_BL,
    TRIANGLE_BR,
    ONEWAY_UP,
    ONEWAY_DOWN,
    ONEWAY_LEFT,
    ONEWAY_RIGHT
};

GLObjectData CreateCircleData(float radius, int segments, const glm::vec3 &color) {
    GLObjectData data;
    data.vertices.reserve(3 * segments);
    glm::vec3 mid;
    mid.x = 0.0f;
    mid.y = 0.0f;
    mid.z = 0.0f;
    for (int i = 0; i < segments; i++) {
        glm::vec3 v1;
        v1.x = radius * cos(2 * PI * i / segments);
        v1.y = radius * sin(2 * PI * i / segments);
        v1.z = 0.0f;
        glm::vec3 v2;
        v2.x = radius * cos(2 * PI * (i + 1) / segments);
        v2.y = radius * sin(2 * PI * (i + 1) / segments);
        v2.z = 0.0f;
        data.vertices.push_back(mid);
        data.vertices.push_back(v1);
        data.vertices.push_back(v2);
        data.colors.push_back(color);
        data.colors.push_back(color * 0.3f);
        data.colors.push_back(color * 0.3f);
    }
    data.color = color;
    return data;
}

ObjectData CreateObjectData(int id) {
    ObjectData data;
    data.radius = 0.04f;
    data.force = 1.0f;
    data.mass = 1.0f;
    data.maxVelocity = 2.0f;
    data.elasticity = 0.8f;
    data.health = 100;
    data.id = id;
    int segments = 20 - id;
    glm::vec3 color;
    switch (id) {
    case 0:
        color = glm::vec3(1.0f, 0.0f, 0.0f);
        break;
    case 1:
        color = glm::vec3(0.0f, 1.0f, 0.0f);
        break;
    case 2:
        color = glm::vec3(0.0f, 0.0f, 1.0f);
        break;
    case 3:
        color = glm::vec3(1.0f, 1.0f, 0.0f);
        break;
    case 4:
        color = glm::vec3(1.0f, 0.0f, 1.0f);
        break;
    case 5:
        color = glm::vec3(0.0f, 1.0f, 1.0f);
        break;
    case 6:
        color = glm::vec3(1.0f, 1.0f, 1.0f);
        break;
    case 7:
        color = glm::vec3(0.7f, 0.2f, 0.5f);
        break;
    case 8:
        color = glm::vec3(0.2f, 0.9f, 0.2f);
        break;
    }
    data.objectData = CreateCircleData(data.radius, segments, color);
    return data;
}

GLObjectData CreateSquareData(float size, const glm::vec3 &color, SectorType type) {
    GLObjectData data;
    data.vertices.reserve(6);
    data.vertices.push_back(UNIT_TL * size);
    data.vertices.push_back(UNIT_TR * size);
    data.vertices.push_back(UNIT_BR * size);
    data.vertices.push_back(UNIT_TL * size);
    data.vertices.push_back(UNIT_BL * size);
    data.vertices.push_back(UNIT_BR * size);
    float tl = (type == ONEWAY_RIGHT || type == ONEWAY_DOWN) ? 0.1f : 1.0f;
    float tr = (type == ONEWAY_LEFT || type == ONEWAY_DOWN) ? 0.1f : 1.0f;
    float bl = (type == ONEWAY_RIGHT || type == ONEWAY_UP) ? 0.1f : 1.0f;
    float br = (type == ONEWAY_LEFT || type == ONEWAY_UP) ? 0.1f : 1.0f;
    data.colors.push_back(tl * color);
    data.colors.push_back(tr * color);
    data.colors.push_back(br * color);
    data.colors.push_back(tl * color);
    data.colors.push_back(bl * color);
    data.colors.push_back(br * color);
    data.color = color;
    return data;
}

GLObjectData CreateTriangleData(float size, const glm::vec3 &color, SectorType type) {
    GLObjectData data;
    data.vertices.reserve(3);
    if (type != TRIANGLE_BR)
        data.vertices.push_back(UNIT_TL * size);
    if (type != TRIANGLE_BL)
        data.vertices.push_back(UNIT_TR * size);
    if (type != TRIANGLE_TL)
        data.vertices.push_back(UNIT_BR * size);
    if (type != TRIANGLE_TR)
        data.vertices.push_back(UNIT_BL * size);
    data.colors.push_back(color);
    data.colors.push_back(color);
    data.colors.push_back(color);
    data.color = color;
    return data;
}

std::vector<std::string> ReadGridFromFile(const char *file) {
    std::vector<std::string> contents;
    std::string line;
    std::ifstream stream(file);
    // Surround borders of the grid with empty
    // sectors so collisions from "outside"
    // will be calculated correctly later
    contents.push_back("  ");
    if (stream.is_open()) {
        while (getline(stream, line)) {
            std::stringstream ss;
            ss << ' ' << line << ' ';
            contents.push_back(ss.str());
        }
        stream.close();
    }
    contents.push_back("  ");
    return contents;
}

std::deque<std::pair<std::string, unsigned int> > ReadSeason(const char *file) {
    std::deque<std::pair<std::string, unsigned int> > result;
    std::string line;
    std::ifstream stream(file);
    if (stream.is_open()) {
        while (getline(stream, line)) {
            size_t pos = line.find(',');
            if (pos != std::string::npos) {
                std::stringstream ss;
                ss << line.substr(pos + 1);
                unsigned int laps;
                ss >> laps;
                result.push_back(std::make_pair(line.substr(0, pos), laps));
                printf("%s -- %d\n", line.substr(0, pos).c_str(), laps);
            }
        }
    }
    return result;
}

std::map<uint64_t, float> ReadRecords(const char *file) {
    std::map<uint64_t, float> highscores;
    std::string line;
    std::ifstream stream(file);
    if (stream.is_open()) {
        while (getline(stream, line)) {
            size_t pos = line.find(',');
            if (pos != std::string::npos) {
                uint64_t hash;
                float time;
                std::stringstream ss;
                ss << line.substr(pos + 1);
                ss >> time;
                std::stringstream ss2;
                ss2 << line.substr(0, pos);
                ss2 >> hash;
                highscores[hash] = time;
            }
        }
        stream.close();
    }
    return highscores;
}

void WriteRecords(const char *file, const std::map<uint64_t, float> &highscores) {
    std::ofstream stream(file);
    if (stream.is_open()) {
        for (std::map<uint64_t, float>::const_iterator it = highscores.begin(); it != highscores.end(); ++it) {
            stream << it->first << "," << it->second << std::endl;
        }
        stream.close();
    }
}
