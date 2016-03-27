#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <glm/glm.hpp>

#define PI 3.14159265f

struct GLObjectData {
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> colors;
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
	return data;
}

GLObjectData CreateSquareData(float size, const glm::vec3 &color) {
	GLObjectData data;
	data.vertices.reserve(6);
	data.vertices.push_back(glm::vec3(-size / 2, size / 2, 0.0f));
	data.vertices.push_back(glm::vec3(size / 2, size / 2, 0.0f));
	data.vertices.push_back(glm::vec3(size / 2, -size / 2, 0.0f));
	data.vertices.push_back(glm::vec3(-size / 2, size / 2, 0.0f));
	data.vertices.push_back(glm::vec3(-size / 2, -size / 2, 0.0f));
	data.vertices.push_back(glm::vec3(size / 2, -size / 2, 0.0f));
	data.colors.push_back(color);
	data.colors.push_back(color);
	data.colors.push_back(color);
	data.colors.push_back(color);
	data.colors.push_back(color);
	data.colors.push_back(color);
	return data;
}

enum TriangleType { TL, TR, BL, BR };

GLObjectData CreateTriangleData(float size, const glm::vec3 &color, TriangleType type) {
	GLObjectData data;
	data.vertices.reserve(3);
	if (type != BR)
		data.vertices.push_back(glm::vec3(-size / 2, size / 2, 0.0f));
	if (type != BL)
		data.vertices.push_back(glm::vec3(size / 2, size / 2, 0.0f));
	if (type != TL)
		data.vertices.push_back(glm::vec3(size / 2, -size / 2, 0.0f));
	if (type != TR)
		data.vertices.push_back(glm::vec3(-size / 2, -size / 2, 0.0f));
	data.colors.push_back(color);
	data.colors.push_back(color);
	data.colors.push_back(color);
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