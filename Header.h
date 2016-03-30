#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <glm/glm.hpp>

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

struct GLObjectData {
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> colors;
};

struct ObjectData {
	float radius;
	float force;
	float mass;
	float maxVelocity;
	float elasticity;
	GLObjectData objectData;
};

enum SectorType { EMPTY, SQUARE, TRIANGLE_TL, TRIANGLE_TR, TRIANGLE_BL, TRIANGLE_BR,
	ONEWAY_UP, ONEWAY_DOWN, ONEWAY_LEFT, ONEWAY_RIGHT };

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

ObjectData CreateObjectData(int id) {
	ObjectData data;
	data.radius = 0.04f;
	data.force = 1.0f;
	data.mass = 1.0f;
	data.maxVelocity = 2.0f;
	data.elasticity = 0.8f;
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

GLObjectData CreateSquareData(float size, const glm::vec3 &color) {
	GLObjectData data;
	data.vertices.reserve(6);
	data.vertices.push_back(UNIT_TL * size * 0.5f);
	data.vertices.push_back(UNIT_TR * size * 0.5f);
	data.vertices.push_back(UNIT_BR * size * 0.5f);
	data.vertices.push_back(UNIT_TL * size * 0.5f);
	data.vertices.push_back(UNIT_BL * size * 0.5f);
	data.vertices.push_back(UNIT_BR * size * 0.5f);
	data.colors.push_back(color);
	data.colors.push_back(color);
	data.colors.push_back(color);
	data.colors.push_back(color);
	data.colors.push_back(color);
	data.colors.push_back(color);
	return data;
}

GLObjectData CreateTriangleData(float size, const glm::vec3 &color, SectorType type) {
	GLObjectData data;
	data.vertices.reserve(3);
	if (type != TRIANGLE_BR)
		data.vertices.push_back(UNIT_TL * size * 0.5f);
	if (type != TRIANGLE_BL)
		data.vertices.push_back(UNIT_TR * size * 0.5f);
	if (type != TRIANGLE_TL)
		data.vertices.push_back(UNIT_BR * size * 0.5f);
	if (type != TRIANGLE_TR)
		data.vertices.push_back(UNIT_BL * size * 0.5f);
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

