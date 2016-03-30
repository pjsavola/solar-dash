#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <GL/glew.h>

#include <glfw3.h>
GLFWwindow* window;

#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>

#include <common/shader.hpp>

#include "playground/Header.h"

#include <deque>
#include <vector>
#include <map>
#include <set>
using namespace std;


// Base class for visible objects
class GLObject {
protected:
	GLObject(const GLObjectData &data) {
		assert(data.vertices.size() == data.colors.size());
		bufferSize = data.vertices.size();

		glGenBuffers(1, &vertexbuffer);
		glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
		glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * bufferSize, &data.vertices[0], GL_STATIC_DRAW);

		glGenBuffers(1, &colorbuffer);
		glBindBuffer(GL_ARRAY_BUFFER, colorbuffer);
		glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * bufferSize, &data.colors[0], GL_STATIC_DRAW);
	}

	~GLObject() {
		glDeleteBuffers(1, &vertexbuffer);
		glDeleteBuffers(1, &colorbuffer);
	}

	void DrawAt(GLuint id, const glm::mat4 &model, const glm::mat4 &view) const {
		const glm::mat4 MVP = view * model;
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
		glUniformMatrix4fv(id, 1, GL_FALSE, &MVP[0][0]);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, colorbuffer);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);

		glDrawArrays(GL_TRIANGLES, 0, bufferSize);
		glDisableVertexAttribArray(0);
		glDisableVertexAttribArray(1);
	}

private:
	unsigned int bufferSize;
	GLuint vertexbuffer;
	GLuint colorbuffer;
};

class AI {

};

class DummyAI : public AI {

};

// Circle shaped object
class Object : public GLObject {
public:
	Object(const ObjectData &data, const glm::vec3 &initialLocation) : GLObject(data.objectData) {
		radius = data.radius;
		force = data.force;
		mass = data.mass;
		maxVelocity = data.maxVelocity;
		elasticity = data.elasticity;
		location = initialLocation;
		velocity = glm::vec3(0.0f);
		acceleration = data.force / data.mass;
		ai = 0;
	}

	void Draw(GLuint id, const glm::mat4 &view) const {
		DrawAt(id, glm::translate(location), view);
	}

	const glm::vec3& GetLocation() const {
		return location;
	}

	const glm::vec3& GetVelocity() const {
		return velocity;
	}

	float GetAcceleration() const {
		return acceleration;
	}

	float GetRadius() const {
		return radius; 
	}

	void SetAI(const AI * const newAI) {
		ai = newAI;
	}

	void Move(float time) {
		// Save the previous location because the object may collide in the new location.
		previousLocation = location;
		location += time * velocity;
	}

	void Accelerate(float time, const glm::vec3 &directions) {
		velocity += time * acceleration * directions;
		if (glm::length(velocity) > maxVelocity) {
			velocity = glm::normalize(velocity) * maxVelocity;
		}
	}

	bool Collides(const Object * const o) const {
		float len = glm::length(o->location - location);
		return len < radius + o->radius;
	}

	// Resolve collision for both objects.
	// Collision is currently assumed to be fully elastic.
	void ResolveCollision(Object * const o) {
		assert(Collides(o));

		const glm::vec3 ncoll = glm::normalize(o->location - location);

		const float u1 = glm::dot(velocity, ncoll);
		const float u2 = glm::dot(o->velocity, ncoll);

		const float v1 = (u1 * (mass - o->mass) + 2 * o->mass * u2) / (mass + o->mass);
		const float v2 = (u2 * (o->mass - mass) + 2 * mass * u1) / (mass + o->mass);

		// Change velocities ...
		velocity += (v1 - u1) * ncoll;
		o->velocity += (v2 - u2) * ncoll;

		// ... but prevent movement.
		location = previousLocation;
		o->location = o->previousLocation;
	}

	// Resolve collision for this object.
	// Collision is not fully elastic.
	void ResolveCollision(const glm::vec3 &normal) {
		velocity = velocity - 2.0f * glm::dot(velocity, normal) * normal;
		velocity *= elasticity;

		// Prevent movement.
		location = previousLocation;
	}

	const AI* GetAI() const {
		return ai;
	}

private:
	float radius;
	float force;
	float mass;
	float acceleration;
	float maxVelocity;
	float elasticity;
	const AI *ai;
	glm::vec3 location;
	glm::vec3 previousLocation;
	glm::vec3 velocity;
};

class Sector {
public:
	Sector(const glm::vec3 &location) : model(glm::translate(location)), location(location) { }

	virtual ~Sector() { }

	virtual void Draw(GLuint id, const glm::mat4 &view) const { }

	virtual bool IsSolid(const glm::vec3 &velocity) const {
		return false;
	}

	virtual SectorType GetType() const {
		return EMPTY;
	}

	const glm::vec3& GetLocation() const {
		return location;
	}

protected:
	const glm::mat4 model;

private:
	const glm::vec3 location;
};

class SolidSector : public Sector, public GLObject {
public:
	void Draw(GLuint id, const glm::mat4 &view) const {
		DrawAt(id, model, view);
	}

	virtual bool IsSolid(const glm::vec3 &velocity) const {
		return true;
	}

	virtual SectorType GetType() const = 0;

	virtual const vector<glm::vec3>& GetCorners() const = 0;

protected:
	SolidSector(const glm::vec3 &location, const GLObjectData &data) :
		GLObject(data), Sector(location) { }
};

class SquareSector : public SolidSector {
public:
	SquareSector(const glm::vec3 &location, float size, const glm::vec3 &color) :
		SolidSector(location, CreateSquareData(size, color)) {
		corners.reserve(4);
		corners.push_back(location + UNIT_TL * size * 0.5f);
		corners.push_back(location + UNIT_TR * size * 0.5f);
		corners.push_back(location + UNIT_BR * size * 0.5f);
		corners.push_back(location + UNIT_BL * size * 0.5f);
	}

	const vector<glm::vec3>& GetCorners() const {
		return corners;
	};

	SectorType GetType() const {
		return SQUARE;
	}

private:
	vector<glm::vec3> corners;
};

class OneWaySector : public SquareSector {
public:
	OneWaySector(const glm::vec3 &location, float size, const glm::vec3 &color, const glm::vec3 &way) :
		SquareSector(location, size, color), way(way) {
		if (way == UNIT_T) type = ONEWAY_UP;
		else if (way == UNIT_B) type = ONEWAY_DOWN;
		else if (way == UNIT_L) type = ONEWAY_LEFT;
		else if (way == UNIT_R) type = ONEWAY_RIGHT;
		else assert(false);
	}

	bool IsSolid(const glm::vec3 &velocity) const {
		return glm::dot(velocity, way) < 0;
	}

	SectorType GetType() const {
		return type;
	}

private:
	const glm::vec3 way;
	SectorType type;
};

class TriangleSector : public SolidSector {
public:
	TriangleSector(const glm::vec3 &location, float size, const glm::vec3 &color, SectorType type) :
		SolidSector(location, CreateTriangleData(size, color, type)), type(type) {
		corners.reserve(3);
		if (type != TRIANGLE_BR)
			corners.push_back(location + UNIT_TL * size * 0.5f);
		if (type != TRIANGLE_BL)
			corners.push_back(location + UNIT_TR * size * 0.5f);
		if (type != TRIANGLE_TL)
			corners.push_back(location + UNIT_BR * size * 0.5f);
		if (type != TRIANGLE_TR)
			corners.push_back(location + UNIT_BL * size * 0.5f);
	}

	const vector<glm::vec3>& GetCorners() const {
		return corners;
	};

	SectorType GetType() const {
		return type;
	}
private:
	vector<glm::vec3> corners;
	const SectorType type;
};

const Sector* CreateSector(char c, vector<Object *> &objects,
	const glm::vec3 &location, float sectorSize, const glm::vec3 &color) {
	if (c >= '1' && c <= '9') {
		int id = c - '1';
		objects.push_back(new Object(CreateObjectData(id), location));
	}
	const Sector *sector = 0;
	switch (c) {
	case '#':
		sector = new SquareSector(location, sectorSize, color);
		break;
	case 'A':
		sector = new TriangleSector(location, sectorSize, color, TRIANGLE_TL);
		break;
	case 'B':
		sector = new TriangleSector(location, sectorSize, color, TRIANGLE_TR);
		break;
	case 'C':
		sector = new TriangleSector(location, sectorSize, color, TRIANGLE_BL);
		break;
	case 'D':
		sector = new TriangleSector(location, sectorSize, color, TRIANGLE_BR);
		break;
	case '>':
		sector = new OneWaySector(location, sectorSize, color * 0.2f, UNIT_R);
		break;
	case '<':
		sector = new OneWaySector(location, sectorSize, color * 0.2f, UNIT_L);
		break;
	case '^':
		sector = new OneWaySector(location, sectorSize, color * 0.2f, UNIT_T);
		break;
	case 'v':
		sector = new OneWaySector(location, sectorSize, color * 0.2f, UNIT_B);
		break;
	}
	if (!sector) {
		sector = new Sector(location);
	}
	return sector;
}

/*
Grid has several purposes:
1. Make map creation easy
2. Help AI in route calculation
3. Compute collisions against wall efficiently
*/
class Grid {
public:
	Grid(float sectorSize) : sectorSize(sectorSize) { }

	void Initialize(const vector<string> &data, vector<Object *> &objects) {

		// Calculate the grid size based on data read from a map file
		unsigned int cols = 0;
		const unsigned int rows = data.size();
		for (vector<string>::const_iterator it = data.begin(); it != data.end(); ++it) {
			const string &row = *it;
			if (cols < row.length()) {
				cols = row.length();
			}
		}

		// Create all objects for the map based on input data
		sectors.reserve(cols);
		for (unsigned int i = 0; i < cols; i++) {
			vector<const Sector *> column;
			column.reserve(rows);
			for (unsigned int j = 0; j < rows; j++) {
				const glm::vec3 location(sectorSize * i, -sectorSize * j, 0.0f);
				const string &row = data.at(j);
				float green = glm::min(0.1f + 1.0f * i / cols, 0.9f);
				float blue = glm::min(0.1f + 1.0f * j / rows, 0.9f);
				const glm::vec3 color(0.0f, green, blue);
				const char c = row.length() > i ? row.at(i) : ' ';
				column.push_back(CreateSector(c, objects, location, sectorSize, color));
			}
			sectors.push_back(column);
		}
		
		// Calculate neighbor sectors for each sector
		for (unsigned int i = 0; i < cols; i++) {
			for (unsigned int j = 0; j < rows; j++) {

				// Solid sectors have no neighbors
				if (sectors.at(i).at(j)->IsSolid(ZERO)) {
					continue;
				}

				vector<pair<int, int>> neighbors;
				int tl = 0;
				int tr = 0;
				int bl = 0;
				int br = 0;
				if (j > 0) {
					if (!sectors.at(i).at(j - 1)->IsSolid(glm::vec3(0.0f, 1.0f, 0.0f))) {
						neighbors.push_back(std::make_pair(i, j - 1));
						tl++;
						tr++;
					}
					else {
						const TriangleSector *triangle = dynamic_cast<const TriangleSector *>(sectors.at(i).at(j - 1));
						if (triangle) {
							if (triangle->GetType() == TRIANGLE_TL) {
								tr++;
								neighbors.push_back(std::make_pair(i, j - 1));
							}
							else if (triangle->GetType() == TRIANGLE_TR) {
								tl++;
								neighbors.push_back(std::make_pair(i, j - 1));
							}
						}
					}
				}
				if (i > 0) {
					if (!sectors.at(i - 1).at(j)->IsSolid(glm::vec3(-1.0f, 0.0f, 0.0f))) {
						neighbors.push_back(std::make_pair(i - 1, j));
						tl++;
						bl++;
					}
					else {
						const TriangleSector *triangle = dynamic_cast<const TriangleSector *>(sectors.at(i - 1).at(j));
						if (triangle) {
							if (triangle->GetType() == TRIANGLE_TL) {
								bl++;
								neighbors.push_back(std::make_pair(i - 1, j));
							}
							else if (triangle->GetType() == TRIANGLE_BL) {
								tl++;
								neighbors.push_back(std::make_pair(i - 1, j));
							}
						}
					}
				}
				if (j < sectors.at(i).size() - 1) {
					if (!sectors.at(i).at(j + 1)->IsSolid(glm::vec3(0.0f, -1.0f, 0.0f))) {
						neighbors.push_back(std::make_pair(i, j + 1));
						bl++;
						br++;
					}
					else {
						const TriangleSector *triangle = dynamic_cast<const TriangleSector *>(sectors.at(i).at(j + 1));
						if (triangle) {
							if (triangle->GetType() == TRIANGLE_BL) {
								br++;
								neighbors.push_back(std::make_pair(i, j + 1));
							}
							else if (triangle->GetType() == TRIANGLE_BR) {
								bl++;
								neighbors.push_back(std::make_pair(i, j + 1));
							}
						}
					}
				}
				if (i < sectors.size() - 1) {
					if (!sectors.at(i + 1).at(j)->IsSolid(glm::vec3(1.0f, 0.0f, 0.0f))) {
						neighbors.push_back(std::make_pair(i + 1, j));
						br++;
						tr++;
					}
					else {
						const TriangleSector *triangle = dynamic_cast<const TriangleSector *>(sectors.at(i + 1).at(j));
						if (triangle) {
							if (triangle->GetType() == TRIANGLE_TR) {
								br++;
								neighbors.push_back(std::make_pair(i + 1, j));
							}
							else if (triangle->GetType() == TRIANGLE_BR) {
								tr++;
								neighbors.push_back(std::make_pair(i + 1, j));
							}
						}
					}
				}
				if (tr == 2 && !sectors.at(i + 1).at(j - 1)->IsSolid(glm::vec3(1.0f, -1.0f, 0.0f))) neighbors.push_back(std::make_pair(i + 1, j - 1));
				if (tl == 2 && !sectors.at(i - 1).at(j - 1)->IsSolid(glm::vec3(-1.0f, -1.0f, 0.0f))) neighbors.push_back(std::make_pair(i - 1, j - 1));
				if (br == 2 && !sectors.at(i + 1).at(j + 1)->IsSolid(glm::vec3(1.0f, 1.0f, 0.0f))) neighbors.push_back(std::make_pair(i + 1, j + 1));
				if (bl == 2 && !sectors.at(i - 1).at(j + 1)->IsSolid(glm::vec3(-1.0f, 1.0f, 0.0f))) neighbors.push_back(std::make_pair(i - 1, j + 1));
				neighborMap[std::make_pair(i, j)] = neighbors;
			}
		}

		CalculateDistances();
	}

	~Grid() {
		for (unsigned int i = 0; i < sectors.size(); i++) {
			for (unsigned int j = 0; j < sectors[i].size(); j++) {
				delete sectors.at(i).at(j);
			}
		}
	}

	void Draw(GLuint id, const glm::mat4 &camera) const {
		for (unsigned int i = 0; i < sectors.size(); i++) {
			for (unsigned int j = 0; j < sectors.at(i).size(); j++) {
				sectors.at(i).at(j)->Draw(id, camera);
			}
		}
	}
	
	glm::vec3 GetCollisionNormal(const glm::vec3 &location, const glm::vec3 &velocity, float radius) const {
		glm::vec3 normal(0.0f);
		const std::pair<int, int> p = GetColRow(location);
		const Sector *sector = GetSector(p.first, p.second);
		if (!sector) return normal;
		std::vector<const SolidSector *> collisionCandidates;
		AddCollisionCandidates(location, velocity, radius, sector, p, collisionCandidates);
		for (std::vector<const SolidSector *>::const_iterator it = collisionCandidates.begin();
			it != collisionCandidates.end(); ++it) {
			const std::vector<glm::vec3> &corners = (*it)->GetCorners();
			float minDistance = radius;

			// Check distance to edges
			for (unsigned int i = 0; i < corners.size(); i++) {
				const glm::vec3 &c1 = corners.at(i);
				const glm::vec3 &c2 = corners.at((i + 1) % corners.size());
				const glm::vec3 edge(c1 - c2);
				glm::vec3 edgeNormal = glm::normalize(glm::vec3(-edge.y, edge.x, 0.0f));
				const glm::vec3 v1 = location - c1;
				float distance = glm::dot(v1, edgeNormal);

				// If we didn't happen to pick the right normal, just change the sign.
				if (distance < 0) {
					edgeNormal = -edgeNormal;
					distance = -distance;
				}

				if (distance < minDistance) {
					const glm::vec3 v2 = location - c2;
					// Feasible sector needed for an edge collision. We guarantee this by checking
					// that angle to point x is under 90 degrees from both corners.
					if (glm::dot(edge, v2) > 0 && glm::dot(-edge, v1) > 0) {
						minDistance = distance;
						normal = edgeNormal;
					}
				}
			}

			if (normal == glm::vec3(0.0f)) {
				// Check distance to corners
				for (std::vector<glm::vec3>::const_iterator it2 = corners.begin(); it2 != corners.end(); ++it2) {
					const glm::vec3 candidate = location - *it2;
					const float distance = glm::length(candidate);
					if (distance < minDistance) {
						minDistance = distance;
						normal = glm::normalize(candidate);
					}
				}
			}
		}
		return normal;
	}

	glm::vec3 GetCollisionNormal(const Object *o) const {
		return GetCollisionNormal(o->GetLocation(), o->GetVelocity(), o->GetRadius());
	}
private:
	// Returns a pair (col, row) which corresponds to the given location
	std::pair<int, int> GetColRow(const glm::vec3 &location) const {
		float x = sectorSize / 2 + location.x;
		float y = sectorSize / 2 - location.y;
		int col = (int) (x / sectorSize);
		int row = (int) (y / sectorSize);
		return std::make_pair(col, row);
	}

	void AddCollisionCandidates(const glm::vec3 &location, const glm::vec3 &velocity, const float radius, const Sector *sector, const std::pair<int, int> &p, std::vector<const SolidSector *> &collisionCandidates) const {
		const glm::vec3 diff = sector->GetLocation() - location;
		int colDelta = 0;
		int rowDelta = 0;
		if (abs(diff.x) + radius > sectorSize * 0.5f) {
			colDelta = diff.x > 0 ? -1 : 1;
		}
		if (abs(diff.y) + radius > sectorSize * 0.5f) {
			rowDelta = diff.y > 0 ? 1 : -1;
		}
		if (sector->IsSolid(velocity)) {
			collisionCandidates.push_back(static_cast<const SolidSector *>(sector));
		}
		if (colDelta != 0) {
			sector = GetSector(p.first + colDelta, p.second);
			if (sector && sector->IsSolid(velocity)) {
				collisionCandidates.push_back(static_cast<const SolidSector *>(sector));
			}
			if (rowDelta != 0) {
				sector = GetSector(p.first + colDelta, p.second + rowDelta);
				if (sector && sector->IsSolid(velocity)) {
					collisionCandidates.push_back(static_cast<const SolidSector *>(sector));
				}
			}
		}
		if (rowDelta != 0) {
			sector = GetSector(p.first, p.second + rowDelta);
			if (sector && sector->IsSolid(velocity)) {
				collisionCandidates.push_back(static_cast<const SolidSector *>(sector));
			}
		}
	}

public:
	void CalculateDistances() {
		std::deque<std::pair<int, int>> worklist;
		for (unsigned int i = 0; i < sectors.size(); i++) {
			for (unsigned int j = 0; j < sectors.at(i).size(); j++) {
				const Sector *sector = sectors.at(i).at(j);
				const OneWaySector *ows = dynamic_cast<const OneWaySector *>(sector);
				if (ows) {
					distanceMap[ows] = 0;
					if (ows->GetType() == ONEWAY_RIGHT) {
						worklist.push_back(std::make_pair(i - 1, j));
					}
					else if (ows->GetType() == ONEWAY_LEFT) {
						worklist.push_back(std::make_pair(i + 1, j));
					}
					else if (ows->GetType() == ONEWAY_UP) {
						worklist.push_back(std::make_pair(i, j + 1));
					}
					else if (ows->GetType() == ONEWAY_DOWN) {
						worklist.push_back(std::make_pair(i, j - 1));
					}
				}
			}
		}

		int d = 1;
		// Separator for increasing distance
		static std::pair<int, int> sep = std::make_pair(-1, -1);
		worklist.push_back(sep);
		while (!worklist.empty()) {
			const std::pair<int, int> p = worklist.front();
			worklist.pop_front();
			if (p == sep) {
				d++;
				if (!worklist.empty()) {
					worklist.push_back(sep);
				}
				continue;
			}
			const Sector *sector = GetSector(p.first, p.second);
			if (!sector) {
				continue;
			}
			if (distanceMap.find(sector) != distanceMap.end()) {
				// Already covered earlier, possibly with less distance
				continue;
			}
			else {
				distanceMap[sector] = d;
				const std::map<std::pair<int, int>, std::vector<std::pair<int, int>>>::const_iterator it = neighborMap.find(p);
				if (it != neighborMap.end()) {
					const std::vector<std::pair<int, int>> &neighbors = it->second;
					for (std::vector<std::pair<int, int>>::const_iterator vit = neighbors.begin();
						vit != neighbors.end(); ++vit) {
						worklist.push_back(*vit);
					}
				}
			}
		}
		//printf("Init done\n");
	}

	void Debug(const Object *o) const {
		std::pair<int, int> p = GetColRow(o->GetLocation());
		const Sector *s = GetSector(p.first, p.second);
		if (s && distanceMap.find(s) != distanceMap.end()) {
			//printf("Distance: %d Velocity: %f\n", distanceMap.find(s)->second, length(o->GetVelocity()));
		}
	}

	std::map<const Object *, glm::vec3> locMap;
	std::map<const Object *, glm::vec3> dirMap;
	std::map<const Object *, float> timeMap;

	glm::vec3 GetAcceleration(const Object *o, float deltaTime) {
		if (dirMap.find(o) != dirMap.end()) {
			float time = timeMap[o];
			glm::vec3 dir = dirMap.find(o)->second;
			time -= deltaTime;
			if (time < 0) {
				timeMap.clear();
				dirMap.clear();
			}
			else {
				timeMap[o] = time;
			}
			return dir;
		}
		float a = o->GetAcceleration();
		const glm::vec3 &v = o->GetVelocity();
		glm::vec3 pos0 = o->GetLocation();
		glm::vec3 pos1;
		float time = glm::max(0.5f, 4.0f * glm::max(v.x, v.y) / a);
		glm::vec3 directions[8];
		int score[8];
		directions[0] = glm::vec3(1.0f, 0.0f, 0.0f);
		directions[1] = glm::vec3(1.0f, 1.0f, 0.0f);
		directions[2] = glm::vec3(1.0f, -1.0f, 0.0f);
		directions[3] = glm::vec3(0.0f, 1.0f, 0.0f);
		directions[4] = glm::vec3(0.0f, -1.0f, 0.0f);
		directions[5] = glm::vec3(-1.0f, 0.0f, 0.0f);
		directions[6] = glm::vec3(-1.0f, 1.0f, 0.0f);
		directions[7] = glm::vec3(-1.0f, -1.0f, 0.0f);
		int sectorsInTime = time * glm::length(o->GetVelocity()) / sectorSize;
		int intervals = glm::max(20, sectorsInTime * 4);
		//printf("Intervals: %d Time: %f\n", intervals, time);
		float intervalTime = time / intervals;
		float intervalTimeSq = 0.5f * a * intervalTime * intervalTime;
		std::pair<int, int> p0 = GetColRow(pos0);
		const Sector *s0 = GetSector(p0.first, p0.second);
		const OneWaySector *ows = dynamic_cast<const OneWaySector *>(s0);
		if (ows) {
			if (ows->GetType() == ONEWAY_RIGHT) return directions[0];
			if (ows->GetType() == ONEWAY_LEFT) return directions[5];
			if (ows->GetType() == ONEWAY_UP) return directions[3];
			if (ows->GetType() == ONEWAY_DOWN) return directions[4];
		}
		int distance = -10;
		if (s0 && distanceMap.find(s0) != distanceMap.end()) distance = distanceMap.find(s0)->second;
		for (int i = 0; i < 8; i++) {
			score[i] = 0;
			const glm::vec3 &direction = directions[i];
			for (int j = 0; j < intervals; j++) {
				pos1 = pos0 + intervalTime * j * o->GetVelocity() + intervalTimeSq * j * j * direction;
				std::pair<int, int> p = GetColRow(pos1);
				const Sector *s = GetSector(p.first, p.second);
				int newDistance = distance;
				glm::vec3 normal(0.0f);
				const glm::vec3 velocity = o->GetVelocity() + j * intervalTime * a * direction;
				if (s) {
					std::map<const Sector *, int>::const_iterator it = distanceMap.find(s);
					if (it != distanceMap.end()) {
						normal = GetCollisionNormal(pos1, velocity, o->GetRadius());
						if (normal == glm::vec3(0.0f)) {
							newDistance = it->second;
						}
					}
				}
				if (normal == glm::vec3(0.0f)) {
					score[i] += distance - newDistance;
					if (dynamic_cast<const OneWaySector *>(s)) {
						score[i] += 10;
						break;
					}
				}
				else {
					if (0.2f > glm::length(o->GetVelocity())) {
						score[i] += (int)(glm::dot(velocity, normal) * 10);
					}
					break;
				}
			}
		}
		int highScore = -1000;
		glm::vec3 result(0.0f);
		std::vector<int> best;
		for (int i = 0; i < 8; i++) {
			if (score[i] > highScore) {
				highScore = score[i];
				best.clear();
				best.push_back(i);
			}
			else if (score[i] == highScore) {
				best.push_back(i);
			}
		}
		if (!best.empty()) {
			int dir = rand() % best.size();
			result = directions[best.at(dir)];
		}

		// If object has not moved at all, randomize acceleration
		if (locMap.find(o) != locMap.end()) {
			if (locMap.find(o)->second == o->GetLocation()) {
				result = directions[rand() % 8];
				dirMap[o] = result;
				timeMap[o] = 0.05f;
			}
		}
		locMap[o] = o->GetLocation();
		return result;
	}

private:
	// Returns the sector from the given column and row.
	// If given column and row are out of bounds, return 0 instead.
	const Sector * GetSector(int col, int row) const {
		if (col < 0 || col >= (int) sectors.size()) return 0;
		if (row < 0 || row >= (int) sectors.at(col).size()) return 0;
		return sectors.at(col).at(row);
	}

	std::map<const Sector *, int> distanceMap;
	std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> neighborMap;
	std::vector<std::vector<const Sector *>> sectors;
	const float sectorSize;
};

class Game {
public:
	Game(const std::vector<std::string> &data) : g(SECTOR_SIZE) {
		g.Initialize(data, objects);
		for (std::vector<Object *>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
			if (!(*it)->GetAI()) {
				human = *it;
				break;
			}
		}
	}

	~Game() {
		for (std::vector<Object *>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
			delete *it;
		}
	}

	void Run() {
		float deltaTime = GetDeltaTime();

		for (std::vector<Object *>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
			(*it)->Move(deltaTime);
		}

		// First check all collisions to solid grid elements
		for (std::vector<Object *>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
			glm::vec3 n = g.GetCollisionNormal(*it);
			if (n != glm::vec3(0.0f)) {
				// Collision
				(*it)->ResolveCollision(n);
			}
		}

		// Then check collisions between objects
		// Revert movement for any colliding pair and start over again until no collisions are found.
		bool collision;
		do {
			collision = false;
			for (std::vector<Object *>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
				for (std::vector<Object *>::const_iterator it2 = it + 1; it2 != objects.end(); ++it2) {
					Object *o1 = *it;
					Object *o2 = *it2;
					if (o1->Collides(o2)) {
						collision = true;
						o1->ResolveCollision(o2);
						break;
					}
				}
				if (collision) break;
			}
		} while (collision);

		GetKeyboardInputs(human, deltaTime);

		g.Debug(human);

		for (std::vector<Object *>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
			if (*it == human) {
				continue;
			}
			glm::vec3 direction = g.GetAcceleration(*it, deltaTime);
			(*it)->Accelerate(deltaTime, direction);
		}
	}

	void Draw(GLuint id) const {
		glm::mat4 camera = glm::translate(cameraOffset - human->GetLocation());
		g.Draw(id, camera);
		for (std::vector<Object *>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
			(*it)->Draw(id, camera);
		}
	}
private:
	void GetKeyboardInputs(Object *o, float deltaTime) {
		glm::vec3 direction(0.0f, 0.0f, 0.0f);
		if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS) {
			direction.y += 1.0f;
		}
		if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS) {
			direction.y -= 1.0f;
		}
		if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS) {
			direction.x += 1.0f;
		}
		if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS) {
			direction.x -= 1.0f;
		}
		cameraOffset -= deltaTime * direction * 0.3f;
		cameraOffset *= (1.0f - deltaTime);
		o->Accelerate(deltaTime, direction);
	}

	float GetDeltaTime() {
		// glfwGetTime is called only once, the first time this function is called
		static double lastTime = glfwGetTime();

		// Compute time difference between current and last frame
		double currentTime = glfwGetTime();
		float deltaTime = float(currentTime - lastTime);

		// For the next frame, the "last time" will be "now"
		lastTime = currentTime;
		return deltaTime;
	}

	Grid g;
	Object *human;
	std::vector<Object *> objects;
	glm::vec3 cameraOffset;
	std::map<const Object *, std::set<const Sector *>> cache;
	std::map<const Object *, double> cacheTimer;
};

int main( void )
{
	// Initialise GLFW
	if( !glfwInit() )
	{
		fprintf( stderr, "Failed to initialize GLFW\n" );
		getchar();
		return -1;
	}

	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_RESIZABLE,GL_FALSE);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	// Open a window and create its OpenGL context
	window = glfwCreateWindow( 700, 700, "Playground", NULL, NULL);
	if( window == NULL ){
		fprintf( stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n" );
		getchar();
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);

	// Initialize GLEW
	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
		getchar();
		glfwTerminate();
		return -1;
	}

	// Ensure we can capture the escape key being pressed below
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);

	// Dark blue background
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	GLuint VertexArrayID;
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);

	// Create and compile our GLSL program from the shaders
	GLuint programID = LoadShaders("SimpleTransform.vertexshader", "SingleColor.fragmentshader");
	// Get a handle for our "MVP" uniform
	GLuint MatrixID = glGetUniformLocation(programID, "MVP");

	std::vector<std::string> data = ReadGridFromFile("viljo.txt");
	Game game(data);

	srand(time(NULL));

	do {
		// Clear the screen
		glClear(GL_COLOR_BUFFER_BIT);
		// Use our shader
		glUseProgram(programID);

		game.Run();
		game.Draw(MatrixID);

		// Swap buffers
		glfwSwapBuffers(window);
		glfwPollEvents();

	} // Check if the ESC key was pressed or the window was closed
	while( glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
		   glfwWindowShouldClose(window) == 0 );

	// Cleanup VBO
	glDeleteProgram(programID);
	glDeleteVertexArrays(1, &VertexArrayID);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();

	return 0;
}

