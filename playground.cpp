#include <stdio.h>
#include <stdlib.h>

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
using namespace glm;

// Base class for visible objects
class GLObject {
protected:
	GLObject(const GLObjectData &data) {
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

	void DrawAt(GLuint id, const glm::mat4 &location, const glm::mat4 &camera) const {
		const glm::mat4 result = location * camera;
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
		glUniformMatrix4fv(id, 1, GL_FALSE, &result[0][0]);
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

glm::vec3 GetLocationFromMatrix(const glm::mat4 &location) {
	glm::vec4 loc4 = location * glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
	return glm::vec3(loc4.x, loc4.y, loc4.z);
}

// Circle shaped object
class Object : public GLObject {
public:
	Object(float radius, int segments, int id, const glm::vec3 &location) : GLObject(CreateCircleData(radius, segments, glm::vec3(1.0f, 0.0f, 0.0f))), id(id) {
		this->radius = radius;
		this->location = glm::mat4(1.0f) * glm::translate(location);
		this->velocity = glm::vec3(0.0f);
		this->force = 1.0f;
		this->mass = 1.0f;
		this->acceleration = force / mass;
	}

	void Draw(GLuint id, const glm::mat4 &camera) const {
		DrawAt(id, location, camera);
	}

	void Move(float time) {
		previousLocation = location;
		location = location * glm::translate(time * velocity);
	}

	glm::vec3 GetLocation() const {
		return GetLocationFromMatrix(location);
	}

	glm::vec3 GetVelocity() const {
		return velocity;
	}

	float GetAcceleration() const {
		return acceleration;
	}

	void Accelerate(float time, const glm::vec3 &direction) {
		velocity += time * acceleration * direction;
		if (length(velocity) > 2.0f) {
			velocity = normalize(velocity) * 2.0f;
		}
	}

	bool Collides(const Object *o) const {
		float len = length((o->location - location) * glm::vec4(1.0f));
		return len < radius + o->radius;
	}

	float GetRadius() const { return radius; }

	void ResolveCollision(Object *o) {
		glm::vec4 coll4 = (o->location - location) * glm::vec4(1.0f);
		glm::vec3 collision(coll4.x, coll4.y, coll4.z);
		glm::vec3 ncoll = normalize(collision);

		float u1 = dot(velocity, ncoll);
		float u2 = dot(o->velocity, ncoll);

		float v1 = (u1 * (mass - o->mass) + 2 * o->mass * u2) / (mass + o->mass);
		float v2 = (u2 * (o->mass - mass) + 2 * mass * u1) / (mass + o->mass);

		velocity += (v1 - u1) * ncoll;
		o->velocity += (v2 - u2) * ncoll;

		location = previousLocation;
		o->location = o->previousLocation;
	}

	void ResolveCollision(const glm::vec3 &normal) {
		velocity = velocity - 2.0f * dot(velocity, normal) * normal;
		velocity *= 0.8f;
		location = previousLocation;
	}

	int GetId() const { return id; }
private:
	float radius;
	float acceleration;
	float force;
	float mass;
	glm::mat4 location;
	glm::mat4 previousLocation;
	glm::vec3 velocity;
	const int id;
};

class Sector {
public:
	Sector(const glm::vec3 &position) {
		location = glm::translate(position);
	}

	virtual ~Sector() { }
	virtual void Draw(GLuint id, const glm::mat4 &camera) const { }
	virtual bool IsSolid(const glm::vec3 &velocity) const { return false; }

	glm::vec3 GetLocation() const {
		return GetLocationFromMatrix(location);
	}
protected:
	glm::mat4 location;
};

class SolidSector : public Sector, public GLObject {
public:
	void Draw(GLuint id, const glm::mat4 &camera) const {
		DrawAt(id, location, camera);
	}

	virtual bool IsSolid(const glm::vec3 &velocity) const { return true; }
	virtual const std::vector<glm::vec3> & GetCorners() const = 0;
protected:
	SolidSector(const glm::vec3 &position, const GLObjectData &data) :
		GLObject(data), Sector(position) { }
};

class SquareSector : public SolidSector {
public:
	SquareSector(glm::vec3 &position, float size, const glm::vec3 &color) :
		SolidSector(position, CreateSquareData(size, color)) {
		corners.reserve(4);
		corners.push_back(glm::vec3(position + glm::vec3(-1.0f, 1.0f, 0.0f) * size * 0.5f)); // TL
		corners.push_back(glm::vec3(position + glm::vec3(1.0f, 1.0f, 0.0f) * size * 0.5f)); // TR
		corners.push_back(glm::vec3(position + glm::vec3(1.0f, -1.0f, 0.0f) * size * 0.5f)); // BR
		corners.push_back(glm::vec3(position + glm::vec3(-1.0f, -1.0f, 0.0f) * size * 0.5f)); // BL
	}

	const std::vector<glm::vec3> & GetCorners() const {
		return corners;
	};
private:
	std::vector<glm::vec3> corners;
};

class OneWaySector : public SquareSector {
public:
	OneWaySector(glm::vec3 &position, float size, const glm::vec3 &color, const glm::vec3 &way) :
		SquareSector(position, size, color), way(way) { }

	bool IsSolid(const glm::vec3 &velocity) const {
		return dot(velocity, way) < 0;
	}

	bool IsUp() const { return dot(way, glm::vec3(0.0f, 1.0f, 0.0f)) > 0; }
	bool IsDown() const { return dot(way, glm::vec3(0.0f, -1.0f, 0.0f)) > 0; }
	bool IsLeft() const { return dot(way, glm::vec3(-1.0f, 0.0f, 0.0f)) > 0; }
	bool IsRight() const { return dot(way, glm::vec3(1.0f, 0.0f, 0.0f)) > 0; }
private:
	const glm::vec3 way;
};
class TriangleSector : public SolidSector {
public:
	TriangleSector(const glm::vec3 &position, float size, const glm::vec3 &color, TriangleType type) :
		SolidSector(position, CreateTriangleData(size, color, type)), type(type) {
		corners.reserve(3);
		if (type != BR)
			corners.push_back(glm::vec3(position + glm::vec3(-1.0f, 1.0f, 0.0f) * size * 0.5f)); // TL
		if (type != BL)
			corners.push_back(glm::vec3(position + glm::vec3(1.0f, 1.0f, 0.0f) * size * 0.5f)); // TR
		if (type != TL)
			corners.push_back(glm::vec3(position + glm::vec3(1.0f, -1.0f, 0.0f) * size * 0.5f)); // BR
		if (type != TR)
			corners.push_back(glm::vec3(position + glm::vec3(-1.0f, -1.0f, 0.0f) * size * 0.5f)); // BL
	}

	const std::vector<glm::vec3> & GetCorners() const {
		return corners;
	};

	TriangleType GetType() const { return type; }
private:
	std::vector<glm::vec3> corners;
	const TriangleType type;
};

/*
Grid has several purposes:
1. Make map creation easy
2. Help AI in route calculation
3. Compute collisions against wall efficiently
*/
class Grid {
public:
	void Initialize(const std::vector<std::string> &data, std::vector<Object *> &objects, float sectorSize, float offset) {
		int cols = 0;
		int rows = data.size();
		for (std::vector<std::string>::const_iterator it = data.begin(); it != data.end(); ++it) {
			const std::string &row = *it;
			if (cols < (int) row.length()) {
				cols = row.length();
			}
		}

		this->offset = offset;
		this->sectorSize = sectorSize;
		float x, y;
		sectors.reserve(cols);
		x = -1.0f + offset;
		for (int i = 0; i < cols; i++, x += sectorSize) {
			std::vector<const Sector *> column;
			column.reserve(rows);
			y = 1.0f - offset;
			for (int j = 0; j < rows; j++, y -= sectorSize) {
				glm::vec3 location(x, y, 0.0f);
				const std::string &row = data.at(j);
				float green = min(0.1f + ((float) i) / cols, 0.9f);
				float blue = min(0.1f + ((float) j) / rows, 0.9f);
				glm::vec3 color(0.0f, green, blue);
				const Sector *sector = 0;
				if ((int) row.length() > i) {
					switch (row.at(i)) {
					case '#':
						sector = new SquareSector(location, sectorSize, color);
						break;
					case '1':
					case '2':
					case '3':
					case '4':
					case '5':
					case '6':
					case '7':
					case '8':
					case '9':
					{
						int id = row.at(i) - '1';
						objects.push_back(new Object(0.04f, 20 - id, id, location));
						break;
					}
					case 'A':
						sector = new TriangleSector(location, sectorSize, color, TL);
						break;
					case 'B':
						sector = new TriangleSector(location, sectorSize, color, TR);
						break;
					case 'C':
						sector = new TriangleSector(location, sectorSize, color, BL);
						break;
					case 'D':
						sector = new TriangleSector(location, sectorSize, color, BR);
						break;
					case '>':
						sector = new OneWaySector(location, sectorSize, color * 0.2f, glm::vec3(1.0f, 0.0f, 0.0f));
						break;
					case '<':
						sector = new OneWaySector(location, sectorSize, color * 0.2f, glm::vec3(-1.0f, 0.0f, 0.0f));
						break;
					}
				}
				if (!sector) {
					sector = new Sector(location);
				}
				column.push_back(sector);
			}
			sectors.push_back(column);
		}
		
		// Calculate neighbor sectors for each sector
		for (int i = 0; i < cols; i++) {
			for (int j = 0; j < rows; j++) {
				if (sectors.at(i).at(j)->IsSolid(glm::vec3(0.0f))) {
					continue;
				}
				std::vector<std::pair<int, int>> neighbors;
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
							if (triangle->GetType() == TL) {
								tr++;
								neighbors.push_back(std::make_pair(i, j - 1));
							}
							else if (triangle->GetType() == TR) {
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
							if (triangle->GetType() == TL) {
								bl++;
								neighbors.push_back(std::make_pair(i - 1, j));
							}
							else if (triangle->GetType() == BL) {
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
							if (triangle->GetType() == BL) {
								br++;
								neighbors.push_back(std::make_pair(i, j + 1));
							}
							else if (triangle->GetType() == BR) {
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
							if (triangle->GetType() == TR) {
								br++;
								neighbors.push_back(std::make_pair(i + 1, j));
							}
							else if (triangle->GetType() == BR) {
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
				//printf("%d,%d: %d\n", i, j, neighbors.size());
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
				glm::vec3 edgeNormal = normalize(glm::vec3(-edge.y, edge.x, 0.0f));
				const glm::vec3 v1 = location - c1;
				float distance = dot(v1, edgeNormal);

				// If we didn't happen to pick the right normal, just change the sign.
				if (distance < 0) {
					edgeNormal = -edgeNormal;
					distance = -distance;
				}

				if (distance < minDistance) {
					const glm::vec3 v2 = location - c2;
					// Feasible sector needed for an edge collision. We guarantee this by checking
					// that angle to point x is under 90 degrees from both corners.
					if (dot(edge, v2) > 0 && dot(-edge, v1) > 0) {
						minDistance = distance;
						normal = edgeNormal;
					}
				}
			}

			if (normal == glm::vec3(0.0f)) {
				// Check distance to corners
				for (std::vector<glm::vec3>::const_iterator it2 = corners.begin(); it2 != corners.end(); ++it2) {
					const glm::vec3 candidate = location - *it2;
					const float distance = length(candidate);
					if (distance < minDistance) {
						minDistance = distance;
						normal = normalize(candidate);
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
		float x = 1.0f - offset + sectorSize / 2 + location.x;
		float y = 1.0f - offset + sectorSize / 2 - location.y;
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
					if (ows->IsRight()) {
						worklist.push_back(std::make_pair(i - 1, j));
					}
					else if (ows->IsLeft()) {
						worklist.push_back(std::make_pair(i + 1, j));
					}
					else if (ows->IsUp()) {
						worklist.push_back(std::make_pair(i, j + 1));
					}
					else if (ows->IsDown()) {
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
				//printf("%d,%d: %d - ", p.first, p.second, d);
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
	/*
	std::set<const Sector *> GetGoodSectors(const Object *o) const {
		std::vector<std::pair<int, int>> route;
		std::pair<int, int> p = GetColRow(o->GetLocation());
		std::set<std::pair<int, int>> visited;
		std::vector<std::pair<int, int>> worklist;
		std::vector<std::vector<std::pair<int, int>>> pathlist;
		std::vector<std::set<const Sector *>> goodSectorList;
		visited.insert(p);
		worklist.push_back(p);
		pathlist.push_back(std::vector<std::pair<int, int>>());
		goodSectorList.push_back(std::set<const Sector *>());
		while (!worklist.empty()) {
			const std::pair<int, int> first = *worklist.begin();
			const std::vector<std::pair<int, int>> path = *pathlist.begin();
			const std::set<const Sector *> goodSectors = *goodSectorList.begin();
			worklist.erase(worklist.begin());
			pathlist.erase(pathlist.begin());
			goodSectorList.erase(goodSectorList.begin());
			const Sector *s = GetSector(first.first, first.second);
			if (s && dynamic_cast<const OneWaySector *>(s)) {
				return goodSectors;
			}
			const std::map<std::pair<int, int>, std::vector<std::pair<int, int>>>::const_iterator it = neighborMap.find(first);
			if (it != neighborMap.end()) {
				const std::vector<std::pair<int, int>> &neighbors = it->second;
				for (std::vector<std::pair<int, int>>::const_iterator vit = neighbors.begin();
					vit != neighbors.end(); ++vit) {
					if (visited.insert(*vit).second) {
						worklist.push_back(*vit);
						std::vector<std::pair<int, int>> newPath = path;
						newPath.push_back(*vit);
						pathlist.push_back(newPath);
						std::set<const Sector *> newGoodSectors = goodSectors;
						newGoodSectors.insert(GetSector(vit->first, vit->second));
						goodSectorList.push_back(newGoodSectors);
					}
				}
			}
		}
		return std::set<const Sector *>();
	}*/

	void Debug(const Object *o) const {
		std::pair<int, int> p = GetColRow(o->GetLocation());
		const Sector *s = GetSector(p.first, p.second);
		if (s && distanceMap.find(s) != distanceMap.end()) {
			printf("Distance: %d Velocity: %f\n", distanceMap.find(s)->second, length(o->GetVelocity()));
		}
	}

	glm::vec3 GetAcceleration(const Object *o, float deltaTime) const {
		float a = o->GetAcceleration();
		const glm::vec3 &v = o->GetVelocity();
		glm::vec3 pos0 = o->GetLocation();
		glm::vec3 pos1;
		float time = max(0.5f, max(v.x, v.y) / a);
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
		int sectorsInTime = time * length(o->GetVelocity()) / sectorSize;
		int intervals = max(20, sectorsInTime * 2);
		float intervalTime = time / intervals;
		float intervalTimeSq = 0.5f * a * intervalTime * intervalTime;
		std::pair<int, int> p0 = GetColRow(pos0);
		const Sector *s0 = GetSector(p0.first, p0.second);
		const OneWaySector *ows = dynamic_cast<const OneWaySector *>(s0);
		if (ows) {
			if (ows->IsRight()) return directions[0];
			if (ows->IsLeft()) return directions[5];
			if (ows->IsUp()) return directions[3];
			if (ows->IsDown()) return directions[4];
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
				bool collision = false;
				int newDistance = distance;
				if (s) {
					collision = true;
					std::map<const Sector *, int>::const_iterator it = distanceMap.find(s);
					if (it != distanceMap.end()) {
						glm::vec3 normal = GetCollisionNormal(pos1, o->GetVelocity() + j * intervalTime * a * direction, o->GetRadius());
						if (normal == glm::vec3(0.0f)) {
							collision = false;
							newDistance = it->second;
						}
					}
				}
				if (!collision) {
					score[i] += distance - newDistance;
					if (dynamic_cast<const OneWaySector *>(s)) {
						score[i] += 10;
						break;
					}
				}
				else {
					score[i] -= (int) (length(o->GetVelocity()) * 10);
					break;
				}
			}
		}
		int highScore = -1000;
		glm::vec3 result(0.0f);
		for (int i = 0; i < 8; i++) {
			//printf("%d ", score[i]);
			if (score[i] > highScore) {
				highScore = score[i];
				result = directions[i];
			}
		}
		//printf("Result: %d\n", highScore);
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
	float sectorSize;
	float offset;
};

class Game {
public:
	Game(const std::vector<std::string> &data) {
		g.Initialize(data, objects, 0.1f, 0.1f);
		for (std::vector<Object *>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
			if ((*it)->GetId() == 0) {
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
			/*
			std::map<const Object *, std::set<const Sector *>>::const_iterator cit = cache.find(*it);
			if (cit == cache.end() || cit->second.empty()) {
				cache[*it] = g.GetGoodSectors(*it);
				cacheTimer[*it] = glfwGetTime();
			}
			else {
				double currentTime = glfwGetTime();
				float diffTime = float(currentTime - cacheTimer[*it]);
				if (diffTime > 2.0f) {
					cache[*it] = g.GetGoodSectors(*it);
					cacheTimer[*it] = currentTime;
				}
			}*/
	        //const std::set<const Sector *> &sectors = cache.find(*it)->second;
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

	std::vector<std::string> data = ReadGridFromFile("map.txt");
	Game game(data);

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

