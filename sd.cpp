#include "sd_gl.h"
#include "sd_data.h"
#include "timer.h"

#include <time.h>
#include <deque>
#include <map>
#include <set>
#include <algorithm>

using namespace std;

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
        health = data.health;
        id = data.id;
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

    int GetHealth() const {
        return health;
    }

    int GetId() const {
        return id;
    }

    // Resolve collision for both objects.
    // Collision is currently assumed to be fully elastic.
    void ResolveCollision(Object * const o) {
        assert(Collides(o));

        glm::vec3 old1 = velocity;
        glm::vec3 old2 = o->velocity;

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

        health -= 2 * elasticity * o->elasticity * glm::length(velocity - old1);
        o->health -= 2 * elasticity * o->elasticity * glm::length(o->velocity - old2);
    }

    // Resolve collision for this object.
    // Collision is not fully elastic.
    void ResolveCollision(const glm::vec3 &normal) {
        glm::vec3 old = velocity;
        velocity = velocity - 2.0f * glm::dot(velocity, normal) * normal;

        health -= 2 * elasticity * glm::length(velocity - old);

        velocity *= elasticity;

        // Prevent movement.
        location = previousLocation;
    }

private:
    float radius;
    float force;
    float mass;
    float acceleration;
    float maxVelocity;
    float elasticity;
    glm::vec3 location;
    glm::vec3 previousLocation;
    glm::vec3 velocity;
    int health;
    int id;
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
    SquareSector(const glm::vec3 &location, float size, const glm::vec3 &color, SectorType type) :
        SolidSector(location, CreateSquareData(size, color, type)) {
        corners.reserve(4);
        corners.push_back(location + UNIT_TL * size);
        corners.push_back(location + UNIT_TR * size);
        corners.push_back(location + UNIT_BR * size);
        corners.push_back(location + UNIT_BL * size);
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
        SquareSector(location, size, color, GetType(way)),
        way(way), type(GetType(way)) {
    }

    bool IsSolid(const glm::vec3 &velocity) const {
        return glm::dot(velocity, way) < 0;
    }

    SectorType GetType() const {
        return type;
    }

private:
    static SectorType GetType(const glm::vec3 &way) {
        if (way == UNIT_T) return ONEWAY_UP;
        else if (way == UNIT_B) return ONEWAY_DOWN;
        else if (way == UNIT_L) return ONEWAY_LEFT;
        else if (way == UNIT_R) return ONEWAY_RIGHT;
        return SQUARE;
    }

    const glm::vec3 way;
    const SectorType type;
};

class TriangleSector : public SolidSector {
public:
    TriangleSector(const glm::vec3 &location, float size, const glm::vec3 &color, SectorType type) :
        SolidSector(location, CreateTriangleData(size, color, type)), type(type) {
        corners.reserve(3);
        if (type != TRIANGLE_BR)
            corners.push_back(location + UNIT_TL * size);
        if (type != TRIANGLE_BL)
            corners.push_back(location + UNIT_TR * size);
        if (type != TRIANGLE_TL)
            corners.push_back(location + UNIT_BR * size);
        if (type != TRIANGLE_TR)
            corners.push_back(location + UNIT_BL * size);
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

class ObjectAI {
public:
    virtual ~ObjectAI() { }
    virtual void Accelerate(float deltaTime) = 0;
};

class AI {
public:

    virtual ~AI() {
        for (map<const Object *, ObjectAI *>::const_iterator it = objects.begin();
             it != objects.end(); ++it) {
            delete it->second;
        }
    }

    virtual void Initialize(const map<const Sector *, vector<const Sector *> > &neighborMap,
                            const map<const OneWaySector *, const Sector *> &previousSectorMap) { }

    virtual void AddObject(Object *o) = 0;

    virtual int GetDistance(const Sector *s) const = 0;

    void Accelerate(float deltaTime) {
        for (map<const Object *, ObjectAI *>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
            it->second->Accelerate(deltaTime);
        }
    }

    Object * GetHuman(const vector<Object *> &allObjects) const {
        for (vector<Object *>::const_iterator it = allObjects.begin(); it != allObjects.end(); ++it) {
            if (objects.find(*it) == objects.end()) {
                return *it;
            }
        }
        return 0;
    }

protected:

    void RegisterAI(const Object *o, ObjectAI *ai) {
        objects[o] = ai;
    }

private:

    map<const Object *, ObjectAI *> objects;
};

const Sector* CreateSector(char c, vector<Object *> &objects,
                           const glm::vec3 &location, float halfSectorSize, const glm::vec3 &color, AI *ai) {
    if (c >= '1' && c <= '9') {
        int id = c - '1';
        objects.push_back(new Object(CreateObjectData(id), location));
        if (c != '1') {
            ai->AddObject(objects.back());
        }
    }
    const Sector *sector = 0;
    switch (c) {
    case '#':
        sector = new SquareSector(location, halfSectorSize, color, SQUARE);
        break;
    case 'A':
        sector = new TriangleSector(location, halfSectorSize, color, TRIANGLE_TL);
        break;
    case 'B':
        sector = new TriangleSector(location, halfSectorSize, color, TRIANGLE_TR);
        break;
    case 'C':
        sector = new TriangleSector(location, halfSectorSize, color, TRIANGLE_BL);
        break;
    case 'D':
        sector = new TriangleSector(location, halfSectorSize, color, TRIANGLE_BR);
        break;
    case '>':
        sector = new OneWaySector(location, halfSectorSize, color * 0.5f, UNIT_R);
        break;
    case '<':
        sector = new OneWaySector(location, halfSectorSize, color * 0.5f, UNIT_L);
        break;
    case '^':
        sector = new OneWaySector(location, halfSectorSize, color * 0.5f, UNIT_T);
        break;
    case 'v':
        sector = new OneWaySector(location, halfSectorSize, color * 0.5f, UNIT_B);
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
    Grid(float sectorSize) : sectorSize(sectorSize), halfSectorSize(sectorSize * 0.5f) { }

    void Initialize(const vector<string> &data, vector<Object *> &objects, AI *ai) {

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
                column.push_back(CreateSector(c, objects, location, halfSectorSize, color, ai));
            }
            sectors.push_back(column);
        }

        // Finally, initialize and feed some maps to AI...
        map<const OneWaySector *, const Sector *> previousSectorMap;
        map<const Sector *, vector<const Sector *> > neighborMap;

        // Calculate neighbor sectors for each sector
        for (unsigned int i = 0; i < cols; i++) {
            for (unsigned int j = 0; j < rows; j++) {

                // Solid sectors have no neighbors
                const Sector * const s = sectors.at(i).at(j);
                if (s->IsSolid(ZERO)) {
                    continue;
                }

                vector<const Sector *> &neighbors =
                    neighborMap.insert(make_pair(s, vector<const Sector *>())).first->second;

                // tl, tr, bl and br are counters for diagonal corners and the counter gets
                // an increment if any of its neighbors is accessible. If counter is incremented
                // twice it means that the corresponding diagonal corner is also a neighbor.
                int tl = 0;
                int tr = 0;
                int bl = 0;
                int br = 0;

                const Sector *sector;

                // top neighbor
                if (j > 0) {
                    sector = sectors.at(i).at(j - 1);
                    if (!sector->IsSolid(UNIT_T)) {
                        neighbors.push_back(sector);
                        tl++;
                        tr++;
                        if (sector->GetType() == ONEWAY_UP) {
                            previousSectorMap[static_cast<const OneWaySector *>(sector)] = s;
                        }
                    }
                    else if (sector->GetType() == TRIANGLE_TL) {
                        tr++;
                        neighbors.push_back(sector);
                    }
                    else if (sector->GetType() == TRIANGLE_TR) {
                        tl++;
                        neighbors.push_back(sector);
                    }
                }

                // left neighbor
                if (i > 0) {
                    sector = sectors.at(i - 1).at(j);
                    if (!sector->IsSolid(UNIT_L)) {
                        neighbors.push_back(sector);
                        tl++;
                        bl++;
                        if (sector->GetType() == ONEWAY_LEFT) {
                            previousSectorMap[static_cast<const OneWaySector *>(sector)] = s;
                        }
                    }
                    else if (sector->GetType() == TRIANGLE_TL) {
                        bl++;
                        neighbors.push_back(sector);
                    }
                    else if (sector->GetType() == TRIANGLE_BL) {
                        tl++;
                        neighbors.push_back(sector);
                    }
                }

                // bottom neighbor
                if (j < sectors.at(i).size() - 1) {
                    sector = sectors.at(i).at(j + 1);
                    if (!sector->IsSolid(UNIT_B)) {
                        neighbors.push_back(sector);
                        bl++;
                        br++;
                        if (sector->GetType() == ONEWAY_DOWN) {
                            previousSectorMap[static_cast<const OneWaySector *>(sector)] = s;
                        }
                    }
                    else if (sector->GetType() == TRIANGLE_BL) {
                        br++;
                        neighbors.push_back(sector);
                    }
                    else if (sector->GetType() == TRIANGLE_BR) {
                        bl++;
                        neighbors.push_back(sector);
                    }
                }

                // right neighbor
                if (i < sectors.size() - 1) {
                    sector = sectors.at(i + 1).at(j);
                    if (!sector->IsSolid(UNIT_R)) {
                        neighbors.push_back(sector);
                        br++;
                        tr++;
                        if (sector->GetType() == ONEWAY_RIGHT) {
                            previousSectorMap[static_cast<const OneWaySector *>(sector)] = s;
                        }
                    }
                    else if (sector->GetType() == TRIANGLE_TR) {
                        br++;
                        neighbors.push_back(sector);
                    }
                    else if (sector->GetType() == TRIANGLE_BR) {
                        tr++;
                        neighbors.push_back(sector);
                    }
                }

                // diagonal neighbors
                if (tr == 2) {
                    sector = sectors.at(i + 1).at(j - 1);
                    if (!sector->IsSolid(UNIT_TR)) {
                        neighbors.push_back(sector);
                    }
                }
                if (tl == 2) {
                    sector = sectors.at(i - 1).at(j - 1);
                    if (!sector->IsSolid(UNIT_TL)) {
                        neighbors.push_back(sector);
                    }
                }
                if (br == 2) {
                    sector = sectors.at(i + 1).at(j + 1);
                    if (!sector->IsSolid(UNIT_BR)) {
                        neighbors.push_back(sector);
                    }
                }
                if (bl == 2) {
                    sector = sectors.at(i - 1).at(j + 1);
                    if (!sector->IsSolid(UNIT_BL)) {
                        neighbors.push_back(sector);
                    }
                }
            }
        }

        ai->Initialize(neighborMap, previousSectorMap);
    }

    ~Grid() {
        for (unsigned int i = 0; i < sectors.size(); i++) {
            for (unsigned int j = 0; j < sectors.at(i).size(); j++) {
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
        glm::vec3 normal = ZERO;
        const pair<int, int> p = GetColRow(location);
        const Sector *sector = GetSector(p.first, p.second);
        if (!sector) {
            return normal;
        }
        vector<const SolidSector *> collisionCandidates;
        AddCollisionCandidates(location, velocity, radius, sector, p, collisionCandidates);
        float minDistance = radius;
        for (vector<const SolidSector *>::const_iterator it = collisionCandidates.begin();
             it != collisionCandidates.end(); ++it) {
            const vector<glm::vec3> &corners = (*it)->GetCorners();

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
        }

        // No edge collision, because they were too far away or sector was not feasible.
        if (normal == ZERO) {
            for (vector<const SolidSector *>::const_iterator it = collisionCandidates.begin();
                 it != collisionCandidates.end(); ++it) {
                const vector<glm::vec3> &corners = (*it)->GetCorners();

                // Check distance to corners
                for (vector<glm::vec3>::const_iterator it2 = corners.begin(); it2 != corners.end(); ++it2) {
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

    void Debug(const Object *o) const {
    }

    const Sector * GetSector(const glm::vec3 &location) const {
        const pair<int, int> p = GetColRow(location);
        return GetSector(p.first, p.second);
    }

    float GetSectorSize() const {
        return sectorSize;
    }

private:
    // Returns a pair (col, row) which corresponds to the given location
    pair<int, int> GetColRow(const glm::vec3 &location) const {
        float x = halfSectorSize + location.x;
        float y = halfSectorSize - location.y;
        int col = (int) (x / sectorSize);
        int row = (int) (y / sectorSize);
        return make_pair(col, row);
    }

    void AddCollisionCandidates(const glm::vec3 &location, const glm::vec3 &velocity, const float radius,
                                const Sector *sector, const pair<int, int> &p, vector<const SolidSector *> &collisionCandidates) const {
        const glm::vec3 diff = sector->GetLocation() - location;
        int colDelta = 0;
        int rowDelta = 0;
        if (glm::abs(diff.x) + radius > halfSectorSize) {
            colDelta = diff.x > 0 ? -1 : 1;
        }
        if (glm::abs(diff.y) + radius > halfSectorSize) {
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

    // Returns the sector from the given column and row.
    // If given column and row are out of bounds, return 0 instead.
    const Sector * GetSector(int col, int row) const {
        if (col < 0 || col >= (int) sectors.size()) return 0;
        if (row < 0 || row >= (int) sectors.at(col).size()) return 0;
        return sectors.at(col).at(row);
    }

    vector<vector<const Sector *> > sectors;
    const float sectorSize;
    const float halfSectorSize;
};

// Object specific AI
class DummyObjectAI : public ObjectAI {
public:
    DummyObjectAI(Object *o, const AI *ai, const Grid &grid) : o(o), ai(ai), grid(grid) {
        previousLocation = ZERO;
        randomMovementTimeLeft = 0.0f;
        minLookaheadTime = 0.5f;
        lookaheadMultiplier = 4.0f;
        minIntervalCount = 20;
        intervalMultiplier = 4;
        unknownDistance = -10;
        oneWaySectorBonus = 10;
        lowSpeedThreshold = 0.2f;
        collisionMultiplier = 10;
        defaultHighScore = -1000;
        randomMovementTime = 0.05f;
    }

    void Accelerate(float deltaTime) {
        o->Accelerate(deltaTime, GetAcceleration(deltaTime));
    }

private:

    // AI parameters
    float minLookaheadTime;
    float lookaheadMultiplier;
    int minIntervalCount;
    int intervalMultiplier;
    int unknownDistance;
    int oneWaySectorBonus;
    float lowSpeedThreshold;
    int collisionMultiplier;
    int defaultHighScore;
    float randomMovementTime;

    glm::vec3 GetAcceleration(float deltaTime) {

        // Check if random movement is active
        if (randomMovementTimeLeft > 0.0f) {
            randomMovementTimeLeft -= deltaTime;
            return randomDirection;
        }

        float a = o->GetAcceleration();
        const glm::vec3 &v = o->GetVelocity();
        const glm::vec3 pos0 = o->GetLocation();
        glm::vec3 pos1;
        float time = glm::max(minLookaheadTime, lookaheadMultiplier * glm::max(v.x, v.y) / a);
        static const int DIRECTIONS = 8;
        static glm::vec3 directions[] =
            { UNIT_T, UNIT_B, UNIT_L, UNIT_R, UNIT_TL, UNIT_TR, UNIT_BL, UNIT_BR };
        int score[DIRECTIONS];
        int sectorsInTime = (int) (time * glm::length(o->GetVelocity()) / grid.GetSectorSize());
        int intervals = glm::max(minIntervalCount, sectorsInTime * intervalMultiplier);
        float intervalTime = time / intervals;
        float intervalTimeSq = 0.5f * a * intervalTime * intervalTime;

        const Sector *s0 = grid.GetSector(pos0);
        if (s0) {
            switch (s0->GetType()) {
            case ONEWAY_RIGHT: return UNIT_R;
            case ONEWAY_LEFT: return UNIT_L;
            case ONEWAY_UP: return UNIT_T;
            case ONEWAY_DOWN: return UNIT_B;
            }
        } else {
            return ZERO;
        }

        int distance = ai->GetDistance(s0);
        distance = distance < 0 ? unknownDistance : distance;

        // Calculate a score for acceleration to each possible direction.
        for (int i = 0; i < DIRECTIONS; i++) {
            score[i] = 0;
            const glm::vec3 &direction = directions[i];
            for (int j = 0; j < intervals; j++) {
                pos1 = pos0 + intervalTime * j * o->GetVelocity() + intervalTimeSq * j * j * direction;
                const Sector *s = grid.GetSector(pos1);
                int newDistance = distance;
                glm::vec3 normal = ZERO;
                const glm::vec3 velocity = o->GetVelocity() + j * intervalTime * a * direction;
                int dist = ai->GetDistance(s);
                if (dist >= 0) {
                    normal = grid.GetCollisionNormal(pos1, velocity, o->GetRadius());
                    if (normal == ZERO) {
                        newDistance = dist;
                    }
                }
                if (normal == ZERO) {
                    score[i] += distance - newDistance;
                    if (dynamic_cast<const OneWaySector *>(s)) {
                        score[i] += oneWaySectorBonus;
                        break;
                    }
                }
                else {
                    // If object moves very slowly, do not penalize for collisions.
                    // Otherwise it's difficult to move in narrow maze.
                    if (lowSpeedThreshold > glm::length(o->GetVelocity())) {
                        score[i] += (int) (glm::dot(velocity, normal) * collisionMultiplier);
                    }
                    break;
                }
            }
        }

        // Randomize direction from those which have the best score
        int highScore = defaultHighScore;
        glm::vec3 result = ZERO;
        vector<int> best;
        for (int i = 0; i < DIRECTIONS; i++) {
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
        if (previousLocation == o->GetLocation()) {
            randomDirection = directions[rand() % DIRECTIONS];
            randomMovementTimeLeft = randomMovementTime;
            result = randomDirection;
        }

        previousLocation = o->GetLocation();
        return result;
    }

    const AI * const ai;
    Object * const o;
    const Grid &grid;

    glm::vec3 previousLocation;
    glm::vec3 randomDirection;
    float randomMovementTimeLeft;
};

// Master AI
class DummyAI : public AI {
public:
    DummyAI(const Grid &g) : grid(g) { }

    // Calculates shortest distances from all accessible sectors to the closest
    // OneWaySector and stores the result to a map
    void Initialize(const map<const Sector *, vector<const Sector *> > &neighborMap,
                    const map<const OneWaySector *, const Sector *> &previousSectorMap) {

        deque<const Sector *> worklist;
        for (map<const OneWaySector *, const Sector *>::const_iterator it = previousSectorMap.begin();
             it != previousSectorMap.end(); ++it) {
            distanceMap[it->first] = 0;
            worklist.push_back(it->second);
        }

        // Distance from OneWaySector
        int currentDistance = 1;

        // Separator for increasing distance
        static const Sector * const sep = 0;
        worklist.push_back(sep);

        while (!worklist.empty()) {
            const Sector *sector = worklist.front();
            worklist.pop_front();
            if (sector == sep) {
                currentDistance++;
                if (!worklist.empty()) {
                    worklist.push_back(sep);
                }
                continue;
            }
            if (distanceMap.find(sector) != distanceMap.end()) {
                // Already covered earlier, possibly with less distance
                continue;
            }
            else {
                distanceMap[sector] = currentDistance;
                const map<const Sector *, vector<const Sector *> >::const_iterator it =
                    neighborMap.find(sector);
                if (it != neighborMap.end()) {
                    const vector<const Sector *> &neighbors = it->second;
                    for (vector<const Sector *>::const_iterator vit = neighbors.begin();
                         vit != neighbors.end(); ++vit) {
                        worklist.push_back(*vit);
                    }
                }
            }
        }
    }

    void AddObject(Object *o) {
        RegisterAI(o, new DummyObjectAI(o, this, grid));
    }

    int GetDistance(const Sector *sector) const {
        int result = -1;
        map<const Sector *, int>::const_iterator it = distanceMap.find(sector);
        if (it != distanceMap.end()) {
            result = it->second;
        }
        return result;
    }

private:

    map<const Sector *, int> distanceMap;
    const Grid &grid;
};

// Sorts objects based on the number of laps completed.
// If number of laps completed is equal, then sort based on
// distance from the goal. Unknown distance is considered to be
// the largest. Finally, if all lap count is equal and all laps
// have been completed, it means that the distance would also be
// equal and then we compare the time spent.
class ObjectComparator {
public:
    ObjectComparator(const map<const Object *, int> &deathDistance, const map<const Object *, vector<float> > &lapTimes, unsigned int laps) : dd(deathDistance), lt(lapTimes), laps(laps) { }
    bool operator() (const Object *a, const Object *b) const {
        int laps1 = GetLaps(a);
        int laps2 = GetLaps(b);
        if (laps1 != laps2) {
            // Different number of laps completed
            return laps1 > laps2;
        } else if (laps1 == laps) {
            // All laps completed
            return lt.find(a)->second.back() < lt.find(b)->second.back();
        }
        // Same number (but not all) laps completed.
        map<const Object *, int>::const_iterator mit1 = dd.find(a);
        map<const Object *, int>::const_iterator mit2 = dd.find(b);
        int dist1 = mit1->second;
        int dist2 = mit2->second;
        if (dist1 == -1) {
            return false;
        } else if (dist2 == -1) {
            return true;
        }
        return dist1 < dist2;
    }

    int GetLaps(const Object *o) const {
        map<const Object *, vector<float> >::const_iterator it = lt.find(o);
        return it == lt.end() ? 0 : it->second.size();
    }
private:
    const map<const Object *, int> &dd;
    const map<const Object *, vector<float> > &lt;
    const unsigned int laps;
};

uint64_t Hash(const char *s, uint64_t seed = 0) {
    unsigned int hash = seed;
    while (*s) {
        hash = hash * 101 + *s++;
    }
    return hash;
}

class Game {
public:
    Game(const vector<string> &data, unsigned int laps, GLFWwindow * const window, const TextRenderer &tr) : g(SECTOR_SIZE), laps(laps), window(window), textRenderer(tr) {
        ai = new DummyAI(g);
        g.Initialize(data, objects, ai);
        human = ai->GetHuman(objects);
        running = false;
        time = 0.0f;
        hash = laps;
        for (vector<string>::const_iterator it = data.begin(); it != data.end(); ++it) {
            hash = Hash(it->c_str(), hash);
        }
    }

    ~Game() {
        for (vector<Object *>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
            delete *it;
        }
        for (vector<Object *>::const_iterator it = deadObjects.begin(); it != deadObjects.end(); ++it) {
            delete *it;
        }
        delete ai;
    }

    void End() {
        DropObjects();
        Standings();
    }

    bool Run() {
        float deltaTime = GetDeltaTime();
        if (deltaTime > 0.4f) {
            // Huge lag for some reason
            return true;
        }

        set<const Object *> toNewLap;

        for (vector<Object *>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
            const Sector *s = g.GetSector((*it)->GetLocation());
            if (s && dynamic_cast<const OneWaySector *>(s)) {
                //glm::vec3 loc = (*it)->GetLocation();
                //printf("%f,%f\n", loc.x, loc.y);
                toNewLap.insert(*it);
            }
            (*it)->Move(deltaTime);
        }

        // First check all collisions to solid grid elements
        for (vector<Object *>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
            glm::vec3 n = g.GetCollisionNormal(*it);
            if (n != ZERO) {
                // Collision
                (*it)->ResolveCollision(n);
            }
        }

        // Then check collisions between objects
        // Revert movement for any colliding pair and start over again until no collisions are found.
        bool collision;
        do {
            collision = false;
            for (vector<Object *>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
                for (vector<Object *>::const_iterator it2 = it + 1; it2 != objects.end(); ++it2) {
                    Object *o1 = *it;
                    Object *o2 = *it2;
                    if (o1->Collides(o2)) {
                        collision = true;
                        o1->ResolveCollision(o2);
                        break;
                    }
                }
                if (collision) {
                    break;
                }
            }
        } while (collision);

        // Check death condition
        for (vector<Object *>::iterator it = objects.begin(); it != objects.end(); ) {
            if ((*it)->GetHealth() <= 0) {
                it = RemoveObject(it);
                if (objects.empty()) {
                    // Game over
                    return false;
                }
            } else {
                ++it;
            }
        }

        // Check if any of the objects which are about to start a new lap actually do it
        for (set<const Object *>::const_iterator it = toNewLap.begin(); it != toNewLap.end(); ++it) {
            const Sector *s = g.GetSector((*it)->GetLocation());
            if (s && !dynamic_cast<const OneWaySector *>(s)) {
                vector<float> &times = lapTimes[*it];
                times.push_back(float(glfwGetTime() - initialTime));
                if (times.size() >= laps) {
                    if (*it == human) {
                        time = times.back();
                    }
                    for (vector<Object *>::iterator vit = objects.begin(); vit != objects.end(); ) {
                        if (*it == *vit) {
                            RemoveObject(vit);
                            if (objects.empty()) {
                                // Game over
                                return false;
                            }
                        } else {
                            ++vit;
                        }
                    }
                }
            }
        }

        if (human) {
            GetKeyboardInputs(human, deltaTime);
            g.Debug(human);
        }

        if (running) {
            ai->Accelerate(deltaTime);
        }
        return true;
    }

    void Draw(GLuint id) const {
        Object *cameraTarget = human ? human : *objects.begin();
        glm::mat4 camera = glm::translate(cameraOffset - cameraTarget->GetLocation());
        g.Draw(id, camera);
        for (vector<Object *>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
            (*it)->Draw(id, camera);
        }

        // Render health
        float x = 25.0;
        stringstream ss;
        for (vector<Object *>::const_iterator it = objects.begin(); it != objects.end(); ++it) {
            ss.str("");
            ss << (*it)->GetHealth();
            textRenderer.RenderText(ss.str(), x, 25.0f, 0.5f, (*it)->GetColor());
            x += 50.0;
        }
    }

    bool IsRunning() const {
        return running;
    }

private:
    void DropObjects() {
        for (vector<Object *>::iterator it = objects.begin(); it != objects.end(); ) {
            it = RemoveObject(it);
        }
    }

    void Standings() const {
        float y = 500.0f;
        glClear(GL_COLOR_BUFFER_BIT);
        textRenderer.RenderText("STANDINGS", 150.0f, y, 1.0f, glm::vec3(1.0f, 1.0f, 1.0f));
        y -= 60.0f;
        int i = 1;
        ObjectComparator cmp(deathDistance, lapTimes, laps);
        vector<const Object *> standings(deadObjects.begin(), deadObjects.end());
        sort(standings.begin(), standings.end(), cmp);
        stringstream ss;
        for (vector<const Object *>::const_iterator it = standings.begin(); it != standings.end(); ++it) {
            int completedLaps = cmp.GetLaps(*it);
            ss << (i++) << ". " << completedLaps << " lap(s), ";
            if (laps == completedLaps) {
                ss << "time: " << lapTimes.find(*it)->second.back();
            } else {
                ss << "dist: " << deathDistance.find(*it)->second;
            }
            textRenderer.RenderText(ss.str(), 150.0f, y, 0.5f, (*it)->GetColor());
            y -= 30.0f;
            ss.str("");
        }
        if (time > 0.0f) {
            map<uint64_t, float> highscores = ReadRecords("highscores.dat");
            map<uint64_t, float>::const_iterator it = highscores.find(hash);
            if (it == highscores.end()) {
                highscores[hash] = time;
                textRenderer.RenderText("New record!", 150.0f, y, 0.3f, glm::vec3(0.8f, 0.8f, 0.8f));
                WriteRecords("highscores.dat", highscores);
            } else if (it->second > time) {
                highscores[hash] = time;
                textRenderer.RenderText("New record!", 150.0f, y, 0.3f, glm::vec3(0.8f, 0.8f, 0.8f));
                WriteRecords("highscores.dat", highscores);
            }
        }
        glBindVertexArray(0);
        glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_FALSE);
        do {
            glfwSwapBuffers(window);
            glfwPollEvents();
        } while (glfwGetKey(window, GLFW_KEY_ENTER) != GLFW_PRESS);
        glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
    }

    vector<Object *>::iterator RemoveObject(vector<Object *>::iterator it) {
        Object *o = *it;
        if (human == o) {
            human = 0;
        }
        deadObjects.push_back(o);
        const Sector *sector = g.GetSector(o->GetLocation());
        int distance = ai->GetDistance(sector);
        deathDistance[o] = distance;
        return objects.erase(it);
    }

    void GetKeyboardInputs(Object *o, float deltaTime) {
        glm::vec3 direction = ZERO;
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

        if (!running) {
            running = direction != ZERO;
            if (running) {
                initialTime = glfwGetTime();
            }
        }
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
    AI *ai;
    Object *human;
    vector<Object *> objects;
    vector<Object *> deadObjects;
    map<const Object *, vector<float> > lapTimes;
    map<const Object *, int> deathDistance;
    glm::vec3 cameraOffset;
    bool running;
    double initialTime;
    const unsigned int laps;
    GLFWwindow * const window;
    const TextRenderer &textRenderer;
    uint64_t hash;
    float time;
};

int Program::Run() const {

    if (!initialized) {
        return -1;
    }

    srand((unsigned int) time(NULL));
    deque<pair<string, unsigned int> > season = ReadSeason("maps/season1.ssn");
    for (deque<pair<string, unsigned int> >::const_iterator it =
             season.begin(); it != season.end(); ++it) {
        vector<string> map = ReadGridFromFile(it->first.c_str());
        const unsigned int laps = it->second;
        Game game(map, laps, window, textRenderer);
        do {
            glClear(GL_COLOR_BUFFER_BIT);
            shader.Use();
            glBindVertexArray(VertexArrayID);
            if (!game.Run()) {
                // over
                break;
            }

            game.Draw(MatrixID);
            glBindVertexArray(0);

            glfwSwapBuffers(window);
            glfwPollEvents();

            if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS) {
                return 0;
            }
            if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS ||
                glfwWindowShouldClose(window) != 0) {
                if (game.IsRunning()) {
                    break;
                }
            }
        } while (true);
        game.End();
    }

    return 0;
}

int main() {
    Program p(700, 700, "Solar Dash");
    return p.Run();
}
