/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/accel.h>
#include <Eigen/Geometry>
#include <stack>
#include <algorithm>
#include <tbb/parallel_for.h>
#include <nori/timer.h>

NORI_NAMESPACE_BEGIN

void Accel::addMesh(Mesh *mesh) {
    // if (m_mesh)
    //     throw NoriException("Accel: only a single mesh is supported!");
    // m_mesh = mesh;
    // m_bbox = m_mesh->getBoundingBox();
    m_meshes.push_back(mesh);
    m_roots.push_back(nullptr);
    if (m_meshes.size() == 1)
    {
        m_bbox = m_meshes[0]->getBoundingBox();
    }
    else
    {
        m_bbox.expandBy(mesh->getBoundingBox());
    }
}

struct Accel::OctreeNode
{
    bool isLeaf = false;
    std::vector<int> triangles;
    OctreeNode* child[8];
    BoundingBox3f bbox;
    OctreeNode(BoundingBox3f b)
    {
        isLeaf = false;
        for (int i = 0; i < 8;i++)
        {
            child[i] = nullptr;
        }
        bbox = b;
    }
    OctreeNode(BoundingBox3f b, std::vector<int> t)
    {
        isLeaf = true;
        triangles = t;
        bbox = b;
    }
};

struct Accel::NodeWithT
{
    OctreeNode* node;
    float t;
};

Accel::OctreeNode *Accel::build(Mesh *mesh, BoundingBox3f box, std::vector<int> triangles, int depth) {
    if (triangles.size() == 0)
        return nullptr;

    if (triangles.size() < 10 || depth > 2.0 / 3 * log2(mesh->getTriangleCount()))
    {
        // for (int i = 0; i < (int)triangles.size(); i++)
        //     std::cout << triangles[i] << ' ';
        // std::cout << std::endl;
        return new OctreeNode(box, triangles);
    }

    std::vector<int> list[8];
    MatrixXf vertices = mesh->getVertexPositions();
    MatrixXu faces = mesh->getIndices();

    for (int triIndex = 0; triIndex < (int)triangles.size(); triIndex++) {
        int triangleIndex = triangles[triIndex];
        Point3i triangle(faces(0, triangleIndex), faces(1, triangleIndex), faces(2, triangleIndex));
        Point3f triMinPoint = vertices.col(triangle[0]), triMaxPoint = vertices.col(triangle[0]);
        for (int j = 1; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                if (triMinPoint[k] > vertices.col(triangle[j])[k])
                {
                    triMinPoint[k] = vertices.col(triangle[j])[k];
                }
                else if (triMaxPoint[k] < vertices.col(triangle[j])[k])
                {
                    triMaxPoint[k] = vertices.col(triangle[j])[k];
                }
            }
        }
        // std::cout << "start" << std::endl;
        for (int i = 0; i < 8; ++i) {
            Point3f boxMinPoint((i % 2 == 0) ? box.min[0] : box.getCenter()[0], (i % 4 < 2) ? box.min[1] : box.getCenter()[1], (i < 4) ? box.min[2] : box.getCenter()[2]);
            Point3f boxMaxPoint((i % 2 == 1) ? box.max[0] : box.getCenter()[0], (i % 4 >= 2) ? box.max[1] : box.getCenter()[1], (i >= 4) ? box.max[2] : box.getCenter()[2]);
            // std::cout << boxMinPoint[0] << ' ' << boxMinPoint[1] << ' ' << boxMinPoint[2] << std::endl;
            // std::cout << boxMaxPoint[0] << ' ' << boxMaxPoint[1] << ' ' << boxMaxPoint[2] << std::endl;
            if (triMaxPoint[0] >= boxMinPoint[0] && triMaxPoint[1] >= boxMinPoint[1] && triMaxPoint[2] >= boxMinPoint[2]
             && triMinPoint[0] <= boxMaxPoint[0] && triMinPoint[1] <= boxMaxPoint[1] && triMinPoint[2] <= boxMaxPoint[2])
            {
                list[i].push_back(triangleIndex);
            }
        }
    }

    OctreeNode *node = new OctreeNode(box);
    tbb::parallel_for(tbb::blocked_range<int>(0, 8), [mesh, node, list, box, depth, this](const tbb::blocked_range<int> &r)
        {
            for (int i = r.begin(); i != r.end(); i++)
            {
                Point3f boxMinPoint((i % 2 == 0) ? box.min[0] : box.getCenter()[0], (i % 4 < 2) ? box.min[1] : box.getCenter()[1], (i < 4) ? box.min[2] : box.getCenter()[2]);
                Point3f boxMaxPoint((i % 2 == 1) ? box.max[0] : box.getCenter()[0], (i % 4 >= 2) ? box.max[1] : box.getCenter()[1], (i >= 4) ? box.max[2] : box.getCenter()[2]);
                node->child[i] = this->build(mesh, BoundingBox3f(boxMinPoint, boxMaxPoint), list[i], depth + 1);
            }
        }
    );
    // for (int i = 0; i < 8; ++i)
    // {
    //     Point3f boxMinPoint((i % 2 == 0) ? box.min[0] : box.getCenter()[0], (i % 4 < 2) ? box.min[1] : box.getCenter()[1], (i < 4) ? box.min[2] : box.getCenter()[2]);
    //     Point3f boxMaxPoint((i % 2 == 1) ? box.max[0] : box.getCenter()[0], (i % 4 >= 2) ? box.max[1] : box.getCenter()[1], (i >= 4) ? box.max[2] : box.getCenter()[2]);
    //     node->child[i] = build(BoundingBox3f(boxMinPoint, boxMaxPoint), list[i], depth + 1);
    // }
    return node;
}

void Accel::build() {
    std::cout << "Building Octree...";
    Timer timer;
    for (int meshIndex = 0; meshIndex < (int)m_meshes.size(); meshIndex++)
    {
        Mesh *mesh = m_meshes[meshIndex];
        std::vector<int> list;
        for (int i = 0; i < (int)mesh->getTriangleCount(); i++)
        {
            list.push_back(i);
        }
        m_roots[meshIndex] = build(mesh, mesh->getBoundingBox(), list, 0);
    }
    std::cout << "done in " << timer.elapsedString() << std::endl;
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

    float u, v, t;
    auto cmp = [](const NodeWithT a, const NodeWithT b)->bool
    {
        if (a.node == nullptr) return false;
        if (b.node == nullptr) return true;
        return a.t > b.t;
    };
    for (int meshIndex = 0; meshIndex < (int)m_meshes.size(); meshIndex++)
    {
        std::stack<OctreeNode*> callStack;
        callStack.push(m_roots[meshIndex]);

        while (!callStack.empty())
        {
            OctreeNode *node = callStack.top();
            callStack.pop();
            if (!node->bbox.rayIntersect(ray)) continue;
            if (!node->isLeaf)
            {
                NodeWithT hitList[8];
                for (int i = 0; i < 8; i++)
                {
                    OctreeNode *child = node->child[i];
                    hitList[i].node = nullptr;
                    float tmax;
                    if (child && child->bbox.rayIntersect(ray, hitList[i].t, tmax) && hitList[i].t < ray.maxt)
                    {
                        hitList[i].node = child;
                    }
                    // if (child && child->bbox.rayIntersect(ray))
                    // {
                    //     callStack.push(child);
                    // }
                }
                std::sort(hitList, hitList + 8, cmp);
                for (int i = 0; i < 8; i++)
                {
                    if (hitList[i].node == nullptr) break;
                    callStack.push(hitList[i].node);
                }
                continue;
            }
            for (int i = 0; i < (int)node->triangles.size(); i++)
            {
                if (m_meshes[meshIndex]->rayIntersect(node->triangles[i], ray, u, v, t))
                {
                    if (shadowRay)
                        return true;
                    ray.maxt = its.t = t;
                    its.uv = Point2f(u, v);
                    its.mesh = m_meshes[meshIndex];
                    f = node->triangles[i];
                    foundIntersection = true;
                }
            }
        }
    }
    
    if (foundIntersection) {
        /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.

           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates, etc..)
        */

        /* Find the barycentric coordinates */
        Vector3f bary;
        bary << 1-its.uv.sum(), its.uv;

        /* References to all relevant mesh buffers */
        const Mesh *mesh   = its.mesh;
        const MatrixXf &V  = mesh->getVertexPositions();
        const MatrixXf &N  = mesh->getVertexNormals();
        const MatrixXf &UV = mesh->getVertexTexCoords();
        const MatrixXu &F  = mesh->getIndices();

        /* Vertex indices of the triangle */
        uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

        Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

        /* Compute the intersection positon accurately
           using barycentric coordinates */
        its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

        /* Compute proper texture coordinates if provided by the mesh */
        if (UV.size() > 0)
            its.uv = bary.x() * UV.col(idx0) +
                bary.y() * UV.col(idx1) +
                bary.z() * UV.col(idx2);

        /* Compute the geometry frame */
        its.geoFrame = Frame((p1-p0).cross(p2-p0).normalized());

        if (N.size() > 0) {
            /* Compute the shading frame. Note that for simplicity,
               the current implementation doesn't attempt to provide
               tangents that are continuous across the surface. That
               means that this code will need to be modified to be able
               use anisotropic BRDFs, which need tangent continuity */

            its.shFrame = Frame(
                (bary.x() * N.col(idx0) +
                 bary.y() * N.col(idx1) +
                 bary.z() * N.col(idx2)).normalized());
        } else {
            its.shFrame = its.geoFrame;
        }
    }

    return foundIntersection;
}

NORI_NAMESPACE_END

