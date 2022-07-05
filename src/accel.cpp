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

NORI_NAMESPACE_BEGIN

void Accel::addMesh(Mesh *mesh) {
    if (m_mesh)
        throw NoriException("Accel: only a single mesh is supported!");
    m_mesh = mesh;
    m_bbox = m_mesh->getBoundingBox();
}

struct Accel::OctreeNode
{
    bool isLeaf = false;
    std::vector<Point3i> triangles;
    OctreeNode* child[8];
    OctreeNode()
    {
        for (int i = 0; i < 8;i++)
        {
            child[i] = nullptr;
        }
    }
    OctreeNode(std::vector<Point3i> t)
    {
        isLeaf = true;
        triangles = t;
    }
};

Accel::OctreeNode *Accel::build(TBoundingBox<Point3f> box, std::vector<Point3i> triangles, int depth) {
    if (triangles.size() == 0)
        return nullptr;

    if (triangles.size() < 10 || depth > 2.0 / 3 * log2(m_mesh->getTriangleCount()))
        return new OctreeNode(triangles);

    std::vector<Point3i> list[8];
    MatrixXf vertices = m_mesh->getVertexPositions();

    for (int triIndex = 0; triIndex < triangles.size(); triIndex++) {
        for (int i = 0; i < 8; ++i) {
            Point3i triangle = triangles[triIndex];
            Point3f triMinPoint = vertices.col(triangle[0]), triMaxPoint = vertices.col(triangle[0]);
            for (int j = 1; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    if (triMinPoint[k] > vertices.col(triangle[j])[k])
                    {
                        triMinPoint[k] = vertices.col(triangle[j])[k];
                    }
                    else if (triMaxPoint[k] > vertices.col(triangle[j])[k])
                    {
                        triMaxPoint[k] = vertices.col(triangle[j])[k];
                    }
                }
            }
            Point3f boxMinPoint((i % 2 == 0) ? box.min[0] : box.getCenter()[0], (i % 4 < 2) ? box.min[1] : box.getCenter()[1], (i / 4 == 0) ? box.min[2] : box.getCenter()[2]);
            Point3f boxMaxPoint((i % 2 == 1) ? box.max[0] : box.getCenter()[0], (i % 4 >= 2) ? box.max[1] : box.getCenter()[1], (i / 4 == 1) ? box.max[2] : box.getCenter()[2]);
            if (triMinPoint[0] >= boxMinPoint[0] && triMinPoint[1] >= boxMinPoint[1] && triMinPoint[2] >= boxMinPoint[2]
             && triMaxPoint[0] <= boxMaxPoint[0] && triMaxPoint[1] <= boxMaxPoint[1] && triMaxPoint[2] <= boxMaxPoint[2])
            {
                list[i].push_back(triangle);
            }
        }
    }

    OctreeNode *node = new OctreeNode();
    for (int i = 0; i < 8; ++i)
    {
        Point3f boxMinPoint((i % 2 == 0) ? box.min[0] : box.getCenter()[0], (i % 4 < 2) ? box.min[1] : box.getCenter()[1], (i / 4 == 0) ? box.min[2] : box.getCenter()[2]);
        Point3f boxMaxPoint((i % 2 == 1) ? box.max[0] : box.getCenter()[0], (i % 4 >= 2) ? box.max[1] : box.getCenter()[1], (i / 4 == 1) ? box.max[2] : box.getCenter()[2]);
        node->child[i] = build(TBoundingBox<Point3f>(boxMinPoint, boxMaxPoint), list[i], depth + 1);
    }
    return node;
}

void Accel::build() {
    std::vector<Point3i> list;
    for (int i = 0; i < m_mesh->getTriangleCount(); i++)
    {
        list.push_back(m_mesh->getIndices().col(i));
    }
    build(m_bbox, list, 0);
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

    /* Brute force search through all triangles */
    for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
        float u, v, t;
        if (m_mesh->rayIntersect(idx, ray, u, v, t)) {
            /* An intersection was found! Can terminate
               immediately if this is a shadow ray query */
            if (shadowRay)
                return true;
            ray.maxt = its.t = t;
            its.uv = Point2f(u, v);
            its.mesh = m_mesh;
            f = idx;
            foundIntersection = true;
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

