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

#pragma once

#include <nori/object.h>

NORI_NAMESPACE_BEGIN

struct EmitterQueryRecord {
    Point3f hitPos;
    Point3f srcPos;
    Normal3f n;
    Vector3f wi;
    float pdf;
    Ray3f shadowRay;

    EmitterQueryRecord(const Point3f& hitPos, const Point3f& srcPos, const Normal3f& n) : hitPos(hitPos), srcPos(srcPos), n(n) {
        wi = (srcPos - hitPos).normalized();
    }
};

/**
 * \brief Superclass of all emitters
 */
class Emitter : public NoriObject {
public:

    /**
     * \brief Return the type of object (i.e. Mesh/Emitter/etc.) 
     * provided by this instance
     * */
    EClassType getClassType() const { return EEmitter; }

    virtual ~Emitter() {}

    virtual Color3f eval(const EmitterQueryRecord& record) const = 0;

    virtual Color3f getRadiance() const = 0;

    virtual float pdf(const Mesh* mesh, const EmitterQueryRecord& eRec) const = 0;

    virtual Color3f sample(const Mesh* mesh, EmitterQueryRecord& eRec, Sampler*) const = 0;
};

NORI_NAMESPACE_END
