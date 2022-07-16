#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class PathEmsIntegrator : public Integrator {
public:
    PathEmsIntegrator(const PropertyList &props) {}

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Color3f li(0), t(1);
        Ray3f r = ray;
        float prob;
        int depth = 1;
        bool wasDiffuse = false;
        while (true)
        {
            Intersection its;
            if (!scene->rayIntersect(r, its)) break;
            if (its.mesh->isEmitter())
            {
                EmitterQueryRecord eRec(r.o, its.p, its.shFrame.n);
                if (!wasDiffuse)
                    li += t * its.mesh->getEmitter()->eval(eRec);
            }
            if (its.mesh->getBSDF()->isDiffuse())
            {
                const Mesh *lightMesh = scene->getRandomEmitter(sampler);
                EmitterQueryRecord lightRec(its.p);
                Color3f Lr = lightMesh->getEmitter()->sample(lightMesh, lightRec, sampler);
                if (scene->rayIntersect(lightRec.shadowRay))
                {
                    Lr = 0;
                }
                float cosTheta = Frame::cosTheta(its.shFrame.toLocal(lightRec.wi));
                // if (cosTheta < 0) cosTheta = 0;
                BSDFQueryRecord bsdfRec(its.toLocal(-r.d), its.toLocal(lightRec.wi), ESolidAngle);
                li += Lr * t * its.mesh->getBSDF()->eval(bsdfRec) * cosTheta / (1.0 / scene->getEmitters().size());
                wasDiffuse = true;
            }
            else
            {
                wasDiffuse = false;
            }
            if (depth >= 3)
            {
                prob = std::min(t.maxCoeff(), 0.99f);
                if (sampler->next1D() > prob) break;
                t /= prob;
            }

            BSDFQueryRecord bsdfRec(its.shFrame.toLocal(-r.d));
            Color3f f = its.mesh->getBSDF()->sample(bsdfRec, sampler->next2D());
            t *= f;
            r = Ray3f(its.p, its.toWorld(bsdfRec.wo));
            depth++;
        }
        return li;
    }

    std::string toString() const {
        return "PathEmsIntegrator[]";
    }

protected:
};

NORI_REGISTER_CLASS(PathEmsIntegrator, "path_ems");
NORI_NAMESPACE_END