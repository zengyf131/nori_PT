#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class PathMatsIntegrator : public Integrator {
public:
    PathMatsIntegrator(const PropertyList &props) {}

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Color3f li(0), t(1);
        Ray3f r = ray;
        float prob;
        int depth = 1;
        while (true)
        {
            Intersection its;
            if (!scene->rayIntersect(r, its)) break;
            if (its.mesh->isEmitter())
            {
                EmitterQueryRecord eRec(r.o, its.p, its.shFrame.n);
                li += t * its.mesh->getEmitter()->eval(eRec);
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
        return "PathMatsIntegrator[]";
    }

protected:
};

NORI_REGISTER_CLASS(PathMatsIntegrator, "path_mats");
NORI_NAMESPACE_END