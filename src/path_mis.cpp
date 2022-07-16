#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class PathMisIntegrator : public Integrator {
public:
    PathMisIntegrator(const PropertyList &props) {}

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Color3f li(0), t(1);
        Ray3f r = ray;
        float prob;
        int depth = 1;
        float w_brdf = 1;
        Intersection its;
        if (!scene->rayIntersect(r, its)) return li;
        while (true)
        {
            if (its.mesh->isEmitter())
            {
                EmitterQueryRecord eRec(r.o, its.p, its.shFrame.n);
                li += t * its.mesh->getEmitter()->eval(eRec) * w_brdf;
            }
            const Mesh *lightMesh = scene->getRandomEmitter(sampler);
            EmitterQueryRecord lightRec(its.p);
            Color3f Lr = lightMesh->getEmitter()->sample(lightMesh, lightRec, sampler) / (1.0 / scene->getEmitters().size());
            float pdf_emitter = lightMesh->getEmitter()->pdf(lightMesh, lightRec);
            if (!scene->rayIntersect(lightRec.shadowRay))
            {
                float cosTheta = Frame::cosTheta(its.shFrame.toLocal(lightRec.wi));
                if (cosTheta < 0) cosTheta = 0;
                BSDFQueryRecord bsdfRec(its.toLocal(-r.d), its.toLocal(lightRec.wi), ESolidAngle);
                float pdf_brdf = its.mesh->getBSDF()->pdf(bsdfRec);
                float w_emitter = (pdf_brdf + pdf_emitter > 0.0f) ? pdf_emitter / (pdf_brdf + pdf_emitter) : pdf_emitter;
                li += Lr * t * its.mesh->getBSDF()->eval(bsdfRec) * cosTheta * w_emitter;
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
            float pdf_brdf = its.mesh->getBSDF()->pdf(bsdfRec);
            Point3f origin = its.p;
            if (!scene->rayIntersect(r, its)) break;
            if (its.mesh->isEmitter()) {
                EmitterQueryRecord newLightRec = EmitterQueryRecord(origin, its.p, its.shFrame.n);
                float new_pdf_emitter = its.mesh->getEmitter()->pdf(its.mesh, newLightRec);
                w_brdf = pdf_brdf + new_pdf_emitter > 0.f ? pdf_brdf / (pdf_brdf + new_pdf_emitter) : pdf_brdf;
            }
            if (bsdfRec.measure == EDiscrete) {
                w_brdf = 1.0f;
            }
            depth++;
        }
        return li;
    }

    std::string toString() const {
        return "PathMisIntegrator[]";
    }

protected:
};

NORI_REGISTER_CLASS(PathMisIntegrator, "path_mis");
NORI_NAMESPACE_END