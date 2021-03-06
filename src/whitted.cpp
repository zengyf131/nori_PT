#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class WhittedIntegrator : public Integrator {
public:
    WhittedIntegrator(const PropertyList &props) {}

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const override {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        Color3f Le(0), Lr(0);
        if (its.mesh->isEmitter())
        {
            EmitterQueryRecord rec(ray.o, its.p, its.shFrame.n);
            Le = its.mesh->getEmitter()->eval(rec);
        }
        if (its.mesh->getBSDF()->isDiffuse())
        {
            const Mesh *lightMesh = scene->getRandomEmitter(sampler);
            EmitterQueryRecord lightRec(its.p);
            Lr = lightMesh->getEmitter()->sample(lightMesh, lightRec, sampler);
            if (scene->rayIntersect(lightRec.shadowRay))
            {
                Lr = 0;
            }
            float cosTheta = Frame::cosTheta(its.shFrame.toLocal(lightRec.wi));
            if (cosTheta < 0) cosTheta = 0;
            BSDFQueryRecord bsdfRec(its.toLocal(-ray.d), its.toLocal(lightRec.wi), ESolidAngle);
            return Le + Lr * its.mesh->getBSDF()->eval(bsdfRec) * cosTheta / (1.0 / scene->getEmitters().size());
        }
        else
        {
            BSDFQueryRecord bsdfRec(its.toLocal(-ray.d));
            Color3f refColor = its.mesh->getBSDF()->sample(bsdfRec, sampler->next2D());
            if (sampler->next1D() < 0.95 && refColor.x() > 0)
            {
                return Li(scene, sampler, Ray3f(its.p, its.toWorld(bsdfRec.wo))) / 0.95 * refColor;
            }
            else
            {
                return Color3f(0);
            }
        }
    }

    std::string toString() const {
        return "WhittedIntegrator[]";
    }

protected:
};

NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");
NORI_NAMESPACE_END