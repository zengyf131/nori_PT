#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class AOIntegrator : public Integrator {
public:
    AOIntegrator(const PropertyList &props) {
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        Frame frame = its.shFrame;
        Point2f samplePoint = sampler->next2D();
        Point3f weightSample = Warp::squareToCosineHemisphere(samplePoint);
        weightSample = frame.toWorld(weightSample).normalized();
        float visibility = 1;
        if (scene->rayIntersect(Ray3f(its.p + weightSample * 0.00001f, weightSample))) {
            visibility = 0;
        }
        return Color3f(float(visibility));
    }

    std::string toString() const {
        return "AOIntegrator[]";
    }

protected:
};

NORI_REGISTER_CLASS(AOIntegrator, "ao");
NORI_NAMESPACE_END