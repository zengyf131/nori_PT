#include <nori/integrator.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

class SimpleIntegrator : public Integrator {
public:
    SimpleIntegrator(const PropertyList &props) {
        position = props.getPoint("position");
        color = props.getColor("energy");
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        Frame frame = its.shFrame;
        Point3f hitPos = its.p;
        Vector3f lightRayDir = (position - hitPos).normalized();
        float visibility = 1;
        if (scene->rayIntersect(Ray3f(hitPos + lightRayDir * 0.00001f, lightRayDir))) {
            visibility = 0;
        }
        float cosTheta = frame.cosTheta(frame.toLocal(lightRayDir));
        float s = std::pow((hitPos - position).norm(), 2.0f);
        Color3f li = (color / (4 * M_PI * M_PI)) * ((std::max(0.0f, cosTheta)) / s) * visibility;
        return li;
    }

    std::string toString() const {
        return "SimpleIntegrator[]";
    }

protected:
    Point3f position;
    Color3f color;
};

NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_NAMESPACE_END