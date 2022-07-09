#include <nori/emitter.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class AreaLight : public Emitter
{
private:

    Color3f radiance;

public:

    AreaLight(const PropertyList &props)
    {
        radiance = props.getColor("radiance");
    }

    Color3f eval(const EmitterQueryRecord &record) const override
    {
        return (record.n.dot(record.wi) < 0) ? radiance : 0;
    }

    Color3f getRadiance() const override
    {
        return radiance;
    }

    Color3f sample(const Mesh *mesh, EmitterQueryRecord &eRec, Sampler *sampler) const override
    {
        SampleResult rec = mesh->sample(sampler);
        eRec.srcPos = rec.p;
        eRec.n = rec.n;
        eRec.wi = (eRec.srcPos - eRec.hitPos).normalized():
        eRec.shadowRay = Ray3f(eRec.hitPos, eRec.wi, Epsilon, (eRec.srcPos - eRec.hitPos).norm() - Epsilon);
        eRec.pdf = pdf(mesh, eRec);
        if (eRec.pdf > 0 && !std::isnan(eRec.pdf) && !std::isinf(eRec.pdf))
        {
            return eval(eRec) / eRec.pdf;
        }
        return Color3f();
    }

    float pdf(const Mesh *mesh, const EmitterQueryRecord &eRec) const override
    {
        float cosTheta = eRec.n.dot(-eRec.wi);
        if (cosTheta > 0)
        {
            return mesh->getPdf().getNormalization() * (eRec.srcPos - eRec.hitPos).squaredNorm() / cosTheta;
        }
        return 0;
    }

    std::string toString() const override
    {
        return "Emitter[]";
    }
};

NORI_REGISTER_CLASS(AreaLight, "area")
NORI_NAMESPACE_END
