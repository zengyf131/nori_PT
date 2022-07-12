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

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class Microfacet : public BSDF {
public:
    Microfacet(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        /* Albedo of the diffuse base material (a.k.a "kd") */
        m_kd = propList.getColor("kd", Color3f(0.5f));

        /* To ensure energy conservation, we must scale the 
           specular component by 1-kd. 

           While that is not a particularly realistic model of what 
           happens in reality, this will greatly simplify the 
           implementation. Please see the course staff if you're 
           interested in implementing a more realistic version 
           of this BRDF. */
        m_ks = 1 - m_kd.maxCoeff();
    }

    static float D(const Vector3f &wh, float alpha)
    {
        float a = std::exp(-Frame::tanTheta(wh) * Frame::tanTheta(wh) / alpha / alpha);
        float b = M_PI * alpha * alpha * std::pow(Frame::cosTheta(wh), 3);
        return a / b;
    }

    static float G1(const Vector3f &wv, const Vector3f &wh, float alpha)
    {
        if (wv.dot(wh) / Frame::cosTheta(wv) <= 0) return 0;
        float b = 1 / alpha / Frame::tanTheta(wv);
        if (b >= 1.6) return 1;
        return (3.535*b + 2.181*b*b) / (1 + 2.276*b + 2.577*b*b);
    }

    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const {
    	Vector3f wh = (bRec.wi + bRec.wo).normalized();
        float d = D(wh, m_alpha);
        float g = G1(bRec.wi, wh, m_alpha) * G1(bRec.wo, wh, m_alpha);
        float f = fresnel(wh.dot(bRec.wi), m_extIOR, m_intIOR);
        return m_kd / M_PI + m_ks * ((d * f * g) / (4 * Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) * Frame::cosTheta(wh)));
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
    	Vector3f wh = (bRec.wi + bRec.wo).normalized();
        float d = D(wh, m_alpha);
        float j = 1 / (4 * abs(wh.dot(bRec.wo)));
        return m_ks * d * Frame::cosTheta(wh) * j + (1 - m_ks) * Frame::cosTheta(bRec.wo) / M_PI;
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {
    	
        if (Frame::cosTheta(bRec.wi) <= 0) return Color3f(0);
        if (_sample.x() > m_ks) // diffuse
        {
            Point2f sample((_sample.x() - m_ks) / (1 - m_ks), _sample.y());
            bRec.wo = Warp::squareToCosineHemisphere(sample);
        }
        else    // specular
        {
            Point2f sample(_sample.x() / m_ks, _sample.y());
            Vector3f wh = Warp::squareToBeckmann(sample, m_alpha);
            bRec.wo = ((2 * wh.dot(bRec.wi) * wh) - bRec.wi).normalized();
        }
        if (Frame::cosTheta(bRec.wo) < 0) return Color3f(0);
        return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
        
        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        // return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
    }

    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }

    std::string toString() const {
        return tfm::format(
            "Microfacet[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "  kd = %s,\n"
            "  ks = %f\n"
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR,
            m_kd.toString(),
            m_ks
        );
    }
private:
    float m_alpha;
    float m_intIOR, m_extIOR;
    float m_ks;
    Color3f m_kd;
};

NORI_REGISTER_CLASS(Microfacet, "microfacet");
NORI_NAMESPACE_END
