#include <iostream>
#include <vector>
#include <memory>
#include <random>
#include <thread>
#include <atomic>
#include <cstdio>
#include <cmath>
#include <cstdint>
#include <future>

#include "vectormath.h"

//--------------------------------------------------------------------------
// Random number utility functions

float randomFloat()
{
    static std::uniform_real_distribution<float> distribution(0.0f, 1.0f);
    static std::mt19937 generator;
    return distribution(generator);
}

float randomFloat(float min, float max)
{
    // Returns a random real in [min,max).
    return min + (max - min) * randomFloat();
}

vec3f randomInUnitDisk()
{
    while (true) {
        vec3f p = vec3f(randomFloat(-1,1), randomFloat(-1,1), 0);
        if (vec_length_sq(p) < 1) {
            return p;
        }
    }
}

vec3f randomVec3f(float min, float max)
{
    return vec3f(randomFloat(min, max), randomFloat(min, max), randomFloat(min, max));
}

vec3f randomVec3f()
{
    return randomVec3f(0.0f, 1.0f);
}

vec3f randomInUnitSphere()
{
    while (true) {
        vec3f p = randomVec3f(-1,1);
        if (vec_length_sq(p) <= 1.0f) {
            return p;
        }
    }
}

vec3f randomUnitVec3f()
{
    return vec_normalize(randomInUnitSphere());
}

vec3f randomInHemisphere(const vec3f& normal)
{
    vec3f v = randomInUnitSphere();
    if (vec_dot(v, normal) > 0.0f) // In the same hemisphere as the normal
        return v;
    else
        return -v;
}

//--------------------------------------------------------------------------
// Camera

struct Camera
{
    vec3f origin;
    vec3f horizontal;
    vec3f vertical;
    vec3f lowerLeft;
    float lensRadius;
    vec3f w, u, v;

    Camera(const vec3f& eye, const vec3f& target, const vec3f& up,
        int aspectRatioNum, int aspectRatioDenom, float vfov, float aperture,
        float focusDistance)
    {
        float theta = vfov / 180.0f * 3.14159f; 
        float h = tanf(theta / 2.0f);
        float viewportHeight = 2.0f * h;
        float viewportWidth = viewportHeight * aspectRatioNum / aspectRatioDenom;

        w = vec_normalize(eye - target);
        u = vec_normalize(vec_cross(up, w));
        v = vec_cross(w, u);

        origin = eye;
        horizontal = focusDistance * viewportWidth * u;
        vertical = focusDistance * viewportHeight * v;
        lowerLeft = origin - horizontal / 2.0f - vertical / 2.0f - focusDistance * w;
        lensRadius = aperture / 2.0f;
    }

    ray getRay(float u_, float v_) const
    {
        vec3f rd = lensRadius * randomInUnitDisk();
        vec3f offset = u * rd.x + v * rd.y;
        return ray(origin + offset, lowerLeft + u_ * horizontal + v_ * vertical - origin - offset);
    }
};


//--------------------------------------------------------------------------
// Material, materials, and shading

struct Material
{
    virtual bool scatter(const ray& incident, const vec3f& p, const vec3f& n, bool front, vec3f& attenuation, vec3f& outgoingDirection) const = 0;
    typedef std::shared_ptr<Material> Ptr;
    virtual ~Material() {}
};

bool isNearZero(const vec3f& v)
{
    // Return true if the vector is close to zero in all dimensions.
    constexpr auto s = 1e-8;
    return (fabsf(v[0]) < s) && (fabsf(v[1]) < s) && (fabsf(v[2]) < s);
}

struct Diffuse : public Material
{
    vec3f color;
    Diffuse(const vec3f& color) :
        color(color)
    {}
    virtual bool scatter(const ray& incident, const vec3f& p, const vec3f& n, bool front, vec3f& attenuation, vec3f& outgoingDirection) const
    {
        outgoingDirection = n + randomUnitVec3f();
        if (vec_dot(outgoingDirection, n) < 0.0f) {
            outgoingDirection = -outgoingDirection;
        }
        if(isNearZero(outgoingDirection)) {
            outgoingDirection = n;
        }
        attenuation = color;
        return true;
    }
    virtual ~Diffuse() {}
};

float reflectance(float cosine, float ref_idx)
{
    // Use Schlick's approximation for reflectance.
    auto r0 = (1-ref_idx) / (1+ref_idx);
    r0 = r0*r0;
    return r0 + (1-r0)*powf((1 - cosine),5);
}

struct Dielectric : public Material
{
    float index;
    Dielectric(float index) :
        index(index)
    {}
    virtual bool scatter(const ray& incident, const vec3f& p, const vec3f& n, bool front, vec3f& attenuation, vec3f& outgoingDirection) const
    {
        vec3f i = vec_normalize(incident.m_direction);
        float eta = front ? (1.0f / index) : index;
        bool internal = vec_refract(eta, i, n, outgoingDirection);
        float cos_theta = std::min(vec_dot(-i, n), 1.0f);
        if(internal || reflectance(cos_theta, eta) > randomFloat()) {
            outgoingDirection = vec_reflect(i, n);
        }
        attenuation = vec3f(1, 1, 1);
        return true;
    }
    virtual ~Dielectric() {}
};

struct Metal : public Material
{
    vec3f color;
    float gloss;
    Metal(const vec3f& color, float gloss) :
        color(color),
        gloss(std::clamp(gloss, 0.0f, 1.0f))
    {}
    virtual bool scatter(const ray& incident, const vec3f& p, const vec3f& n, bool front, vec3f& attenuation, vec3f& outgoingDirection) const
    {
        vec3f i = vec_normalize(incident.m_direction);
        vec3f r = vec_reflect(i, n);
        outgoingDirection = vec_reflect(i, n) + gloss * randomInUnitSphere();
        if(vec_dot(outgoingDirection, n) < 0) {
            outgoingDirection = r;
        }
        attenuation = color;
        return vec_dot(outgoingDirection, n) > 0;
    }
    virtual ~Metal() {}
};

struct ShadeParams
{
    vec3f p;
    vec3f n;
    bool f;
    float t;
    Material *m;        // Control flow that stores a Material* here must guarantee a Material would not be deleted before the ShadeParams is dtor'd

    void setNormal(const vec3f& incident, const vec3f& n_)
    {
        f = (vec_dot(incident, n_) < 0);
        n = f ? n_ : -n_;
    }
};

vec3f shadeBackground(const ray &r)
{
    vec3f dir = vec_normalize(r.m_direction);
    float t = 0.5f * (dir.y + 1.0f);
    return (1.0f - t) * vec3f(1.0f, 1.0f, 1.0f) + t * vec3f(0.5f, 0.7f, 1.0f);
}

//--------------------------------------------------------------------------
// Things that can be hit with a ray

struct Hittable
{
    virtual bool hit(const ray& r, float tmin, float tmax, ShadeParams *hit) const = 0;
    virtual ~Hittable() {};
    typedef std::shared_ptr<Hittable> Ptr;
};

struct Group : public Hittable
{
    std::vector<Hittable::Ptr> children;

    Group(std::vector<Hittable::Ptr> children) :
        children(children)
    {
    }
    virtual bool hit(const ray& r, float tmin, float tmax, ShadeParams *hit) const
    {
        float success = false;

        for(const auto& child: children) {
            ShadeParams candidate;
            if(child->hit(r, tmin, tmax, &candidate)) {
                *hit = candidate;
                tmax = candidate.t;
                success = true;
            }
        }

        return success;
    }

    virtual ~Group() {};
};

struct Sphere : public Hittable
{
    vec3f center;
    float radius;
    Material::Ptr mtl;

    Sphere() :
        center(vec3f(0, 0, 0)),
        radius(1)
    {}

    Sphere(vec3f center, float radius, Material::Ptr mtl) :
        center(center),
        radius(radius),
        mtl(mtl)
    {}

    virtual bool hit(const ray& r, float tmin, float tmax, ShadeParams *hit) const
    {
        vec3f oc = r.m_origin - center;
        float a = vec_dot(r.m_direction, r.m_direction);
        float half_b = vec_dot(oc, r.m_direction);
        float c = vec_dot(oc, oc) - radius * radius;
        float discriminant = half_b * half_b - a * c;
        if (discriminant < 0) {
            return false;
        }

        float sd = sqrtf(discriminant);

        // Find the nearest root that lies in the acceptable range.
        float root = (-half_b - sd) / a;
        if (root < tmin || tmax < root) {
            root = (-half_b + sd) / a;
            if (root < tmin || tmax < root) {
                return false;
            }
        }

        vec3f point = r.at(root);
        vec3f normal = (point - center) / radius;

        hit->t = root;
        hit->p = point;
        hit->m = mtl.get();
        hit->setNormal(r.m_direction, normal);

        return true;
    }

    virtual ~Sphere() {};
};

//--------------------------------------------------------------------------
// Image I/O

void writePixel(FILE *fp, const vec3f& color)
{
    uint8_t rgb8[3];
    for(int i = 0; i < 3; i++) {
        rgb8[i] = static_cast<int>(std::clamp(color[i] * 255.999f, 0.0f, 255.999f));
    }
    fwrite(rgb8, 3, 1, fp);
}

//--------------------------------------------------------------------------
// Cast a ray and calculate a color

/*
 * Materials must NOT be deleted from a ShadeParams by any function
 * called from cast()
 */
vec3f cast(const ray& r, Hittable::Ptr thingie, int depth)
{
    if(depth < 0) {
        return vec3f(0, 0, 0);
    }

    ShadeParams params;
    bool hit = thingie->hit(r, 0.001f, FLT_MAX, &params);
    if(hit) {
        vec3f bounceColor;
        vec3f bounce;
        bool again = params.m->scatter(r, params.p, params.n, params.f, bounceColor, bounce);
        if(again) {
            return bounceColor * cast(ray(params.p, bounce), thingie, depth - 1);
        }
        return vec3f(0, 0, 0);
    }

    return shadeBackground(r);
}

//--------------------------------------------------------------------------

Hittable::Ptr generateRandomScene()
{
    std::vector<Hittable::Ptr> shapes;

    auto ground_material = std::make_shared<Diffuse>(vec3f(0.5, 0.5, 0.5));
    shapes.push_back(std::make_shared<Sphere>(vec3f(0.0f,-1000.0f,0.0f), 1000.0f, ground_material));

    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            auto choose_mat = randomFloat();
            vec3f center(a + 0.9f*randomFloat(), 0.2f, b + 0.9f*randomFloat());

            if (vec_length(center - vec3f(4, 0.2f, 0)) > 0.9f) {
                std::shared_ptr<Material> sphere_material;

                if (choose_mat < 0.8f) {
                    // diffuse
                    auto albedo = randomVec3f() * randomVec3f();
                    sphere_material = std::make_shared<Diffuse>(albedo);
                    shapes.push_back(std::make_shared<Sphere>(center, 0.2f, sphere_material));
                } else if (choose_mat < 0.95f) {
                    // metal
                    auto albedo = randomVec3f(0.5, 1);
                    auto fuzz = randomFloat(0, 0.5);
                    sphere_material = std::make_shared<Metal>(albedo, fuzz);
                    shapes.push_back(std::make_shared<Sphere>(center, 0.2f, sphere_material));
                } else {
                    // glass
                    sphere_material = std::make_shared<Dielectric>(1.5f);
                    shapes.push_back(std::make_shared<Sphere>(center, 0.2f, sphere_material));
                }
            }
        }
    }

    auto material1 = std::make_shared<Dielectric>(1.5f);
    shapes.push_back(std::make_shared<Sphere>(vec3f(0, 1, 0), 1.0f, material1));

    auto material2 = std::make_shared<Diffuse>(vec3f(0.4f, 0.2f, 0.1f));
    shapes.push_back(std::make_shared<Sphere>(vec3f(-4, 1, 0), 1.0f, material2));

    auto material3 = std::make_shared<Metal>(vec3f(0.7f, 0.6f, 0.5f), 0.0f);
    shapes.push_back(std::make_shared<Sphere>(vec3f(4, 1, 0), 1.0f, material3));

    return std::make_shared<Group>(shapes);
}

int main(int argc, char **argv)
{
    int maxBounceDepth = 50;
    int sampleCount = 500;
    int aspectRatioNum = 3;
    int aspectRatioDenom = 2;
    int imageWidth = 1200;
    int imageHeight = imageWidth * aspectRatioDenom / aspectRatioNum;
    vec3f *image = new vec3f[imageWidth * imageHeight];
    float *imageWeight = new float[imageWidth * imageHeight];

    vec3f eye(13, 2, 3);
    vec3f target(0,0,0);
    vec3f up(0,1,0);
    float distToFocus = 10.0;
    float aperture = .1f;
    Camera cam(eye, target, up, aspectRatioNum, aspectRatioDenom, 20, aperture, distToFocus);

    Hittable::Ptr scene = generateRandomScene();

    for(int j = imageHeight - 1; j >= 0; j--) {
        for(int i = 0; i < imageWidth; i++) {
            image[i + j * imageWidth] = vec3f(0, 0, 0);
            imageWeight[i + j * imageWidth] = 0;
        }
    }

    std::atomic_int rowsRemaining = imageHeight;

    std::vector<std::future<void>> lines;
    for(int j = 0; j < imageHeight; j++) {
        auto f = [&,j]{
            for(int i = 0; i < imageWidth; i++) {

                vec3f &pixel = image[i + j * imageWidth];
                float &weight = imageWeight[i + j * imageWidth];

                for(int sampleIndex = 0; sampleIndex < sampleCount; sampleIndex++) {
                    float u = (i + randomFloat()) / (imageWidth - 1);
                    float v = (j + randomFloat()) / (imageHeight - 1);

                    pixel += cast(cam.getRay(u, v), scene, maxBounceDepth);
                    weight += 1;
                }
            }
            rowsRemaining --;
        };
        auto line = std::async(std::launch::async, f);
        lines.push_back(std::move(line));
    }

    do {
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
        std::cout << "\rScanlines remaining: " << rowsRemaining << "  " << std::flush;
    } while(rowsRemaining > 0);

    for(auto& l: lines) {
        l.wait();
    }

    for(int j = imageHeight - 1; j >= 0; j--) {
        for(int i = 0; i < imageWidth; i++) {
            vec3f &pixel = image[i + j * imageWidth];
            float &weight = imageWeight[i + j * imageWidth];
            pixel /= weight;
            pixel.x = sqrtf(pixel.x);
            pixel.y = sqrtf(pixel.y);
            pixel.z = sqrtf(pixel.z);
        }
    }

    FILE *fp = fopen("image.ppm", "wb");
    fprintf(fp, "P6 %d %d 255\n", imageWidth, imageHeight);
    for(int j = imageHeight - 1; j >= 0; j--) {
        for(int i = 0; i < imageWidth; i++) {
            vec3f &pixel = image[i + j * imageWidth];
            writePixel(fp, pixel);
        }
    }
    fclose(fp);

    std::cout << "\n";
}
