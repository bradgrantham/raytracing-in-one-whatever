#include <iostream>
#include <vector>
#include <memory>
#include <cstdio>
#include <cstdint>
#include <random>

#include "vectormath.h"

float random_float()
{
    static std::uniform_real_distribution<float> distribution(0.0f, 1.0f);
    static std::mt19937 generator;
    return distribution(generator);
}

float random_float(float min, float max)
{
    // Returns a random real in [min,max).
    return min + (max-min)*random_float();
}

vec3f random_vec3f(float min, float max)
{
    return vec3f(random_float(min, max), random_float(min, max), random_float(min, max));
}

vec3f random_vec3f()
{
    return random_vec3f(0.0f, 1.0f);
}

vec3f randomInUnitSphere()
{
    while (true) {
        vec3f p = random_vec3f(-1,1);
        if (vec_length_sq(p) >= 1) continue;
        return p;
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

struct camera
{
    vec3f horizontal;
    vec3f vertical;
    vec3f lowerLeft;

    camera(int aspectRatioNum, int aspectRatioDenom, float viewportHeight, float focalLength)
    {
        float viewportWidth = viewportHeight * aspectRatioNum / aspectRatioDenom;
        horizontal = vec3f(viewportWidth, 0, 0);
        vertical = vec3f(0, viewportHeight, 0);
        lowerLeft = vec3f(0, 0, 0) - horizontal / 2.0f - vertical / 2.0f - vec3f(0, 0, focalLength);
    }

    ray getRay(float u, float v) const
    {
        return ray(vec3f(0, 0, 0), lowerLeft + u * horizontal + v * vertical - vec3f(0, 0, 0));
    }
};

struct shadepoint
{
    vec3f p;
    vec3f n;
    float t;
};

struct hittable
{
    virtual bool hit(const ray& r, float tmin, float tmax, shadepoint *hit) const = 0;
    virtual ~hittable() {};
};

struct group : public hittable
{
    std::vector<std::shared_ptr<hittable>> children;

    group(std::vector<std::shared_ptr<hittable>> children) :
        children(children)
    {
    }
    virtual bool hit(const ray& r, float tmin, float tmax, shadepoint *hit) const
    {
        float success = false;
        hit->t = tmax;
        for(const auto& child: children) {
            shadepoint candidate;
            if(child->hit(r, tmin, hit->t, &candidate)) {
                *hit = candidate;
                success = true;
            }
        }
        return success;
    }

    virtual ~group() {};
};

struct sphere : public hittable
{
    vec3f center;
    float radius;

    sphere() :
        center(vec3f(0, 0, 0)),
        radius(1)
    {}

    sphere(vec3f center, float radius) :
        center(center),
        radius(radius)
    {}

    virtual bool hit(const ray& r, float tmin, float tmax, shadepoint *hit) const
    {
        vec3f oc = r.m_origin - center;
        float a = vec_dot(r.m_direction, r.m_direction);
        float half_b = vec_dot(oc, r.m_direction);
        float c = vec_dot(oc, oc) - radius * radius;
        float discriminant = half_b * half_b - a * c;
        if (discriminant < 0) {
            return false;
        }

        float sqrtd = sqrtf(discriminant);

        // Find the nearest root that lies in the acceptable range.
        float root = (-half_b - sqrtd) / a;
        if (root < tmin || tmax < root) {
            root = (-half_b + sqrtd) / a;
            if (root < tmin || tmax < root)
                return false;
        }

        hit->t = root;
        hit->p = r.at(hit->t);
        hit->n = (hit->p - center) / radius;

        return true;
    }

    virtual ~sphere() {};
};

void writePixel(FILE *fp, const vec3f& color)
{
    uint8_t rgb8[3];
    for(int i = 0; i < 3; i++) {
        rgb8[i] = static_cast<int>(std::clamp(color[i] * 255.999f, 0.0f, 255.999f));
    }
    fwrite(rgb8, 3, 1, fp);
}

vec3f cast(const ray& r, hittable *thingie, int depth)
{
    if(depth < 0) {
        return vec3f(0, 0, 0);
    }

    shadepoint point;
    bool hit = thingie->hit(r, 0.001f, FLT_MAX, &point);
    if(hit) {
        vec3f bounce = randomInHemisphere(point.n);
        vec3f bounceColor = cast(ray(point.p, bounce), thingie, depth - 1);
        return 0.5 * bounceColor;
    }

    vec3f dir = vec_normalize(r.m_direction);
    float t = 0.5f * (dir.y + 1.0f);
    return (1.0f - t) * vec3f(1.0f, 1.0f, 1.0f) + t * vec3f(0.5f, 0.7f, 1.0f);
}

int main(int argc, char **argv)
{
    int maxBounceDepth = 50;
    int sampleCount = 16;
    int aspectRatioNum = 16;
    int aspectRatioDenom = 9;
    float viewportHeight = 2.0f;
    auto focalLength = 1.0;

    int imageWidth = 512;
    int imageHeight = imageWidth * aspectRatioDenom / aspectRatioNum;

    FILE *fp = fopen("image.ppm", "wb");
    fprintf(fp, "P6 %d %d 255\n", imageWidth, imageHeight);


    camera cam(aspectRatioNum, aspectRatioDenom, viewportHeight, focalLength);

    std::vector<std::shared_ptr<hittable>> shapes;

    shapes.push_back(std::make_shared<sphere>(vec3f(0, -100.5f, -1), 100.0f));

    shapes.push_back(std::make_shared<sphere>(vec3f(0, 0, -1), 0.5));

    group scene = group(shapes);

    for(int j = imageHeight - 1; j >= 0; j--) {
        for(int i = 0; i < imageWidth; i++) {
            std::cout << "\rScanlines remaining: " << j << ' ' << std::flush;

            vec3f color(0, 0, 0);
            for(int sampleIndex = 0; sampleIndex < sampleCount; sampleIndex++) {
                float u = (i + random_float()) * 1.0f / (imageWidth - 1);
                float v = (j + random_float()) * 1.0f / (imageHeight - 1);

                vec3f sample = cast(cam.getRay(u, v), &scene, maxBounceDepth);
                color += sample;
            }
            color /= sampleCount;
            color.x = sqrtf(color.x);
            color.y = sqrtf(color.y);
            color.z = sqrtf(color.z);

            writePixel(fp, color);
        }
    }
    std::cout << "\n";
}
