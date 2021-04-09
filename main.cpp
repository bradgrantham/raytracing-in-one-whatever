#include <iostream>
#include <cstdio>
#include <cstdint>

#include "vectormath.h"

struct shadepoint
{
    vec3f p;
    vec3f n;
    float t;
};

struct hittable
{
    virtual bool hit(const ray& r, float tmin, float tmax, shadepoint *hit) const = 0;
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
};

void writePixel(FILE *fp, const vec3f& color)
{
    uint8_t rgb8[3];
    for(int i = 0; i < 3; i++) {
        rgb8[i] = static_cast<int>(color[i] * 255.999);
    }
    fwrite(rgb8, 3, 1, fp);
}

vec3f evaluateRay(const ray& r, hittable *thingie)
{
    shadepoint point;
    bool hit = thingie->hit(r, -FLT_MAX, FLT_MAX, &point);
    if(hit) {
        return 0.5 * vec3f(point.n.x + 1, point.n.y + 1, point.n.z + 1);
    }

    vec3f dir = vec_normalize(r.m_direction);
    float t = 0.5f * (dir.y + 1.0f);
    return (1.0f - t) * vec3f(1.0f, 1.0f, 1.0f) + t * vec3f(0.5f, 0.7f, 1.0f);
}

int main(int argc, char **argv)
{
    int aspectRatioNum = 16;
    int aspectRatioDenom = 9;
    int imageWidth = 512;
    int imageHeight = imageWidth * aspectRatioDenom / aspectRatioNum;


    FILE *fp = fopen("image.ppm", "wb");
    fprintf(fp, "P6 %d %d 255\n", imageWidth, imageHeight);

    float viewportHeight = 2.0f;
    float viewportWidth = viewportHeight * aspectRatioNum / aspectRatioDenom;
    auto focalLength = 1.0;

    vec3f origin(0, 0, 0);
    vec3f horizontal(viewportWidth, 0, 0);
    vec3f vertical(0, viewportHeight, 0);
    auto lowerLeft = origin - horizontal / 2.0f - vertical / 2.0f - vec3f(0, 0, focalLength);

    sphere scene({0,0,-1}, 0.5);

    for(int j = imageHeight - 1; j >= 0; j--) {
        for(int i = 0; i < imageWidth; i++) {
            std::cout << "\rScanlines remaining: " << j << ' ' << std::flush;

            float u = i * 1.0f / (imageWidth - 1);
            float v = j * 1.0f / (imageHeight - 1);

            ray r(origin, lowerLeft + u * horizontal + v * vertical - origin);

            vec3f color = evaluateRay(r, &scene);

            writePixel(fp, color);
        }
    }
    std::cout << "\n";
}
