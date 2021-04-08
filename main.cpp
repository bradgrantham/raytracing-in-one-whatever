#include <iostream>
#include <cstdio>
#include <cstdint>

#include "vectormath.h"

void writePixel(FILE *fp, const vec3f& color)
{
    uint8_t rgb8[3];
    for(int i = 0; i < 3; i++) {
        rgb8[i] = static_cast<int>(color[i] * 255.999);
    }
    fwrite(rgb8, 3, 1, fp);
}

bool sphereIsHit(const vec3f& center, float radius, const ray& r)
{
    vec3f oc = r.m_origin - center;
    float a = vec_dot(r.m_direction, r.m_direction);
    float b = 2.0f * vec_dot(oc, r.m_direction);
    float c = vec_dot(oc, oc) - radius * radius;
    float discriminant = b * b - 4 * a * c;
    return (discriminant > 0);
}

vec3f evaluateRay(const ray& r)
{
    if (sphereIsHit(vec3f(0, 0, -1), 0.5, r)) {
        return vec3f(1, 0, 0);
    }
    vec3f dir = vec_normalize(r.m_direction);
    auto t = 0.5f * (dir.y + 1.0f);
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

    for(int j = imageHeight - 1; j >= 0; j--) {
        for(int i = 0; i < imageWidth; i++) {
            std::cout << "\rScanlines remaining: " << j << ' ' << std::flush;

            float u = i * 1.0f / (imageWidth - 1);
            float v = j * 1.0f / (imageHeight - 1);

            ray r(origin, lowerLeft + u * horizontal + v * vertical - origin);

            vec3f color = evaluateRay(r);

            writePixel(fp, color);
        }
    }
    std::cout << "\n";
}
