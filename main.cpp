#include <iostream>
#include <cstdio>
#include <cstdint>

#include "vectormath.h"

void write_pixel(FILE *fp, const vec3f& color)
{
    uint8_t rgb8[3];
    for(int i = 0; i < 3; i++) {
        rgb8[i] = static_cast<int>(color[i] * 255.999);
    }
    fwrite(rgb8, 3, 1, fp);
}

int main(int argc, char **argv)
{
    int imageWidth = 512;
    int imageHeight = 512;

    FILE *fp = fopen("image.ppm", "wb");
    fprintf(fp, "P6 %d %d 255\n", imageWidth, imageHeight);

    for(int j = imageHeight - 1; j >= 0; j--) {
        for(int i = 0; i < imageWidth; i++) {
            std::cout << "\rScanlines remaining: " << j << ' ' << std::flush;
            float u = i * 1.0f / (imageWidth - 1);
            float v = j * 1.0f / (imageHeight - 1);

            vec3f color(u, v, .25f);
            write_pixel(fp, color);
        }
    }
    std::cout << "\n";
}
