#include <iostream>
#include <cstdio>
#include <cstdint>

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

            float color[3];
            color[0] = u;
            color[1] = v;
            color[2] = 0.25f;

            uint8_t rgb8[3];
            rgb8[0] = static_cast<int>(color[0] * 255.999);
            rgb8[1] = static_cast<int>(color[1] * 255.999);
            rgb8[2] = static_cast<int>(color[2] * 255.999);

            fwrite(rgb8, 3, 1, fp);
        }
    }
    std::cout << "\n";
}

