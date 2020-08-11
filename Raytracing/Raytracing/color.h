#ifndef COLOR_H
#define COLOR_H

#include "vec3.h"
#include <iostream>
#include <algorithm>

void write_color(std::ostream& out, color pixel_color) {
    // Write the translated [0,255] value of each color component.
    out << static_cast<int>(256 * std::clamp(pixel_color.x(), 0.0, 0.999)) << ' '
        << static_cast<int>(256 * std::clamp(pixel_color.y(), 0.0, 0.999)) << ' '
        << static_cast<int>(256 * std::clamp(pixel_color.z(), 0.0, 0.999)) << '\n';
}

void write_color(std::ostream& out, color pixel_color, int samples_per_pixel) {
    auto r = pixel_color.x();
    auto g = pixel_color.y();
    auto b = pixel_color.z();

    // Divide the color by the number of samples.
    auto scale = 1.0 / samples_per_pixel;
    r *= scale;
    g *= scale;
    b *= scale;

    // Write the translated [0,255] value of each color component.
    out << static_cast<int>(256 * std::clamp(r, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * std::clamp(g, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * std::clamp(b, 0.0, 0.999)) << '\n';
}

#endif