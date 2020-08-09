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

#endif