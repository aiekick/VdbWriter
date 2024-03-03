#include "VDBWriter.hpp"
#include <cstdint>   // int32_t
#include <chrono>    // std::chrono
#include <iostream>  // std::cout

#define JULIA_REVOLUTE
//#define USE_ANIMATED_WAVE
//#define USE_VDB_WRITER

#ifdef JULIA_REVOLUTE
static double mix(const double& x, const double& y, const double& a) { return x * (1.0 - a) + y * a; }
int           main() {
    const int32_t SIZE       = 375;
    const double  ZOOM_XZ    = 7.5;
    const double  ZOOM_Y     = 7.5;
    const int32_t ITERATIONS = 5;
    const int32_t FRAMES     = 10;

    double time = 0.0f;
    double kk, hh, px, pz, an, cx, cy, path, rev_x, rev_y, tmp_x, tmp_y, rev_x_squared, rev_y_squared, df;
    double rot2D[4] = {1, 0, 0, 1};  // c,-s,s,c for t=0

    std::array<int32_t, 3> offset       = {};
    bool                   first_offset = true;
    vdb::VDBWriter         vdb;
    int32_t cube_color;
    double  time_step = 6.28318 / (double)FRAMES;
    for (int32_t f = 0; f < FRAMES; ++f) {
        vdb.setKeyFrame(f);
        auto* densityLayerPtr = vdb.getFloatLayer(0, "density");
        auto* sdfLayerPtr     = vdb.getDoubleLayer(1, "sdf");
        time += time_step;
        for (int32_t i = -SIZE; i < SIZE; ++i) {
            px = ((double)i * 2.0 / (double)SIZE - 1.0) * ZOOM_XZ;
            for (int32_t k = -SIZE; k < SIZE; ++k) {
                pz       = ((double)k * 2.0 / (double)SIZE - 1.0) * ZOOM_XZ;
                an       = std::atan2(px, pz);
                cx       = mix(0.2, -0.5, std::sin(an * 2.0));
                cy       = mix(0.5, 0.0, std::sin(an * 3.0));
                path     = sqrt(px * px + pz * pz) - 3.0;
                rot2D[0] = std::cos(an + time);
                rot2D[1] = -std::sin(an + time);
                rot2D[2] = std::sin(an + time);
                rot2D[3] = std::cos(an + time);
                for (int32_t j = -SIZE; j < SIZE; ++j) {
                    tmp_y = ((double)j * 2.0 / (double)SIZE - 1.0) * ZOOM_Y;
                    tmp_x = path;
                    rev_x = rot2D[0] * tmp_x + rot2D[1] * tmp_y;  // rx2d
                    rev_y = rot2D[2] * tmp_x + rot2D[3] * tmp_y;  // ry2d
                    kk    = 1.0;
                    hh    = 1.0;
                    for (int32_t idx = 0; idx < ITERATIONS; ++idx) {
                        rev_x_squared = rev_x * rev_x;
                        rev_y_squared = rev_y * rev_y;
                        hh *= 4.0 * kk;
                        kk = rev_x_squared + rev_y_squared;
                        if (kk > 4.0) {
                            break;
                        }
                        rev_y = 2.0 * rev_x * rev_y + cy;
                        rev_x = rev_x_squared - rev_y_squared + cx;
                    }
                    df = std::abs(sqrt(kk / hh) * std::log10(kk)) - 0.01;
                    if (df < 0.0) {
                        if (first_offset) {
                            first_offset = false;
                            offset[0]    = i;
                            offset[1]    = k;
                            offset[2]    = j;
                        }
                        cube_color = (int32_t)((std::sin(rev_x + rev_y) * 0.5 + 0.5) * 6.0) + 249;
                        densityLayerPtr->addVoxel(i + SIZE - offset[0], k + SIZE - offset[1], j + SIZE - offset[2], 1.0f);
                        sdfLayerPtr->addVoxel(i + SIZE - offset[0], k + SIZE - offset[1], j + SIZE - offset[2], df);
                    }
                }
            }
        }
    }
    vdb.saveToFile("julia_revolute.vdb");
}
#endif

#ifdef USE_ANIMATED_WAVE
int main() {
    const int32_t SIZE      = 150;
    const double  D_SIZE    = (double)SIZE;
    const int32_t OFFSET    = SIZE;
    const float   Z_SCALE   = 0.5f;
    const int32_t FRAMES    = 10;
    const float   len_ratio = 1.0f / (SIZE * SIZE);
    vdb::VDBWriter          vdb;
    float                   time        = 0.0f;
    std::array<float, 3>    colorF      = {120.0f, 240.8f, 20.2f};
    for (int32_t f = 0; f < FRAMES; ++f) {
        vdb.setKeyFrame(f);
        auto* floatLayerPtr = vdb.getFloatLayer(0, "density");
        auto* vec3sLayerPtr = vdb.getVec3sLayer(1, "color");
        for (int32_t i = -SIZE; i < SIZE; ++i) {
            for (int32_t j = -SIZE; j < SIZE; ++j) {
                float   len        = (i * i + j * j) * len_ratio;
                int32_t pz         = (int32_t)((std::sin(len * 10.0 + time) * 0.5 + 0.5) * (std::abs(50.0f - 25.0f * len)) * Z_SCALE);
                int32_t cube_color = (int32_t)(len * 100.0) % 255 + 1;
                auto    px         = i + SIZE;
                auto    py         = j + SIZE;
                floatLayerPtr->addVoxel(px, py, pz, 1.0f);
                float r  = sin(len * 10.0f) * 0.5f + 0.5f;
                float g = sin(len * 7.0f) * 0.5f + 0.5f;
                float b  = sin(len * 5.0f) * 0.5f + 0.5f;
                vec3sLayerPtr->addVoxel(px, py, pz, r, g, b);
            }
        }
        time += 0.5f;
    }
    vdb.saveToFile("wave.vdb");
}
#endif

#ifdef USE_VDB_WRITER
int main() {
    vdb::VDBWriter vdb;
    auto*          floatLayerPtr = vdb.getFloatLayer(0, "density");
    const uint32_t R = 128;
    const uint32_t D = R * 2;
    for (uint32_t z = 0; z < D; ++z) {
        for (uint32_t y = 0; y < D; ++y) {
            for (uint32_t x = 0; x < D; ++x) {
                const auto& px             = x - R;
                const auto& py             = y - R;
                const auto& pz             = z - R;
                const auto& length_squared = px * px + py * py + pz * pz;
                if (length_squared < R * R) {
                    floatLayerPtr->addVoxel(x, y, z, 1.0f);
                    float       fx  = (float)px;
                    float       fy  = (float)py;
                    float       fz  = (float)pz;
                    const auto& len = sqrtf(fx * fx + fy * fy + fz * fz);
                    fx /= len;
                    fy /= len;
                    fz /= len;
                }
            }
        }
    }
    vdb.saveToFile("sphere.vdb");
}
#endif
