# VdbWriter

VdbWriter is a specific and limited writer for openvdb.

He is based on https://github.com/jangafx/simple-vdb-writer

## Features :

* cpp11 minimal support
* one header file only
* dependency free (no openvdb dependencies)
* can write multi layers format :
  * float
  * double
  * vec3s
  * vec3d
  * vec3i
  * vec3ui
  * other format can be used by his template, but not sure if softwares can open it
* can set the current frame for animation recording +> will save in many files
 
## App

the main.cpp file show you how to generate a quick file :

```cpp
#include "VdbWriter.h"
int main() {
    const int32_t  SIZE      = 150;
    const double   D_SIZE    = (double)SIZE;
    const int32_t  OFFSET    = SIZE;
    const float    Z_SCALE   = 0.5f;
    const int32_t  FRAMES    = 10;
    const float    len_ratio = 1.0f / (SIZE * SIZE);
    vdb::VdbWriter vdb;
    float          r, g, b;
    float          time = 0.0f;
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
                r = sin(len * 10.0f) * 0.5f + 0.5f;
                g = sin(len * 7.0f) * 0.5f + 0.5f;
                b = sin(len * 5.0f) * 0.5f + 0.5f;
                vec3sLayerPtr->addVoxel(px, py, pz, r, g, b);
            }
        }
        time += 0.5f;
    }
    vdb.saveToFile("wave.vdb");
}
```

### Output in blender3D :

![alt](https://github.com/aiekick/VdbWriter/blob/master/doc/wave_blender.jpg)

### Another Output in blender3D with my julia revolute

![alt](https://github.com/aiekick/VdbWriter/blob/master/doc/julia_revolute_blender.jpg)
