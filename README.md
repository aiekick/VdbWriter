# VdbWriter

VdbWriter is a specific and limited writer for openvdb.

He is based on https://github.com/jangafx/simple-vdb-writer
 
## App

the main.cpp file show you how to generate a quick file :

```cpp
#include "VdbWriter.h"
int main() 
{
    const int32_t  SIZE      = 189;
    const int32_t  OFFSET    = SIZE;
    const float    Z_SCALE   = 1.0f;
    const int32_t  FRAMES    = 30;
    const float    len_ratio = 1.0f / (SIZE * SIZE);
	
    vdb::VdbWriter vdb;
    
    vdb.StartTimeLogging();
    float time = 0.0f;
    for (int32_t k = 0; k < FRAMES; ++k) {
        vdb.SetKeyFrame(k);
        for (int32_t i = -SIZE; i < SIZE; ++i) {
            for (int32_t j = -SIZE; j < SIZE; ++j) {
                float   len        = (i * i + j * j) * len_ratio;
                int32_t pz         = (int32_t)((std::sin(len * 10.0 + time) * 0.5 + 0.5) * (std::abs(50.0f - 25.0f * len)) * Z_SCALE);
                vdb.addVoxelDensity(i + OFFSET, j + OFFSET, pz, 1.0f);  // blender3D use the z as up axis
            }
        }
        time += 0.5f;
    }
    vdb.StopTimeLogging();
    vdb.SaveToFile("output.vdb");
    vdb.PrintStats();
}
```

### Output in blender3D :

![alt](https://github.com/aiekick/VdbWriter/blob/master/doc/wave_blender.jpg)

### Another Output in blender3D with my julia revolute

![alt](https://github.com/aiekick/VdbWriter/blob/master/doc/julia_revolute_blender.jpg)
