#include "VdbWriter.h"
#include <vector>
#include <iomanip>
#include <iostream>

namespace vdb {

void write_ptr(FILE* fp, const void* data, size_t elementSize, size_t count = 1) { fwrite(data, elementSize, count, fp); }
void write_ptr(FILE* fp, void* data, size_t elementSize, size_t count = 1) { fwrite(data, elementSize, count, fp); }
template<typename T>
void write_data(FILE* fp, T data) {
    write_ptr(fp, &data, sizeof(T));
}
template <typename T>
void write_data_arr(FILE* fp, T* data, size_t count) {
    write_ptr(fp, data, sizeof(T), count);
}
template <typename T>
void write_data_arr(FILE * fp, const T* data, size_t count) {
    write_ptr(fp, data, sizeof(T), count);
}

void write_string(FILE* fp, std::string str) { write_ptr(fp, str.data(), sizeof(uint8_t), str.size()); }
void write_vec3i(FILE* fp, const std::array<int32_t, 3>& data) { write_data_arr<int32_t>(fp, data.data(), 3); }

void write_name(FILE* fp, const std::string& name) {
    write_data<uint32_t>(fp, (uint32_t)name.size());
    write_string(fp, name);
}

void write_meta_string(FILE* fp, const std::string& name, const std::string& str) {
    write_name(fp, name);
    write_name(fp, "string");
    write_name(fp, str);
}

void write_meta_bool(FILE* fp, const std::string& name, bool flag) {
    write_name(fp, name);
    write_name(fp, "bool");
    write_data<uint32_t>(fp, 1);  // one byte is used to store the bool
    write_data<uint8_t>(fp, flag ? 1 : 0);
}

void write_meta_vec3i(FILE* fp, const std::string& name, const std::array<int32_t, 3>& data) {
    write_name(fp, name);
    write_name(fp, "vec3i");
    write_data<uint32_t>(fp, 12U);  // 12 bytes (4 bytes * 3) are used to store the vec3i
    write_vec3i(fp, data);
}

uint32_t getBitIndex4(uint32_t x, uint32_t y, uint32_t z) {
    x &= (uint32_t)(4096 - 1);
    y &= (uint32_t)(4096 - 1);
    z &= (uint32_t)(4096 - 1);
    uint32_t idx_3d[3] = {x >> 7, y >> 7, z >> 7};
    uint64_t idx       = idx_3d[2] | (idx_3d[1] << 5) | (idx_3d[0] << 10);
    return (uint32_t)idx;
}

uint32_t getBitIndex3(uint32_t x, uint32_t y, uint32_t z) {
    x &= (uint32_t)(128 - 1);
    y &= (uint32_t)(128 - 1);
    z &= (uint32_t)(128 - 1);
    uint32_t idx_3d[3] = {x >> 3, y >> 3, z >> 3};
    uint64_t idx       = idx_3d[2] | (idx_3d[1] << 4) | (idx_3d[0] << 8);
    return (uint32_t)idx;
}

uint32_t getBitIndex0(uint32_t x, uint32_t y, uint32_t z) {
    x &= (uint32_t)(8 - 1);
    y &= (uint32_t)(8 - 1);
    z &= (uint32_t)(8 - 1);
    uint32_t idx_3d[3] = {x >> 0, y >> 0, z >> 0};
    uint64_t idx       = idx_3d[2] | (idx_3d[1] << 3) | (idx_3d[0] << 6);
    return (uint32_t)idx;
}

void VdbWriter::saveToFile(const std::string& vFilePathName) {
    if (!vFilePathName.empty()) {
        auto dot_p = vFilePathName.find_last_of('.');
        if (dot_p != std::string::npos) {
            auto   base_file_path_name = vFilePathName.substr(0, dot_p);
            size_t idx                 = 1;
            for (auto& vdb : m_Vdbs) {
                std::stringstream str;
                if (m_Vdbs.size() > 1) { // many frames
                    str << base_file_path_name << "_" << std::setfill('0') << std::setw(4) << idx++ << ".vdb";
                }
                else {
                    str << base_file_path_name << ".vdb";
                }
                if (m_OpenFileForWriting(str.str())) {
                    m_WriteVdb(m_File, &vdb.second);
                    m_CloseFile();
                }
            }
        }
    }
}

void VdbWriter::printStats() const {
    std::cout << "---- Stats ------------------------------" << std::endl;
    std::cout << "Volume : " << maxVolume.Size().x << " x " << maxVolume.Size().y << " x " << maxVolume.Size().z << std::endl;
    /*std::cout << "count cubes : " << cubes.size() << std::endl;
    std::map<KeyFrame, size_t> frame_counts;
    for (const auto& cube : cubes) {
        for (auto& key_xyzi : cube.xyzis) {
            frame_counts[key_xyzi.first] += key_xyzi.second.numVoxels;
        }
    }*/
    /*size_t voxels_total = 0U;
    if (frame_counts.size() > 1U) {
        std::cout << "count key frames : " << frame_counts.size() << std::endl;
        std::cout << "-----------------------------------------" << std::endl;
        for (const auto& frame_count : frame_counts) {
            std::cout << " o--\\-> key frame : " << frame_count.first << std::endl;
            std::cout << "     \\-> voxels count : " << frame_count.second << std::endl;
            if (m_FrameTimes.find(frame_count.first) != m_FrameTimes.end()) {
                std::cout << "      \\-> elapsed time : " << m_FrameTimes.at(frame_count.first) << " secs" << std::endl;
            }
            voxels_total += frame_count.second;
        }
        std::cout << "-----------------------------------------" << std::endl;
    } else if (!frame_counts.empty()) {
        voxels_total = frame_counts.begin()->second;
    }
    std::cout << "voxels total : " << voxels_total << std::endl;*/
    std::cout << "total elapsed time : " << m_TotalTime << " secs" << std::endl;
    std::cout << "-----------------------------------------" << std::endl;
}

void VdbWriter::addLayer(const size_t& layer, const std::string layerName) { m_Labels[layer] = layerName; }

void VdbWriter::addVoxelFloat(const uint32_t& x, const uint32_t& y, const uint32_t& z, float value, size_t layer) {
    maxVolume.Combine(dvec3((float)x, (float)y, (float)z));
    auto  bit_index_4 = getBitIndex4(x, y, z);
    auto  bit_index_3 = getBitIndex3(x, y, z);
    auto  bit_index_0 = getBitIndex0(x, y, z);
    auto& nodes5Ref   = m_Vdbs[m_KeyFrame].nodes[layer];
    auto& nodes4Ref   = nodes5Ref.nodes[bit_index_4];
    auto& nodes3Ref   = nodes4Ref.nodes[bit_index_3];
    nodes5Ref.mask[bit_index_4 >> 6] |= static_cast<uint64_t>(1) << (bit_index_4 & (64 - 1));
    nodes4Ref.mask[bit_index_3 >> 6] |= static_cast<uint64_t>(1) << (bit_index_3 & (64 - 1));
    nodes3Ref.mask[bit_index_0 >> 6] |= static_cast<uint64_t>(1) << (bit_index_0 & (64 - 1));
    nodes3Ref.data[bit_index_0] = value;
}

void VdbWriter::startTimeLogging() {
    m_TimeLoggingEnabled = true;
    m_StartTime          = std::chrono::steady_clock::now();
    m_LastKeyFrameTime   = m_StartTime;
};

void VdbWriter::stopTimeLogging() {
    if (m_TimeLoggingEnabled) {
        const auto now           = std::chrono::steady_clock::now();
        m_FrameTimes[m_KeyFrame] = std::chrono::duration_cast<std::chrono::milliseconds>(now - m_LastKeyFrameTime).count() * 1e-3;
        if (m_KeyFrameTimeLoggingFunctor) {
            m_KeyFrameTimeLoggingFunctor(m_KeyFrame, m_FrameTimes.at(m_KeyFrame));
        }
        m_TotalTime          = std::chrono::duration_cast<std::chrono::milliseconds>(now - m_StartTime).count() * 1e-3;
        m_TimeLoggingEnabled = false;
    }
}

void VdbWriter::setKeyFrameTimeLoggingFunctor(const KeyFrameTimeLoggingFunctor& vKeyFrameTimeLoggingFunctor) {
    m_KeyFrameTimeLoggingFunctor = vKeyFrameTimeLoggingFunctor;
}

void VdbWriter::setKeyFrame(uint32_t vKeyFrame) {
    if (m_KeyFrame != vKeyFrame) {
        if (m_TimeLoggingEnabled) {
            const auto now           = std::chrono::steady_clock::now();
            const auto elapsed       = now - m_LastKeyFrameTime;
            m_FrameTimes[m_KeyFrame] = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count() * 1e-3;
            if (m_KeyFrameTimeLoggingFunctor) {
                m_KeyFrameTimeLoggingFunctor(m_KeyFrame, m_FrameTimes.at(m_KeyFrame));
            }
            m_LastKeyFrameTime = now;
        }
        m_KeyFrame = vKeyFrame;
    }
}

// Routines for writing the actual format
void VdbWriter::m_WriteNode5Header(FILE* fp, const Node5& node) {
    // Origin of the 5-node
    write_vec3i(fp, {0, 0, 0});
    // Child masks
    write_data_arr<uint64_t>(fp, node.mask, 512);
    // Value masks are zero for now
    static std::vector<uint64_t> mask_arr_uint64_empty(512);  // we dont use a std::array for not increase the bin size
    write_data_arr<uint64_t>(fp, mask_arr_uint64_empty.data(), mask_arr_uint64_empty.size());
    // Write uncompressed node values, 6 means no compression
    write_data<uint8_t>(fp, 6);
    static std::vector<float> mask_arr_float_empty(32768);  // we dont use a std::array for not increase the bin size
    write_data_arr<float>(fp, mask_arr_float_empty.data(), mask_arr_float_empty.size());
}

void VdbWriter::m_WriteNode4Header(FILE* fp, const Node4& node) {
    // Child masks
    write_data_arr<uint64_t>(fp, node.mask, 64);
    // Value masks are zero for now
    static std::vector<uint64_t> mask_arr_uint64_empty(64);  // we dont use a std::array for not increase the bin size
    write_data_arr<uint64_t>(fp, mask_arr_uint64_empty.data(), mask_arr_uint64_empty.size());
    // Write uncompressed node values, 6 means no compression
    write_data<uint8_t>(fp, 6);
    static std::vector<float> mask_arr_float_empty(4096);  // we dont use a std::array for not increase the bin size
    write_data_arr<float>(fp, mask_arr_float_empty.data(), mask_arr_float_empty.size());
}

template <typename T>
int count_trailing_zeros(T value) {
    int count = 0;
    while ((value & 1) == 0 && value != 0) {
        ++count;
        value >>= 1;
    }
    return count;
}

void VdbWriter::m_WriteTree(FILE* fp, VDB* vdb, size_t layer) {
    // We need to write a 1, apparently
    write_data<uint32_t>(fp, 1);
    // Root node background value
    write_data<float>(fp, 0);
    // Number of tiles
    write_data<uint32_t>(fp, 0);
    // Number of 5-nodes
    write_data<uint32_t>(fp, 1);
    const auto& nodes5Ref = vdb->nodes[layer];
    m_WriteNode5Header(fp, nodes5Ref);
    // Iterate 4-nodes
    size_t word5_idx = 0;
    for (auto word5 : nodes5Ref.mask) {
        const auto& base_bit_4_idx = ((uint32_t)word5_idx++) * 64;
        for (; word5 != 0; word5 &= word5 - 1) {
            const auto& bit_4_index = base_bit_4_idx + (uint32_t)count_trailing_zeros(word5);
            const auto& nodes4Ref   = nodes5Ref.nodes.at(bit_4_index);
            m_WriteNode4Header(fp, nodes4Ref);
            // Iterate 3-nodes
            size_t word4_idx = 0;
            for (auto word4 : nodes4Ref.mask) {
                const auto& base_bit_3_idx = ((uint32_t)word4_idx++) * 64;
                for (; word4 != 0; word4 &= word4 - 1) {
                    const auto& bit_3_index = base_bit_3_idx + (uint32_t)count_trailing_zeros(word4);
                    const auto& nodes3Ref   = nodes4Ref.nodes.at(bit_3_index);
                    write_data_arr<uint64_t>(fp, nodes3Ref.mask, 8);
                }
            }
        }
    }
    // Iterate 4-nodes
    word5_idx = 0;
    for (auto word5 : nodes5Ref.mask) {
        const auto& base_bit_4_idx = ((uint32_t)word5_idx++) * 64;
        for (; word5 != 0; word5 &= word5 - 1) {
            const auto& bit_4_index = base_bit_4_idx + (uint32_t)count_trailing_zeros(word5);
            const auto& nodes4Ref   = nodes5Ref.nodes.at(bit_4_index);
            // Iterate 3-nodes
            size_t word4_idx = 0;
            for (auto word4 : nodes4Ref.mask) {
                const auto& base_bit_3_idx = ((uint32_t)word4_idx++) * 64;
                for (; word4 != 0; word4 &= word4 - 1) {
                    const auto& bit_3_index = base_bit_3_idx + (uint32_t)count_trailing_zeros(word4);
                    const auto& nodes3Ref   = nodes4Ref.nodes.at(bit_3_index);
                    write_data_arr<uint64_t>(fp, nodes3Ref.mask, 8);
                    // Write uncompressed node values, 6 means no compression
                    write_data<uint8_t>(fp, 6);
                    write_data_arr<float>(fp, nodes3Ref.data, 512);
                }
            }
        }
    }
}

void VdbWriter::m_WriteMetadata(FILE* fp, const std::string& layerName) {
    // Number of entries
    write_data<uint32_t>(fp, 5);
    write_meta_string(fp, "class", "unknown");
    write_meta_string(fp, "file_compression", "none");
    write_meta_vec3i(fp, "file_bbox_max", {(int32_t)maxVolume.upperBound.x, (int32_t)maxVolume.upperBound.y, (int32_t)maxVolume.upperBound.z});
    write_meta_vec3i(fp, "file_bbox_min", {(int32_t)maxVolume.lowerBound.x, (int32_t)maxVolume.lowerBound.y, (int32_t)maxVolume.lowerBound.z});
    write_meta_string(fp, "name", layerName);
}

void VdbWriter::m_WriteTransform(FILE* fp) {
    write_name(fp, "UniformScaleTranslateMap");
    // write Translation
    const auto& center = maxVolume.GetCenter();
    write_data<double>(fp, -center.x);
    write_data<double>(fp, -center.y);
    write_data<double>(fp, -center.z);
    // write ScaleValues
    write_data<double>(fp, 1.0);
    write_data<double>(fp, 1.0);
    write_data<double>(fp, 1.0);
    // write VoxelSize
    write_data<double>(fp, 1.0);
    write_data<double>(fp, 1.0);
    write_data<double>(fp, 1.0);
    // write ScaleValuesInverse
    write_data<double>(fp, 1.0);
    write_data<double>(fp, 1.0);
    write_data<double>(fp, 1.0);
    // write InvScaleSqr
    write_data<double>(fp, 1.0);
    write_data<double>(fp, 1.0);
    write_data<double>(fp, 1.0);
    // write InvTwiceScale;
    write_data<double>(fp, 0.5);
    write_data<double>(fp, 0.5);
    write_data<double>(fp, 0.5);

    /*write_name(fp, "AffineMap");
    //write_data_arr<double>(fp, &mat[0][0], 16);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            write_data<double>(fp, mat[i][j]);
        }
    }*/
}

void VdbWriter::m_WriteGrid(FILE* fp, VDB* vdb, size_t layer, const std::string& layerName) {
    // Grid name
    write_name(fp, layerName);
    // Grid type
    write_name(fp, "Tree_float_5_4_3");
    // No instance parent
    write_data<uint32_t>(fp, 0);
    // Grid descriptor stream position
    uint64_t stream_pos = m_GetFilePos();
    write_data<uint64_t>(fp, stream_pos + sizeof(uint64_t) * 3);
    write_data<uint64_t>(fp, 0);
    write_data<uint64_t>(fp, 0);
    // No compression
    write_data<uint32_t>(fp, 0);
    m_WriteMetadata(fp, layerName);
    m_WriteTransform(fp);
    m_WriteTree(fp, vdb, layer);
}

void VdbWriter::m_WriteVdb(FILE* fp, VDB* vdb) {
    // Magic number
    std::array<uint8_t, 8> header = {0x20, 0x42, 0x44, 0x56, 0x0, 0x0, 0x0, 0x0};
    write_ptr(fp, header.data(), sizeof(uint8_t), header.size());
    // File version
    write_data<uint32_t>(fp, 224);
    // Library version (we're just gonna pretend we're OpenVDB 8.1)
    write_data<uint32_t>(fp, 8);  // major
    write_data<uint32_t>(fp, 1);  // minor
    // We do not have grid offsets
    write_data<uint8_t>(fp, 0);
    // Temporary UUID
    write_string(fp, "00000000-0000-0000-0000-000000000000");
    // No metadata for now
    write_data<uint32_t>(fp, 0);
    // One grid
    write_data<uint32_t>(fp, (uint32_t)m_Labels.size());
    for (const auto& label : m_Labels) {
        m_WriteGrid(fp, vdb, label.first, label.second);
    }
}

bool VdbWriter::m_OpenFileForWriting(const std::string& vFilePathName) {
#if _MSC_VER
    lastError = fopen_s(&m_File, vFilePathName.c_str(), "wb");
#else
    m_File    = fopen(vFilePathName.c_str(), "wb");
    lastError = m_File ? 0 : errno;
#endif
    if (lastError != 0)
        return false;
    return true;
}

void VdbWriter::m_CloseFile() { fclose(m_File); }

long VdbWriter::m_GetFilePos() const { return ftell(m_File); }

void VdbWriter::m_SetFilePos(const long& vPos) {
    //  SEEK_SET	Beginning of file
    //  SEEK_CUR	Current position of the file pointer
    //	SEEK_END	End of file
    fseek(m_File, vPos, SEEK_SET);
}

}  // namespace vdb