#pragma once

#include <unordered_map>
#include <functional>
#include <iostream>
#include <cstdint>
#include <sstream>
#include <cstring>
#include <string>
#include <cstdio>
#include <chrono>
#include <array>
#include <map>

namespace vdb {

template <typename T>
inline T mini(const T& a, T& b) {
    return a < b ? a : b;
}
template <typename T>
inline T maxi(const T& a, T& b) {
    return a > b ? a : b;
}
template <typename T>
inline T clamp(const T& n) {
    return n >= T(0) && n <= T(1) ? n : T(n > T(0));
}  // clamp n => 0 to 1
template <typename T>
inline T clamp(const T& n, const T& b) {
    return n >= T(0) && n <= b ? n : T(n > T(0)) * b;
}  // clamp n => 0 to b
template <typename T>
inline T clamp(const T& n, const T& a, const T& b) {
    return n >= a && n <= b ? n : n < a ? a : b;
}  // clamp n => a to b

// specialized
struct dvec3 {
    double x, y, z;
    dvec3() { x = 0.0, y = 0.0, z = 0.0; }
    dvec3(const double& vxyz) { x = vxyz, y = vxyz, z = vxyz; }
    dvec3(const double& vx, const double& vy, const double& vz) { x = vx, y = vy, z = vz; }
    void operator+=(const double v) {
        x += v;
        y += v;
        z += v;
    }
    void operator-=(const double v) {
        x -= v;
        y -= v;
        z -= v;
    }
    void operator+=(const dvec3 v) {
        x += v.x;
        y += v.y;
        z += v.z;
    }
    void operator-=(const dvec3 v) {
        x -= v.x;
        y -= v.y;
        z -= v.z;
    }
    void operator*=(double v) {
        x *= v;
        y *= v;
        z *= v;
    }
    void operator/=(double v) {
        x /= v;
        y /= v;
        z /= v;
    }
    void operator*=(dvec3 v) {
        x *= v.x;
        y *= v.y;
        z *= v.z;
    }
    void operator/=(dvec3 v) {
        x /= v.x;
        y /= v.y;
        z /= v.z;
    }
};
inline dvec3 operator+(const dvec3& v, const double& f) { return dvec3(v.x + f, v.y + f, v.z + f); }
inline dvec3 operator+(const dvec3& v, dvec3 f) { return dvec3(v.x + f.x, v.y + f.y, v.z + f.z); }
inline dvec3 operator-(const dvec3& v, const double& f) { return dvec3(v.x - f, v.y - f, v.z - f); }
inline dvec3 operator-(const dvec3& v, dvec3 f) { return dvec3(v.x - f.x, v.y - f.y, v.z - f.z); }
inline dvec3 operator*(const dvec3& v, const double& f) { return dvec3(v.x * f, v.y * f, v.z * f); }
inline dvec3 operator*(const dvec3& v, dvec3 f) { return dvec3(v.x * f.x, v.y * f.y, v.z * f.z); }
inline dvec3 operator/(const dvec3& v, const double& f) { return dvec3(v.x / f, v.y / f, v.z / f); }
inline dvec3 operator/(dvec3& v, const double& f) { return dvec3(v.x / f, v.y / f, v.z / f); }
inline dvec3 operator/(const double& f, dvec3& v) { return dvec3(f / v.x, f / v.y, f / v.z); }
inline dvec3 operator/(const dvec3& v, dvec3 f) { return dvec3(v.x / f.x, v.y / f.y, v.z / f.z); }

// specialized
struct dAABBCC  // copy of b2AABB struct
{
    dvec3 lowerBound;  ///< the lower left vertex
    dvec3 upperBound;  ///< the upper right vertex

    dAABBCC() : lowerBound(0.0), upperBound(0.0) {}
    dAABBCC(dvec3 vlowerBound, dvec3 vUpperBound) {
        lowerBound = vlowerBound;
        upperBound = vUpperBound;
    }
    /// Add a vector to this vector.
    void operator+=(const dvec3& v) {
        lowerBound += v;
        upperBound += v;
    }

    /// Subtract a vector from this vector.
    void operator-=(const dvec3& v) {
        lowerBound -= v;
        upperBound -= v;
    }

    /// Multiply this vector by a scalar.
    void operator*=(double a) {
        lowerBound *= a;
        upperBound *= a;
    }

    /// Divide this vector by a scalar.
    void operator/=(double a) {
        lowerBound /= a;
        upperBound /= a;
    }

    /// Get the center of the AABB.
    const dvec3 GetCenter() const { return (lowerBound + upperBound) * 0.5; }

    /// Get the extents of the AABB (half-widths).
    const dvec3 GetExtents() const { return (upperBound - lowerBound) * 0.5; }

    /// Get the perimeter length
    double GetPerimeter() const {
        double wx = upperBound.x - lowerBound.x;
        double wy = upperBound.y - lowerBound.y;
        double wz = upperBound.z - lowerBound.z;
        return 2.0 * (wx + wy + wz);
    }

    /// Combine a point into this one.
    void Combine(dvec3 pt) {
        lowerBound.x = mini<double>(lowerBound.x, pt.x);
        lowerBound.y = mini<double>(lowerBound.y, pt.y);
        lowerBound.z = mini<double>(lowerBound.z, pt.z);
        upperBound.x = maxi<double>(upperBound.x, pt.x);
        upperBound.y = maxi<double>(upperBound.y, pt.y);
        upperBound.z = maxi<double>(upperBound.z, pt.z);
    }

    /// Does this aabb contain the provided vec2.
    bool ContainsPoint(const dvec3& pt) const {
        bool result = true;
        result      = result && lowerBound.x <= pt.x;
        result      = result && lowerBound.y <= pt.y;
        result      = result && lowerBound.z <= pt.z;
        result      = result && pt.x <= upperBound.x;
        result      = result && pt.y <= upperBound.y;
        result      = result && pt.z <= upperBound.z;
        return result;
    }

    bool Intersects(const dAABBCC& other) {
        bool result = true;
        result      = result || lowerBound.x <= other.lowerBound.x;
        result      = result || lowerBound.y <= other.lowerBound.y;
        result      = result || lowerBound.z <= other.lowerBound.z;
        result      = result || other.upperBound.x <= upperBound.x;
        result      = result || other.upperBound.y <= upperBound.y;
        result      = result || other.upperBound.z <= upperBound.z;
        return result;
    }

    const dvec3 Size() const { return dvec3(upperBound - lowerBound); }
};

/// Add a float to a dAABBCC.
inline dAABBCC operator+(const dAABBCC& v, float f) { return dAABBCC(v.lowerBound + f, v.upperBound + f); }

/// Add a dAABBCC to a dAABBCC.
inline dAABBCC operator+(const dAABBCC& v, dAABBCC f) { return dAABBCC(v.lowerBound + f.lowerBound, v.upperBound + f.upperBound); }

/// Substract a float from a dAABBCC.
inline dAABBCC operator-(const dAABBCC& v, float f) { return dAABBCC(v.lowerBound - f, v.upperBound - f); }

/// Substract a dAABBCC to a dAABBCC.
inline dAABBCC operator-(const dAABBCC& v, dAABBCC f) { return dAABBCC(v.lowerBound - f.lowerBound, v.upperBound - f.upperBound); }

/// Multiply a float with a dAABBCC.
inline dAABBCC operator*(const dAABBCC& v, float f) { return dAABBCC(v.lowerBound * f, v.upperBound * f); }

/// Multiply a dAABBCC with a dAABBCC.
inline dAABBCC operator*(const dAABBCC& v, dAABBCC f) { return dAABBCC(v.lowerBound * f.lowerBound, v.upperBound * f.upperBound); }

/// Divide a dAABBCC by a float.
inline dAABBCC operator/(const dAABBCC& v, float f) { return dAABBCC(v.lowerBound / f, v.upperBound / f); }

/// Divide a dAABBCC by a float.
inline dAABBCC operator/(dAABBCC& v, float f) { return dAABBCC(v.lowerBound / f, v.upperBound / f); }

/// Divide a dAABBCC by a dAABBCC.
inline dAABBCC operator/(const dAABBCC& v, dAABBCC f) { return dAABBCC(v.lowerBound / f.lowerBound, v.upperBound / f.upperBound); }

struct Node3 {
    uint64_t mask[8]   = {};
    float    data[512] = {};
};

struct Node4 {
    uint64_t                  mask[64] = {};
    std::map<uint32_t, Node3> nodes;
};

struct Node5 {
    uint64_t                  mask[512] = {};
    std::map<uint32_t, Node4> nodes;
};

struct VDB {
    Node5 nodes;
};

typedef std::array<float, 3> fVec3;
typedef std::array<float, 4> fVec4;

typedef std::array<std::array<double, 4>, 4> Mat4x4;

typedef uint32_t KeyFrame;
typedef dAABBCC  Volume;

typedef std::function<void(const KeyFrame& vKeyFrame, const double& vValue)> KeyFrameTimeLoggingFunctor;

class VdbWriter {
private:
    std::unordered_map<KeyFrame, VDB> m_Vdbs;
    FILE*   m_File    = nullptr;
    int32_t lastError = 0;
    Volume  maxVolume = Volume(1e7, -1e7);
    KeyFrame m_KeyFrame = 0;

    bool m_TimeLoggingEnabled = false;  // for log elapsed time between key frames and total

    std::chrono::steady_clock::time_point m_StartTime;
    std::chrono::steady_clock::time_point m_LastKeyFrameTime;
    std::map<KeyFrame, double>            m_FrameTimes;
    double                                m_TotalTime;

    KeyFrameTimeLoggingFunctor m_KeyFrameTimeLoggingFunctor;

public:
    void addVoxelDensity(const uint32_t& x, const uint32_t& y, const uint32_t& z, float value);
    void addVoxelNormal(const uint32_t& x, const uint32_t& y, const uint32_t& z, fVec3 normal);
    void addVoxelColor(const uint32_t& x, const uint32_t& y, const uint32_t& z, fVec4 color);
    void startTimeLogging();
    void stopTimeLogging();
    void setKeyFrameTimeLoggingFunctor(const KeyFrameTimeLoggingFunctor& vKeyFrameTimeLoggingFunctor);
    void setKeyFrame(uint32_t vKeyFrame);
    void saveToFile(const std::string& vFilePathName);
    void printStats() const;

private:
    void m_WriteNode5Header(FILE* fp, const Node5& node);
    void m_WriteNode4Header(FILE* fp, const Node4& node);
    void m_WriteTree(FILE* fp, VDB* vdb);
    void m_WriteMetadata(FILE* fp);
    void m_WriteTransform(FILE* fp, Mat4x4 mat);
    void m_WriteGrid(FILE* fp, VDB* vdb, Mat4x4 mat);
    void m_WriteVdb(FILE* fp, VDB* vdb, Mat4x4 mat);
    bool m_OpenFileForWriting(const std::string& vFilePathName);
    void m_CloseFile();
    long m_GetFilePos() const;
    void m_SetFilePos(const long& vPos);
};
}  // namespace vdb
