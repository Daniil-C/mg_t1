#ifndef __GEOM_H__
#define __GEOM_H__

#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>

template <size_t N, typename T> 
class vec {
    T data[N];
public:
    vec() { 
        for (size_t i = 0; i < N; i++)
            data[i] = T();
    }

    vec(T n1, T n2, T n3) {
        data[0] = n1;
        data[1] = n2;
        data[2] = n3;
    }

    vec(T n1, T n2, T n3, T n4) {
        data[0] = n1;
        data[1] = n2;
        data[2] = n3;
        data[3] = n4;
    }
    
    T& operator[](const size_t i) {
        assert(i < N);
        return data[i];
    }
    const T& operator[](const size_t i) const {
        assert(i < N);
        return data[i];
    }

    float norm() {
        float sum = 0;
        for (size_t i = 0; i < N; i++) {
            sum += data[i] * data[i];
        }
        return std::sqrt(sum);
    }
    vec<N,T> & normalize(T l = 1) {
        *this = (*this) * (l / norm());
        return *this;
    }
};

template<size_t N, typename T>
T operator*(const vec<N, T>& l, const vec<N, T>& r) {
    T res = T();
    for (size_t i = 0; i < N; i++)
        res += l[i] * r[i];
    return res;
}

template<size_t N, typename T>
vec<N, T> operator+(vec<N, T> l, const vec<N, T>& r) {
    for (size_t i = 0; i < N; i++)
        l[i] += r[i];
    return l;
}

template<size_t N, typename T>
vec<N, T> operator-(vec<N, T> l, const vec<N, T>& r) {
    for (size_t i = 0; i < N; i++)
        l[i] -= r[i];
    return l;
}

template<size_t N, typename T, typename U> 
vec<N, T> operator*(const vec<N, T> &l, const U& r) {
    vec<N, T> res;
    for (size_t i = 0; i < N; i++)
        res[i] = l[i] * r;
    return res;
}

template<size_t N, typename T>
vec<N,T> operator-(const vec<N, T> &l) {
    return l * T(-1);
}

template <typename T> 
vec<3, T> cross(vec<3, T> v1, vec<3, T> v2) {
    return vec<3, T>(v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0]);
}

#endif