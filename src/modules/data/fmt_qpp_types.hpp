#ifndef FMT_QPP_TYPES_HPP
#define FMT_QPP_TYPES_HPP

#include <symm/index.hpp>
#include <data/data.hpp>
#include <mathf/lace3d.hpp>

#include <fmt/format.h>
#include <fmt/ostream.h>

namespace fmt {

// Generic Matrix
template <typename Scalar, int N, int M>
struct formatter<qpp::generic_matrix<Scalar, N, M>> : ostream_formatter {};

// Index
template <>
struct formatter<qpp::index> : ostream_formatter {};

// Bool
template <>
struct formatter<qpp::Bool> : formatter<bool> {
    auto format(const qpp::Bool& b, format_context& ctx) const {
        return formatter<bool>::format(static_cast<bool>(b), ctx);
    }
};

// Vector3 
template <typename T>
struct formatter<qpp::vector3<T>> : ostream_formatter {};

// Basic Types
template <>
struct formatter<qpp::basic_types> : formatter<int> {
    auto format(qpp::basic_types bt, format_context& ctx) const {
        return formatter<int>::format(static_cast<int>(bt), ctx);
    }
};

} // namespace fmt

#endif
