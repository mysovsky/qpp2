#ifndef FMT_QPP_TYPES_HPP
#define FMT_QPP_TYPES_HPP

#include <fmt/format.h>
#include <mathf/lace3d.hpp>

namespace fmt {

//Function for formating custom class(here is generic_matrix) in fmt (for new version of fmt)
template <typename Scalar, int N, int M>
struct formatter<qpp::generic_matrix<Scalar, N, M>> : formatter<std::string> {
    template <typename FormatContext>
    auto format(const qpp::generic_matrix<Scalar, N, M>& mat, FormatContext& ctx) const {
        std::string str;
        if constexpr (M == 1) {
            str = mat.to_string_vec();
        } else {
            str = mat.to_string_matr();
        }
        return formatter<std::string>::format(str, ctx);
    }
};

} 
#endif
