#ifndef QPP_STRFUN
#define QPP_STRFUN

#include <data/types.hpp>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <mathf/lace3d.hpp>
#include <io/simplefun.hpp>
#include <string_view>
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)

#include <pybind11/pybind11.h>

namespace py = pybind11;
#endif
namespace qpp {

  template<class REAL>
  qpp::vector3<REAL> vec_from_string(STRING_EX &_inst, int idx = 0, int idy = 1, int idz = 2 ) {

    std::vector<STRING_EX> vfs_l = split(_inst);
    REAL vx, vy, vz = 0.0;
    s2t(vfs_l[idx], vx);
    s2t(vfs_l[idy], vy);
    s2t(vfs_l[idz], vz);
    return qpp::vector3<REAL>(vx, vy, vz);

  }

  void replace_string_inplace(STRING_EX& subject, const STRING_EX& search,
                              const STRING_EX& replace);

  //https://www.bfilipek.com/2018/07/string-view-perf-followup.html
  std::vector<std::string_view> split_sv(std::string_view strv, std::string_view delims = " ");

  char *vec_str_to_char(const STRING_EX & s);
  const char *vec_str_to_char_ref(const STRING_EX & s);

  template<typename charT>
  struct ci_equal {
      ci_equal( const std::locale& loc ) : loc_(loc) {}
      bool operator() (charT ch1, charT ch2) {
        return std::tolower(ch1, loc_) == std::tolower(ch2, loc_);
      }
    private:
      const std::locale& loc_;
  };

  template<typename T>
  bool find_string_ci( const T& str1, const T& str2, const std::locale& loc = std::locale()) {
    typename T::const_iterator it = std::search( str1.begin(), str1.end(),
                                                 str2.begin(), str2.end(),
                                                 ci_equal<typename T::value_type>(loc) );
    if ( it != str1.end() ) return true/*it - str1.begin()*/;
    else return false/*-1*/; // not found
    return false;
  }

  int common_begin(const STRING_EX & s1, const STRING_EX & s2);

  // fixme - place it in different module!
  template<class T>
  //bool oneof(const T & t, std::initializer_list<T> list ){
  bool oneof(const T & t, const std::vector<T> & list ){
    for (const T & elem : list)
      if (t==elem)
	return true;
    return false;
  }

  STRING_EX atomic_name_to_symbol(const STRING_EX & nm);
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)

  int fileno(py::object &f);
#endif  
}

#endif
