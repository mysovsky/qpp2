#ifndef QPP_SIMPLEFUN
#define QPP_SIMPLEFUN

#include <data/types.hpp>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iostream>
//#include <mathf/lace3d.hpp>
#include <string_view>

namespace qpp {

  // -------------------- Simple tokenizer -----------------------------------

  class tokenizer {
      std::basic_istream<CHAR_EX,TRAITS> * _input;
      STRING_EX _buff, _dump, _sepr;
      int _line_number;
      STRING_EX _filename;
      bool _created_here;

    public:

      tokenizer(ISTREAM & input, const STRING_EX & __filename="") {
        _input = & input;
        _dump = " \t";
        _line_number = 0;
        _filename = __filename;
        _created_here = false;
      }
      /*
        tokenizer(const STRING_EX & str)
        {
          _input = new std::basic_stringstream<CHAR,TRAITS>(str);
          _dump = " \t";
          _line_number = 0;
          _filename = "";
          _created_here = true;
        }
        */
      tokenizer(const STRING_EX & __filename){
        _input = new
                 std::basic_ifstream<CHAR_EX,TRAITS>(__filename.c_str());
        _dump = " \t";
        _line_number = 0;
        _filename = __filename;
        _created_here = true;
      }

      ~tokenizer(){
        if (_created_here)
          delete _input;
      }

      void dump(const STRING_EX & smb){
        _dump = smb;
      }

      void separate(const STRING_EX & smb){
        _sepr = smb;
      }

      STRING_EX get(){
        int i;
        if (_buff == "" ) {
            std::getline(*_input, _buff);

            _line_number++;
          }

        do {
            i = _buff.find_first_not_of(_dump);

            //debug
            //std::cout << "i = " << i << "\""  << _buff << "\"\n";
            if (i != std::string::npos) {
                _buff = _buff.substr(i);
                break;
              }
            else if ( !_input -> eof() ){
                std::getline(*_input, _buff);
                _line_number++;
              }
            else{
                _buff = "";
                break;
              }
          } while ( true );

        if ( _input -> eof() && _buff.size()==0 )
          return "";

        //debug
        //std::cout << "\"" << _buff << "\"\n";

        STRING_EX rez;
        i = _buff.find_first_of(_sepr + _dump);
        if (i==0){
            //debug
            //std::cout << "here1\n";

            rez = _buff.substr(0,1);
            _buff = _buff.substr(1);

          }
        else if (i != std::string::npos) {
            //debug
            //std::cout << "here2\n";

            rez =  _buff.substr(0,i);
            _buff = _buff.substr(i);
          }
        else {
            //debug
            //std::cout << "here3\n";

            rez =  _buff;
            _buff = "";
          }
        return rez;
      }

      void back(STRING_EX s){
        _buff = s + " " + _buff;
      }

      bool eof() const{
        return _input -> eof() && _buff == "";
      }

      int line() const{
        return _line_number;
      }

      STRING_EX file() const{
        return _filename;
      }

  };

  // -----------------------------------------------------------

  STRING_EX tolower(const STRING_EX & s);
  // Make lowercase

  // -----------------------------------------------------------

  bool icompare(const STRING_EX & s1, const STRING_EX s2);
  // Case insensitive comparison of two strings
  // -----------------------------------------------------------

  void split(const STRING_EX &s, std::vector<STRING_EX> &elems, const STRING_EX & delims = " \t");
  // fixme - not efficient!

  std::vector<STRING_EX> split(const STRING_EX &s, const STRING_EX & delims=" \t");

  bool is_identifier(const STRING_EX &s);

  // --------------------------------------------------------------------//

  int strnf(const STRING_EX & s);

  // -------------------------------- string to type T convertor ----------------------------

  template<typename T>
  bool s2t(const STRING_EX & s, T & val) {
    std::basic_stringstream<CHAR_EX,TRAITS> ss(s);
    ss >> val;

    //std::cout << "ss eof= " << ss.eof() << "\n";
    return (!ss.fail()) && ss.eof();
  }

  template<>
  bool s2t<bool>(const STRING_EX & s, bool & val);

  // -------------------------------------------------------------

  template<typename T>
  STRING_EX t2s(const T & val) {

    std::basic_stringstream<CHAR_EX,TRAITS> ss;
    ss << val;
    return ss.str();

  }

  template<>
  STRING_EX t2s<bool>(const bool & val);
  template<>
  STRING_EX t2s<std::complex<float>>(const std::complex<float> & val);
  template<>
  STRING_EX t2s<std::complex<double>>(const std::complex<double> & val);
  

  // -------------------------------------------------------------


  STRING_EX extract_base_name(STRING_EX const & path);

};
#endif
