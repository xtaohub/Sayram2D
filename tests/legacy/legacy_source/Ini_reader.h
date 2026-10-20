/*
 * File:    ini_reader.h
 * Author:  Xin Tao <xtao@ustc.edu.cn>
 * Date: 2020-06-28 11:16 
 *
 * Copyright (c) Xin Tao
 *
 */

#ifndef INI_READER_H_
#define INI_READER_H_

#include <string>
#include <sstream>
#include "ini.h"

class Ini_reader{
 public:

  mINI::INIStructure ini; // this is made public, so we might call it directly. 
  
  Ini_reader(const std::string& filename):
    file_(filename) 
  {
    file_.read(ini); 
  }

  
  template<typename T>
  void read(const std::string& section, const std::string& key, T* valuep) {
    
    T& value = *valuep;

    bool found_section = ini.has(section);
    if (!found_section){
      throw section_not_found(section); 
    }

    bool found_key = ini[section].has(key); 
    if (!found_key){
      throw key_not_found(key);
    }    

    string_as_T<T>(ini.get(section).get(key), value);
    
  }

  // a short cut: set section first, then read key in this section
  void set_section(const std::string& section) { section_ = section; }

  template<typename T>
  void read(const std::string& key, T* valuep) {
    read(section_, key, valuep); 
  }

  struct section_not_found {
    std::string section;
    section_not_found(const std::string& section_ = string())
      :section(section_) {}
  };
  
  struct key_not_found {
    string key;
    key_not_found(const string & key_ = string())
      :key(key_)
    {}
  };
  

protected:
  template < class T > static void string_as_T(const string & s, T& t);
  
private:

  mINI::INIFile file_;
  std::string section_; 

};

// These are from github ConfigFile class 
/* static */
template < class T > void Ini_reader::string_as_T (const string & s, T& t) {
  // Convert from a string to a T
  // Type T must support >> operator
  std::istringstream ist(s);
  ist >> t;
}

template <> inline void Ini_reader::string_as_T < string > (const string & sin, string& sout) {
  // Convert from a string to a string
  // In other words, do nothing
  sout = sin.c_str();
}

/* static */
template <> inline void Ini_reader::string_as_T < bool > (const string & s, bool& b) {
  using std::cout;
  using std::endl;

  // Convert from a string to a bool
  // Interpret "false", "F", "no", "n", "0" as false
  // Interpret "true", "T", "yes", "y", "1", "-1", or anything else as true
  b = true;
  string sup = s;
  for (string::iterator p = sup.begin(); p != sup.end(); ++p)
    *p = toupper(*p);           // make string all caps

  if (sup == string("FALSE") || sup == string("F")
      || sup == string("NO") || sup == string("N")
      || sup == string("0") || sup == string("NONE")){

    b = false;
  }

}


#endif /* INI_READER_H_ */
