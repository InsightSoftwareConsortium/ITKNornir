// -*- Mode: c++; tab-width: 8; c-basic-offset: 2; indent-tabs-mode: t -*-
// NOTE: the first line of this file sets up source code indentation rules
// for Emacs; it is also a hint to anyone modifying this file.

/*
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


// File         : the_text.hxx
// Author       : Pavel A. Koshevoy
// Created      : Sun Aug 29 14:53:00 MDT 2004
// Copyright    : (C) 2004-2008 University of Utah
// License      : GPLv2
// Description  : text convenience class.

#ifndef itkIRText_h
#define itkIRText_h

// system includes:
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>
#include <list>


//----------------------------------------------------------------
// the_text_t
// 
class the_text_t
{
public:
  the_text_t(const char * text = "");
  the_text_t(const char * text, const size_t & size);
  the_text_t(const the_text_t & text);
  the_text_t(const std::list<char> & text);
  ~the_text_t();
  
  // assignment operator:
  inline the_text_t & operator = (const the_text_t & text)
  {
    if (this != &text)
    {
      assign(text.text_, text.size_);
    }
    
    return *this;
  }
  
  // clear the string:
  inline void clear()
  {
    delete [] text_;
    text_ = nullptr;
    size_ = 0;
  }
  
  // shorthand:
  inline bool is_empty() const
  { return size_ == 0; }
  
  // assign a new string to this text:
  inline void assign(const char * text)
  { assign(text, strlen(text)); }
  
  void assign(const char * text, const size_t & text_size);
  
  // append a new string to this text:
  inline void append(const char * text)
  { append(text, strlen(text)); }
  
  void append(const char * text, const size_t & text_size);

  // replace
  void replace(const char find, const char replace);
  
  // equality/inequality tests:
  inline bool operator == (const the_text_t & text) const
  { return ((size_ == text.size_) && (strcmp(text_, text.text_) == 0)); }
  
  inline bool operator != (const the_text_t & text) const
  { return !(*this == text); }
  
  inline bool operator == (const char * text) const
  { return (*this == the_text_t(text)); }
  
  inline bool operator != (const char * text) const
  { return !(*this == text); }
  
  inline bool operator < (const the_text_t & text) const
  { return (strcmp(text_, text.text_) < 0); }
  
  inline bool operator > (const the_text_t & text) const
  { return (strcmp(text_, text.text_) > 0); }
  
  // arithmetic:
  inline the_text_t & operator += (const the_text_t & text)
  {
    append(text.text_, text.size_);
    return *this;
  }
  
  inline the_text_t operator + (const the_text_t & text) const
  {
    the_text_t text_sum(*this);
    text_sum += text;
    return text_sum;
  }
  
  // access operators:
  template <typename index_t>
  inline const char & operator [] (const index_t & index) const
  { return text_[index]; }
  
  template <typename index_t>
  inline char & operator [] (const index_t & index)
  { return text_[index]; }
  
  // accessors:
  inline const char * text() const
  { return text_; }
  
  inline const size_t & size() const
  { return size_; }
  
  // conversion operator:
  inline operator const char * () const
  { return text_; }
  
  inline static the_text_t pad(const char * str,
                               const size_t width = 0,
                               const char pad_char = ' ',
                               const bool pad_left = true)
  {
    the_text_t txt(str);
    
    if (width > txt.size())
    {
      the_text_t padding;
      padding.fill(pad_char, width - txt.size());
      txt = pad_left ? padding + txt : txt + padding;
    }
    
    return txt;
  }
  
  // helpers:
  template <class number_t>
  static the_text_t number(const number_t & number,
                           const size_t width = 0,
                           const char pad_char = ' ',
                           const bool pad_left = true)
  {
    std::ostringstream os;
    os << number;
    
    std::string str = os.str();
    return pad(str.c_str(), width, pad_char, pad_left);
  }
  
  inline static the_text_t number(const size_t & number,
                                  const size_t width = 0,
                                  const char pad_char = ' ',
                                  const bool pad_left = true)
  {
#ifdef _WIN32
#ifndef snprintf
#define snprintf _snprintf_s
#endif
#endif
    
    static char buffer[256];
    snprintf(buffer, sizeof(buffer), "%llu", (long long unsigned int)(number));
    return pad(buffer, width, pad_char, pad_left);
  }
  
  short int                toShort(bool * ok = 0, int base = 10) const;
  unsigned short int        toUShort(bool * ok = nullptr, int base = 10) const;
  
  int                        toInt(bool * ok = nullptr, int base = 10) const;
  unsigned int                toUInt(bool * ok = nullptr, int base = 10) const;
  
  long int                toLong(bool * ok = nullptr, int base = 10) const;
  unsigned long int        toULong(bool * ok = nullptr, int base = 10) const;
  
  float                        toFloat(bool * ok = nullptr) const;
  double                toDouble(bool * ok = nullptr) const;

  void to_ascii();
  void to_lower();
  void to_upper();
  void fill(const char & c, const size_t size);
  
  void fill(const char & c)
  { fill(c, size_); }
  
  // return true if the given text matches the head/tail of this text:
  bool match_head(const the_text_t & t, bool ignore_case = false) const;
  bool match_tail(const the_text_t & t, bool ignore_case = false) const;
  
  bool match_text(const the_text_t & t,
                  const size_t & index,
                  bool ignore_case = false) const;
  
  // remove leading/tailing white space, replace internal white space
  // with a single space:
  the_text_t simplify_ws() const;
  
  // split the text into a set of tokens, return the number of tokens:
  size_t split(std::vector<the_text_t> & tokens,
               const char & separator,
               const bool & empty_ok = false) const;

  // splitAt splits string into two parts upon a character.
  std::vector<the_text_t> splitAt(const char split_char, 
				  unsigned int n) const;
  
  // count the number of occurrences of a given symbol in the text:
  size_t contains(const char & symbol) const;
  
  // extract a portion of the string:
  void extract(the_text_t & to,
               const size_t & from,
               const size_t & size) const
  {
    assert(from + size < size_ + 1);
    to.assign(&text_[from], size);
  }
  
  inline the_text_t extract(const size_t & from,
                            const size_t & size) const
  {
    the_text_t to;
    extract(to, from, size);
    return to;
  }
  
  inline the_text_t reverse() const
  {
    the_text_t rev(*this);
    for (size_t i = 0; i < size_; i++)
    {
      rev.text_[i] = text_[size_ - i - 1];
    }
    
    return rev;
  }
  
  inline the_text_t cut(const char & separator,
                        size_t f0,
                        size_t f1 = 0) const
  {
    const char sep_str[2] = { separator, '\0' };
    
    std::vector<the_text_t> fields;
    split(fields, separator, true);
    size_t num_fields = fields.size();
    
    if (f1 < f0)
    {
      f1 = f0;
    }
    else if (f1 >= num_fields)
    {
      f1 = num_fields - 1;
    }
    
    the_text_t out;
    for (size_t f = f0; f <= f1; f++)
    {
      out += fields[f];
      if (f + 1 <= f1)
      {
        out += sep_str;
      }
    }
    
    return out;
  }
  
private:
  // the text itself:
  char * text_;
  
  // the length of the text:
  size_t size_;
};

extern std::ostream &
operator << (std::ostream & out, const the_text_t & text);

extern std::istream &
operator >> (std::istream & in, the_text_t & text);

extern std::istream &
getline(std::istream & in, the_text_t & text);

//----------------------------------------------------------------
// to_binary
//
// return a 0 and 1 string representation of a byte
// 
extern the_text_t
to_binary(const unsigned char & byte, bool lsb_first = true);


#endif // itkIRText_h
