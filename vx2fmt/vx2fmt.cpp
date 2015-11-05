/*

This file is part of the Osmium History Splitter

Copyright 2015 Matt Amos <zerebubuth@gmail.com>.

Boost Software License - Version 1.0 - August 17th, 2003

Permission is hereby granted, free of charge, to any person or organization
obtaining a copy of the software and accompanying documentation covered by
this license (the "Software") to use, reproduce, display, distribute,
execute, and transmit the Software, and to prepare derivative works of the
Software, and to permit third-parties to whom the Software is furnished to
do so, all subject to the following:

The copyright notices in the Software and this entire statement, including
the above license grant, this restriction and the following disclaimer,
must be included in all copies of the Software, in whole or in part, and
all derivative works of the Software, unless such copies or derivative
works are solely in the form of machine-executable object code generated by
a source language processor.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.

*/

#include "vx2fmt.hpp"

#include <cstdint>
#include <cstring>
#include <stdexcept>

namespace vx2fmt {

const std::string header_magic("VX2FMT");

struct uvarint {
  uvarint() : m_value(0) {}
  explicit uvarint(uint64_t v) : m_value(v) {}

  inline bool dump(std::ostream &out) const {
    uint8_t buffer[10];
    uint64_t v = m_value;

    uint8_t *ptr = buffer;
    while (v > 127) {
      *ptr++ = uint8_t(128 | (v & 127));
      v >>= 7;
    }
    *ptr++ = uint8_t(v);

    out.write(reinterpret_cast<const char *>(&buffer[0]), std::distance(buffer, ptr));
    return out.good();
  }

  inline bool load(std::istream &in) {
    uint8_t c = 0;
    uint64_t value = 0;

    do {
      in.read(reinterpret_cast<char *>(&c), 1);
      if (!in.good()) { return false; }

      value = (value << 7) | uint64_t(c & 127);
    } while (c > 127);

    m_value = value;
    return true;
  }

  inline uint64_t value() const { return m_value; }

private:
  uint64_t m_value;
};

inline std::ostream &operator<<(std::ostream &out, const uvarint &i) {
  bool ok = i.dump(out);
  if (!ok) { throw std::runtime_error("Write to stream failed."); }
  return out;
}

inline std::istream &operator>>(std::istream &in, uvarint &i) {
  bool ok = i.load(in);
  if (!ok) { throw std::runtime_error("Failed to read from stream."); }
  return in;
}

struct svarint : public uvarint {
  static inline uint64_t encode(int64_t v) {
    if (v < 0) {
      return (uint64_t(-v) << 1) | uint64_t(1);
    } else {
      return uint64_t(v) << 1;
    }
  }

  static inline int64_t decode(uint64_t v) {
    if (v & 1) {
      return -int64_t(v >> 1);
    } else {
      return int64_t(v >> 1);
    }
  }

  svarint() : uvarint() {}
  svarint(int64_t v) : uvarint(encode(v)) {}

  inline int64_t value() const { return decode(uvarint::value()); }
};

inline std::ostream &operator<<(std::ostream &out, const svarint &i) {
  bool ok = i.dump(out);
  if (!ok) { throw std::runtime_error("Write to stream failed."); }
  return out;
}

inline std::istream &operator>>(std::istream &in, svarint &i) {
  bool ok = i.load(in);
  if (!ok) { throw std::runtime_error("Failed to read from stream."); }
  return in;
}

template <int PrefixBits>
inline bool write_string_prefix(std::ostream &out, const std::string &s) {
  uvarint sz((uint64_t(s.size()) << PrefixBits) |
             ((uint64_t(1) << PrefixBits) - 1));
  bool ok = sz.dump(out);
  if (ok) { out.write(s.data(), s.size()); }
  return ok;
}

template <int PrefixBits>
inline bool read_string_prefix(std::istream &in, std::string &s, int &prefix) {
  uvarint sz;
  bool ok = sz.load(in);
  if (ok) {
    size_t bytes = sz.value() >> PrefixBits;
    s.resize(bytes);
    in.read(&s[0], bytes);
    prefix = int(sz.value() & ((1 << PrefixBits) - 1));
  }
  return ok;
}

bool write_string(std::ostream &out, const std::string &s) {
  return write_string_prefix<0>(out, s);
}

bool read_string(std::istream &in, std::string &s) {
  int dummy = 0;
  return read_string_prefix<0>(in, s, dummy);
}

bool write_header(std::ostream &out, const dict_t &headers) {
  out.write(header_magic.data(), header_magic.size());
  if (!out) { return false; }

  return write_tags(out, headers);
}

bool read_header(std::istream &in, dict_t &headers) {
  char magic[6];
  in.read(magic, 6);
  if (!in) { return false; }
  if (strncmp(magic, header_magic.c_str(), 6) != 0) { return false; }

  return read_tags(in, headers);
}

bool write_tags(std::ostream &out, const dict_t &tags) {
  bool ok = true;

  for (const auto &val : tags) {
    ok = write_string_prefix<1>(out, val.first);
    if (!ok) { return ok; }

    ok = write_string(out, val.second);
    if (!ok) { return ok; }
  }

  ok = uvarint().dump(out);
  return ok;
}

bool read_tags(std::istream &in, dict_t &tags) {
  bool ok = true;
  int flag = 0;
  std::string key, val;

  while ((ok = read_string_prefix<1>(in, key, flag)) && (flag == 1)) {
    ok = read_string(in, val);

    if (ok) {
      tags.insert(std::make_pair(key, val));
    }
  }

  return ok;
}

} // namespace vx2fmt
