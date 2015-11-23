#ifndef OSMIUM_HISTORY_SPLITTER_MMAPPED_ARRAY_HPP
#define OSMIUM_HISTORY_SPLITTER_MMAPPED_ARRAY_HPP

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

#include <vector>
#include <memory>
#include <algorithm>
#include <iterator>

#include <boost/optional.hpp>
#include <boost/filesystem.hpp>

#include "chunked_array.hpp"
#include "util.hpp"

namespace hsplitter {

template <typename T>
struct mmap_ptr {
  typedef T& reference;
  typedef T* pointer;

  mmap_ptr() : m_ptr(nullptr), m_len(0) {}
  mmap_ptr(const std::string &filename);
  mmap_ptr(mmap_ptr<T> &&p) : m_ptr(nullptr), m_len(0) { swap(p); }
  mmap_ptr(const mmap_ptr<T> &) = delete;

  ~mmap_ptr();

  mmap_ptr &operator=(mmap_ptr<T> &&p) { swap(p); return *this; }
  mmap_ptr &operator=(const mmap_ptr<T> &) = delete;

  void swap(mmap_ptr &other) {
    std::swap(m_ptr, other.m_ptr);
    std::swap(m_len, other.m_len);
  }

  reference operator*() const { return *reinterpret_cast<pointer>(m_ptr); }
  pointer operator->() const { return reinterpret_cast<pointer>(m_ptr); }

  explicit operator bool() const { return m_ptr != nullptr; }

  pointer begin() const { return reinterpret_cast<pointer>(m_ptr); }
  pointer end() const { return reinterpret_cast<pointer>(m_ptr) + (m_len / sizeof(T)); }

private:
  void *m_ptr;
  size_t m_len;
};

template <typename T>
mmap_ptr<T>::mmap_ptr(const std::string &filename)
  : m_ptr(nullptr)
  , m_len(0) {
  int fd = open(filename.c_str(), O_RDONLY | O_NOATIME);
  if (fd == -1) {
    throw std::runtime_error((boost::format("Unable to open(%1%): %2%")
                              % filename % strerror(errno)).str());
  }

  struct stat file_status;
  int status = fstat(fd, &file_status);
  if (status == -1) {
    throw std::runtime_error((boost::format("Unable to fstat(%1%): %2%")
                              % filename % strerror(errno)).str());
  }

  assert(file_status.st_size > 0);

  void *ptr = mmap(nullptr, file_status.st_size, PROT_READ, MAP_SHARED, fd, 0);
  if (ptr == MAP_FAILED) {
    throw std::runtime_error((boost::format("Unable to mmap(%1%): %2%")
                              % filename % strerror(errno)).str());
  }

  status = close(fd);
  if (status == -1) {
    throw std::runtime_error((boost::format("Unable to close(%1%): %2%")
                              % filename % strerror(errno)).str());
  }

  m_ptr = ptr;
  m_len = file_status.st_size;
}

template <typename T>
mmap_ptr<T>::~mmap_ptr() {
  if (m_ptr != nullptr) {
    int status = munmap(m_ptr, m_len);
    if (status != 0) {
      std::cerr << "Failed to munmap pointer: " << m_ptr << ": " << strerror(errno) << std::endl;
    }
    m_ptr = nullptr;
    m_len = 0;
  }
}

struct tmpfile {
  tmpfile() : m_filename() {}
  tmpfile(tmpfile &&t) : m_filename(std::move(t.m_filename)) {}
  tmpfile(const tmpfile &) = delete;

  ~tmpfile() {
    if (m_filename) {
      boost::filesystem::remove(*m_filename);
      m_filename = boost::none;
    }
  }

  tmpfile &operator=(tmpfile &&t) { m_filename = std::move(t.m_filename); return *this; }
  tmpfile &operator=(const tmpfile &) = delete;
  tmpfile &operator=(const std::string &filename) { m_filename = filename; return *this; }

  std::string filename() const {
    assert(bool(m_filename));
    return *m_filename;
  }

private:
  boost::optional<std::string> m_filename;
};

typedef util::pair_snd_iterator<const std::pair<uint32_t, uint32_t> *> pair_snd_iterator;

struct mmapped_subarray {
  typedef uint32_t key_type;
  typedef uint32_t mapped_type;
  typedef pair_snd_iterator const_iterator;
  typedef std::pair<uint32_t, uint32_t> pair_t;

  mmapped_subarray()
    : m_fd(-1), m_array(), m_mmap(), m_empty(false), m_tmpfile(), m_last() {
    m_array.reserve(262144);
  }
  mmapped_subarray(mmapped_subarray &&t)
    : m_fd(t.m_fd)
    , m_array(std::move(t.m_array))
    , m_mmap(std::move(t.m_mmap))
    , m_empty(t.m_empty)
    , m_tmpfile(std::move(t.m_tmpfile))
    , m_last(std::move(t.m_last)) {
  }
  mmapped_subarray(const mmapped_subarray &) = delete;

  mmapped_subarray &operator=(mmapped_subarray &&t) {
    m_fd = t.m_fd;
    m_array = std::move(t.m_array);
    m_mmap = std::move(t.m_mmap);
    m_empty = t.m_empty;
    m_tmpfile = std::move(t.m_tmpfile);
    m_last = std::move(t.m_last);
    return *this;
  }
  mmapped_subarray &operator=(const mmapped_subarray &) = delete;

  void insert(key_type k, mapped_type v) {
    assert(!bool(m_mmap));
    if (m_fd == -1) { open_tmpfile(); }
    auto p = std::make_pair(k, v);

    // sanity check - must append _in_order_
    assert(!bool(m_last) || (*m_last < p));
    m_last = p;

    if (m_array.size() == m_array.capacity()) {
      flush_array();
    }

    m_array.emplace_back(std::move(p));
  }

  util::iter_pair_range<const_iterator> equal_range(key_type k) const {
    assert(bool(m_mmap) || m_empty);
    if (m_mmap) {
      const pair_t *begin = m_mmap.begin();
      const pair_t *end = m_mmap.end();
      const auto &lb = std::lower_bound(begin, end, std::make_pair(k, mapped_type(0)));
      const auto &ub = std::upper_bound(lb, end, std::make_pair(k, std::numeric_limits<mapped_type>::max()));
      return util::iter_pair_range<const_iterator>(std::make_pair(const_iterator(lb), const_iterator(ub)));
    } else {
      return util::iter_pair_range<const_iterator>(std::make_pair(const_iterator(nullptr), const_iterator(nullptr)));
    }
  }

  void freeze() {
    assert(!bool(m_mmap));
    if (m_fd >= 0) {
      flush_array();

      m_mmap = mmap_ptr<pair_t>(m_tmpfile.filename());

      int status = close(m_fd);
      if (status == -1) {
        throw std::runtime_error("Unable to close temporary file.");
      }
      m_fd = -1;

    } else {
      m_empty = true;
    }
  }

private:
  int m_fd;
  std::vector<pair_t> m_array;
  mmap_ptr<pair_t> m_mmap;
  bool m_empty;
  tmpfile m_tmpfile;
  boost::optional<pair_t> m_last;

  void open_tmpfile() {
    const char *tmpdir = ::getenv("TMPDIR");
    std::ostringstream ostr;
    if (tmpdir == nullptr) {
      ostr << "/tmp";
    } else {
      ostr << tmpdir;
    }
    ostr << "/splitter_mmapped_subarray_XXXXXX";
    std::string filename = ostr.str();
    char *buf = (char *)::alloca(filename.size() + 1);
    strncpy(buf, filename.c_str(), filename.size() + 1);
    int fd = ::mkstemp(buf);
    if (fd == -1) {
      throw std::runtime_error((boost::format("Unable to make temporary file to "
                                              "write subarray: %1%")
                                % strerror(errno)).str());
    }
    filename = buf;
    m_tmpfile = filename;
    m_fd = fd;
  }

  void flush_array() {
    assert(m_fd >= 0);
    ssize_t sz = sizeof(pair_t) * m_array.size();
    ssize_t len = ::write(m_fd, m_array.data(), sz);
    assert(len == sz);
    m_array.clear();
  }
};

struct mmapped_array {
  typedef int64_t key_type;
  typedef uint32_t mapped_type;
  typedef pair_snd_iterator const_iterator;
  static const int shift_size = 32;

  mmapped_array() {}
  mmapped_array(mmapped_array &&t) : m_array(std::move(t.m_array)) {}
  mmapped_array(const mmapped_array &) = delete;

  mmapped_array &operator=(mmapped_array &&t) { m_array = std::move(t.m_array); return *this; }
  mmapped_array &operator=(const mmapped_array &) = delete;

  void insert(key_type k, mapped_type v) {
    assert(k >= 0);
    uint64_t idx = k >> shift_size;
    uint32_t prt = k & ((uint64_t(1) << shift_size) - 1);
    //std::cerr << "insert(" << k << ", " << v << ") -> " << idx << ", " << prt << std::endl;
    if (idx >= m_array.size()) {
      m_array.resize(idx + 1);
    }
    mmapped_subarray &sub = m_array[idx];
    sub.insert(prt, v);
  }

  void freeze() {
    for (mmapped_subarray &sub : m_array) {
      sub.freeze();
    }
  }

  util::iter_pair_range<const_iterator> equal_range(key_type k) const {
    assert(k >= 0);
    uint32_t idx = k >> shift_size;
    uint32_t prt = k & ((uint64_t(1) << shift_size) - 1);
    if (idx >= m_array.size()) {
      return util::iter_pair_range<const_iterator>(std::make_pair(pair_snd_iterator(nullptr), pair_snd_iterator(nullptr)));

    } else {
      const mmapped_subarray &sub = m_array[idx];
      return sub.equal_range(prt);
    }
  }

  std::vector<mmapped_subarray> m_array;
};

} // namespace hsplitter

#endif /* OSMIUM_HISTORY_SPLITTER_MMAPPED_ARRAY_HPP */
