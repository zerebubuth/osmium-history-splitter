#ifndef OSMIUM_HISTORY_SPLITTER_UTIL_HPP
#define OSMIUM_HISTORY_SPLITTER_UTIL_HPP

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

namespace hsplitter {
namespace util {

template <typename ConstIterator>
struct pair_snd_iterator {
  typedef ConstIterator const_iterator;
  explicit pair_snd_iterator(const_iterator i) : m_itr(i) {}
  pair_snd_iterator(const pair_snd_iterator &other) : m_itr(other.m_itr) {}

  typedef const uint32_t &reference;

  inline reference operator*() { return m_itr->second; }
  inline pair_snd_iterator &operator++() { ++m_itr; return *this; }
  inline pair_snd_iterator operator++(int) { pair_snd_iterator rv = *this; ++(*this); return rv; }

  inline bool operator==(const pair_snd_iterator &other) const { return m_itr == other.m_itr; }
  inline bool operator!=(const pair_snd_iterator &other) const { return m_itr != other.m_itr; }

private:
  const_iterator m_itr;
};

template <typename Iterator>
struct iter_pair_range {
  typedef std::pair<Iterator, Iterator> value_type;
  value_type m_pair;

  explicit iter_pair_range(value_type v) : m_pair(v) {}

  Iterator begin() const { return m_pair.first; }
  Iterator end()   const { return m_pair.second; }
};

} // namespace util
} // namespace hsplitter

#endif /* OSMIUM_HISTORY_SPLITTER_UTIL_HPP */
