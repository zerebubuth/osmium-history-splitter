#ifndef OSMIUM_HISTORY_SPLITTER_CHUNKED_ARRAY_HPP
#define OSMIUM_HISTORY_SPLITTER_CHUNKED_ARRAY_HPP

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
#include <cassert>

namespace hsplitter {
namespace container {

template <typename OuterIter, typename InnerIter, size_t InnerSize>
struct nested_iterator
  : public std::iterator<std::random_access_iterator_tag,
                         typename InnerIter::value_type> {
  typedef std::iterator<std::random_access_iterator_tag,
                        typename InnerIter::value_type,
                        typename InnerIter::difference_type,
                        typename InnerIter::pointer,
                        typename InnerIter::reference> parent;
  typedef typename parent::pointer pointer;
  typedef typename parent::reference reference;
  typedef typename parent::difference_type difference_type;
  typedef typename parent::value_type value_type;
  typedef nested_iterator<OuterIter, InnerIter, InnerSize> self;

  nested_iterator(OuterIter itr, OuterIter begin, OuterIter end)
    : m_outer(itr), m_outer_begin(begin), m_outer_end(end) {
    if (m_outer != m_outer_end) {
      m_inner = std::begin(*m_outer);
    }
  }

  nested_iterator(const self &) = default;
  nested_iterator(self &&) = default;
  self &operator=(const self &) = default;
  self &operator=(self &&) = default;

  reference operator*() { return *m_inner; }
  pointer operator->() { return m_inner.operator->(); }

  inline self &operator++() {
    if (m_outer != m_outer_end) {
      ++m_inner;
      if (m_inner == std::end(*m_outer)) {
        ++m_outer;
        if (m_outer != m_outer_end) {
          m_inner = std::begin(*m_outer);
        }
      }
    }
    return *this;
  }
  inline self operator++(int) {
    self result = *this;
    ++(*this);
    return result;
  }

  friend void swap(self &a, self &b) {
    std::swap(a.m_outer, b.m_outer);
    std::swap(a.m_outer_begin, b.m_outer_begin);
    std::swap(a.m_outer_end, b.m_outer_end);
    std::swap(a.m_inner, b.m_inner);
  }

  inline bool operator==(const self &other) const {
    return ((m_outer == other.m_outer) &&
            ((m_outer == m_outer_end) || (m_inner == other.m_inner)) &&
            (m_outer_begin == other.m_outer_begin) &&
            (m_outer_end == other.m_outer_end));
  }
  inline bool operator!=(const self &other) const {
    return !operator==(other);
  }

  inline self &operator--() {
    if ((m_inner == std::begin(*m_outer) &&
         m_outer != m_outer_begin) ||
        (m_outer == m_outer_end)) {
      --m_outer;
      m_inner = std::end(*m_outer);
    }
    --m_inner;
    return *this;
  }
  inline self operator--(int) {
    self result = *this;
    ++(*this);
    return result;
  }

  self operator+(difference_type i) const {
    self result = *this;
    result += i;
    return result;
  }
  self operator-(difference_type i) const {
    self result = *this;
    result -= i;
    return result;
  }
  difference_type operator-(const self &other) const {
    if (other > *this) {
      return -(other.operator-(*this));

    } else if (other.m_outer == m_outer_end) {
      // must be end/end comparison, as other <= this.
      return 0;

    } else {
      InnerIter inner = m_inner;
      OuterIter outer = m_outer;

      if (outer == m_outer_end) {
        --outer;
        inner = std::end(*outer);
      }

      if (outer == other.m_outer) {
        return inner - other.m_inner;

      } else {
        difference_type d1 = std::end(*other.m_outer) - other.m_inner;
        difference_type d2 = InnerSize * ((outer - other.m_outer) - 1);
        difference_type d3 = inner - std::begin(*outer);

        return d1 + d2 + d3;
      }
    }
  }

  self &operator+=(difference_type i) {
    if (i < 0) {
      operator-=(-i);

    } else if (i > 0) {
      difference_type left_in_block = std::end(*m_outer) - m_inner;
      assert(left_in_block >= 0);
      if (i < left_in_block) {
        m_inner += i;

      } else {
        i -= left_in_block;
        m_outer += (i / InnerSize) + 1;
        m_inner = std::begin(*m_outer) + (i % InnerSize);
      }
    }

    return *this;
  }
  self &operator-=(difference_type i) {
    if (i < 0) {
      operator+=(-i);

    } else if (i > 0) {
      if (m_outer == m_outer_end) {
        --m_outer;
        m_inner = std::end(*m_outer);
      }

      difference_type left_in_block = m_inner - std::begin(*m_outer);
      assert(left_in_block >= 0);

      if (i < left_in_block) {
        m_inner -= i;

      } else {
        i -= left_in_block;
        m_outer -= (i / InnerSize) + 1;
        m_inner = std::end(*m_outer) - (i % InnerSize);
        if (m_inner == std::end(*m_outer)) {
          ++m_outer;
          m_inner = std::begin(*m_outer);
        }
      }
    }

    return *this;
  }

  inline bool operator<(const self &other) const {
    return (m_outer < other.m_outer ||
            (m_outer == other.m_outer &&
             m_inner < other.m_inner));
  }
  inline bool operator>(const self &other) const {
    return (m_outer > other.m_outer ||
            (m_outer == other.m_outer &&
             m_inner > other.m_inner));
  }
  inline bool operator<=(const self &other) const {
    return !operator>(other);
  }
  inline bool operator>=(const self &other) const {
    return !operator<(other);
  }

  inline reference operator[](difference_type i) {
    self result = *this;
    result += i;
    return *result;
  }

private:
  OuterIter m_outer, m_outer_begin, m_outer_end;
  InnerIter m_inner;
};

template <typename It1, typename It2, size_t S>
inline nested_iterator<It1, It2, S> operator+(std::ptrdiff_t i, const nested_iterator<It1, It2, S> &it) {
  return it + i;
}

template <typename T, size_t chunk_size = 16384>
struct chunked_array {
  typedef std::vector<T> inner_type;
  typedef std::vector<inner_type> outer_type;
  typedef nested_iterator<typename outer_type::iterator, typename inner_type::iterator, chunk_size> iterator;
  typedef nested_iterator<typename outer_type::const_iterator, typename inner_type::const_iterator, chunk_size> const_iterator;
  typedef T value_type;

  inline iterator begin() {
    return iterator(std::begin(m_container), std::begin(m_container), std::end(m_container));
  }
  inline iterator end() {
    return iterator(std::end(m_container), std::begin(m_container), std::end(m_container));
  }
  inline const_iterator begin() const {
    return const_iterator(std::begin(m_container), std::begin(m_container), std::end(m_container));
  }
  inline const_iterator end() const {
    return const_iterator(std::end(m_container), std::begin(m_container), std::end(m_container));
  }
  inline const_iterator cbegin() const { return begin(); }
  inline const_iterator cend() const { return end(); }

  inline void push_back(const T &t) {
    auto &inner = ensure_capacity();
    assert(inner.size() < chunk_size);
    inner.push_back(t);
  }
  inline void push_back(T &&t) {
    auto &inner = ensure_capacity();
    assert(inner.size() < chunk_size);
    inner.push_back(std::move(t));
  }

  inline size_t capacity() const { return chunk_size * m_container.size(); }
  inline size_t size() const {
    if (m_container.empty()) {
      return 0;

    } else {
      return chunk_size * (m_container.size() - 1) + m_container.back().size();
    }
  }
  inline bool empty() { return m_container.empty(); }

  T &back()              { return m_container.back().back(); }
  const T &back() const  { return m_container.back().back(); }
  T &front()             { return m_container.front().front(); }
  const T &front() const { return m_container.front().front(); }

  void clear() { m_container.clear(); }

private:
  outer_type m_container;

  inline inner_type &ensure_capacity() {
    if (m_container.empty() ||
        (m_container.back().size() == chunk_size)) {
      m_container.emplace_back();
      m_container.back().reserve(chunk_size);
    }
    return m_container.back();
  }
};

} // namespace container
} // namespace hsplitter

namespace std {
template <typename It1, typename It2, size_t S>
struct iterator_traits<hsplitter::container::nested_iterator<It1, It2, S> > {
  typedef random_access_iterator_tag iterator_category;
  typedef typename hsplitter::container::nested_iterator<It1, It2, S>::pointer pointer;
  typedef typename hsplitter::container::nested_iterator<It1, It2, S>::reference reference;
  typedef typename hsplitter::container::nested_iterator<It1, It2, S>::value_type value_type;
  typedef typename hsplitter::container::nested_iterator<It1, It2, S>::difference_type difference_type;
};
} // namespace std

#endif /* OSMIUM_HISTORY_SPLITTER_CHUNKED_ARRAY_HPP */
