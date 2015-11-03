#include <osmium/io/any_input.hpp>
#include <osmium/io/any_compression.hpp>
#include <osmium/io/input_iterator.hpp>

#include "splitter.hpp"

#include <iostream>
#include <map>
#include <set>
#include <cassert>


namespace {
template <typename ObjType, typename TilesType, typename Container>
void insert_tiles_from(const osmium::OSMObject &o, const TilesType &tiles,
                       Container &target) {
  const auto &obj = static_cast<const ObjType &>(o);
  auto range = tiles.equal_range(obj.id());
  if (range.begin() != range.end()) {
    target.insert(range.begin(), range.end());
  }
}

struct tile_file {
  hsplitter::tile_t id;
  void write(const osmium::OSMEntity &entity) {
    const auto &obj = static_cast<const osmium::OSMObject &>(entity);
    std::cout << id << "\t"
              << obj.type() << "\t"
              << obj.id() << "\n";
  }
};

struct tile_grid {
  tile_file &tile(hsplitter::tile_t id) {
    m_file.id = id;
    return m_file;
  }
  tile_file m_file;
};

struct pair_snd_iterator {
  typedef std::vector<std::pair<uint32_t, uint32_t> >::const_iterator const_iterator;
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

struct tile_map_subarray {
  typedef uint32_t key_type;
  typedef uint32_t mapped_type;
  typedef pair_snd_iterator const_iterator;
  typedef std::vector<std::pair<uint32_t, uint32_t> > cont_t;

  tile_map_subarray() {}
  tile_map_subarray(tile_map_subarray &&t)
    : m_is_sorted(true), m_last_key(0), m_array(std::move(t.m_array)) {
  }
  tile_map_subarray(const tile_map_subarray &) = delete;

  tile_map_subarray &operator=(tile_map_subarray &&t) {
    m_is_sorted = t.m_is_sorted;
    m_last_key = t.m_last_key;
    m_array = std::move(t.m_array);
    return *this;
  }
  tile_map_subarray &operator=(const tile_map_subarray &) = delete;

  void insert(key_type k, mapped_type v) {
    if (m_is_sorted && m_last_key > k) { m_is_sorted = false; }
    auto p = std::make_pair(k, v);
    if (m_array.empty() || p != m_array.back()) {
      m_array.push_back(p);
      m_last_key = k;
    }
  }

  iter_pair_range<const_iterator> equal_range(key_type k) const {
    if (!m_is_sorted) {
      std::sort(m_array.begin(), m_array.end());
      auto itr = std::unique(m_array.begin(), m_array.end());
      m_array.erase(itr, m_array.end());
      m_is_sorted = true;
      m_last_key = m_array.empty() ? 0 : m_array.back().first;
    }
    const auto &lb = std::lower_bound(m_array.begin(), m_array.end(), std::make_pair(k, mapped_type(0)));
    const auto &ub = std::upper_bound(lb, m_array.end(), std::make_pair(k, std::numeric_limits<mapped_type>::max()));
    return iter_pair_range<const_iterator>(std::make_pair(const_iterator(lb), const_iterator(ub)));
  }

  mutable bool m_is_sorted;
  mutable uint64_t m_last_key;
  mutable cont_t m_array;
};

struct tile_map_array {
  typedef int64_t key_type;
  typedef uint32_t mapped_type;
  typedef pair_snd_iterator const_iterator;

  tile_map_array() {}
  tile_map_array(tile_map_array &&t) : m_array(std::move(t.m_array)) {}
  tile_map_array(const tile_map_array &) = delete;

  tile_map_array &operator=(tile_map_array &&t) { m_array = std::move(t.m_array); return *this; }
  tile_map_array &operator=(const tile_map_array &) = delete;

  void insert(key_type k, mapped_type v) {
    assert(k >= 0);
    uint32_t idx = k >> 32;
    uint32_t prt = k & ((uint64_t(1) << 32) - 1);
    if (idx >= m_array.size()) {
      m_array.resize(idx + 1);
    }
    tile_map_subarray &sub = m_array[idx];
    sub.insert(prt, v);
  }

  iter_pair_range<const_iterator> equal_range(key_type k) const {
    assert(k >= 0);
    uint32_t idx = k >> 32;
    uint32_t prt = k & ((uint64_t(1) << 32) - 1);
    if (idx >= m_array.size()) {
      return iter_pair_range<const_iterator>(std::make_pair(pair_snd_iterator(m_empty.end()), pair_snd_iterator(m_empty.end())));

    } else {
      const tile_map_subarray &sub = m_array[idx];
      return sub.equal_range(prt);
    }
  }

  std::vector<tile_map_subarray> m_array;
  static const tile_map_subarray::cont_t m_empty;
};

const tile_map_subarray::cont_t tile_map_array::m_empty;

} // anonymous namespace

int main(int argc, char *argv[]) {
  using tileset_t = tile_map_array;

  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " INFILE\n";
    return 1;
  }

  tileset_t node_tiles, way_tiles, extra_node_tiles, rel_tiles, extra_rel_tiles;

  std::cerr << "First pass..." << std::endl;
  { // first pass through file
    osmium::io::File infile(argv[1]);
    osmium::io::Reader reader(infile, osmium::osm_entity_bits::nwr);
    auto input_range = osmium::io::make_input_iterator_range<osmium::OSMObject>(reader);

    auto itr = input_range.begin();
    std::cerr << " >> nodes" << std::endl;
    node_tiles = hsplitter::tiles_for_nodes<tileset_t>(itr, input_range.end());

    std::cerr << " >> ways" << std::endl;
    auto pair_way = hsplitter::tiles_for_ways(itr, input_range.end(), node_tiles);
    way_tiles = std::move(pair_way.first);
    extra_node_tiles = std::move(pair_way.second);

    std::cerr << " >> relations" << std::endl;
    auto pair_rel = hsplitter::tiles_for_relations(itr, input_range.end(), node_tiles,
                                                   way_tiles, extra_node_tiles);
    rel_tiles = std::move(pair_rel.first);
    extra_rel_tiles = std::move(pair_rel.second);

    reader.close();
  }

  tile_grid grid;

  std::cerr << "Second pass..." << std::endl;
  { // second pass through file
    osmium::io::File infile(argv[1]);
    osmium::io::Reader reader(infile, osmium::osm_entity_bits::nwr);
    auto input_range = osmium::io::make_input_iterator_range<osmium::OSMObject>(reader);

    //osmium::io::Header header = reader.header();

    for (const auto &entity : input_range) {
      std::set<hsplitter::tile_t> tiles;

      if (entity.type() == osmium::item_type::node) {
        insert_tiles_from<osmium::Node>(entity, node_tiles, tiles);
        insert_tiles_from<osmium::Node>(entity, extra_node_tiles, tiles);

      } else if (entity.type() == osmium::item_type::way) {
        insert_tiles_from<osmium::Way>(entity, way_tiles, tiles);

      } else if (entity.type() == osmium::item_type::relation) {
        insert_tiles_from<osmium::Relation>(entity, rel_tiles, tiles);
        insert_tiles_from<osmium::Relation>(entity, extra_rel_tiles, tiles);
      }

      for (auto tile_id : tiles) {
        auto &file = grid.tile(tile_id);
        file.write(entity);
      }
    }

    reader.close();
  }

  return 0;
}
