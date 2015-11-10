#include <osmium/io/any_input.hpp>
#include <osmium/io/any_output.hpp>
#include <osmium/io/any_compression.hpp>
#include <osmium/io/input_iterator.hpp>
#include <osmium/thread/pool.hpp>
#include <osmium/thread/util.hpp>
#include <osmium/util/memory_mapping.hpp>

#include "splitter.hpp"
#include "tile_grid.hpp"

#include <boost/format.hpp>

#include <iostream>
#include <map>
#include <set>
#include <cassert>
#include <fstream>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

using hsplitter::tile_file;
using hsplitter::tile_grid;

namespace hsplitter {
size_t g_evictions = 0, g_flushes = 0;
} // namespace hsplitter

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

  tile_map_subarray()
    : m_unsorted_count(0), m_array() {}
  tile_map_subarray(tile_map_subarray &&t)
    : m_unsorted_count(t.m_unsorted_count),
      m_array(std::move(t.m_array)) {
  }
  tile_map_subarray(const tile_map_subarray &) = delete;

  tile_map_subarray &operator=(tile_map_subarray &&t) {
    m_unsorted_count = t.m_unsorted_count;
    m_array = std::move(t.m_array);
    return *this;
  }
  tile_map_subarray &operator=(const tile_map_subarray &) = delete;

  void insert(key_type k, mapped_type v) {
    auto p = std::make_pair(k, v);
    if (m_array.empty() || p != m_array.back()) {
      if ((m_unsorted_count > 0) ||
          (!m_array.empty() && (m_array.back() > p))) {
        ++m_unsorted_count;
      }
      m_array.push_back(p);

      if (m_unsorted_count > 16777216) {
        size_t before = m_array.size();
        sort_array();
        size_t after = m_array.size();
        std::cerr << "too many unsorted, sorting! " << before << " -> " << after << std::endl;
      }
    }
  }

  void sort_array() const {
    //std::cerr << "sorting!" << std::endl;
    std::sort(m_array.begin(), m_array.end());
    auto itr = std::unique(m_array.begin(), m_array.end());
    m_array.erase(itr, m_array.end());
    m_unsorted_count = 0;
  }

  iter_pair_range<const_iterator> equal_range(key_type k) const {
    if (m_unsorted_count > 0) {
      sort_array();
    }
    const auto &lb = std::lower_bound(m_array.begin(), m_array.end(), std::make_pair(k, mapped_type(0)));
    const auto &ub = std::upper_bound(lb, m_array.end(), std::make_pair(k, std::numeric_limits<mapped_type>::max()));
    return iter_pair_range<const_iterator>(std::make_pair(const_iterator(lb), const_iterator(ub)));
  }

  mutable uint64_t m_unsorted_count;
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
    uint64_t idx = k >> 24;
    uint32_t prt = k & ((uint64_t(1) << 24) - 1);
    //std::cerr << "insert(" << k << ", " << v << ") -> " << idx << ", " << prt << std::endl;
    if (idx >= m_array.size()) {
      m_array.resize(idx + 1);
    }
    tile_map_subarray &sub = m_array[idx];
    sub.insert(prt, v);
  }

  iter_pair_range<const_iterator> equal_range(key_type k) const {
    assert(k >= 0);
    uint32_t idx = k >> 24;
    uint32_t prt = k & ((uint64_t(1) << 24) - 1);
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
    node_tiles = hsplitter::tiles_for_nodes<tileset_t, 14>(itr, input_range.end());

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

  tile_grid grid(1000, 100*1024);
  std::unordered_set<hsplitter::tile_t> all_tiles;
  osmium::io::Header header;
  size_t count = 0;

  std::cerr << "Second pass..." << std::endl;
  { // second pass through file
    osmium::io::File infile(argv[1]);
    osmium::io::Reader reader(infile, osmium::osm_entity_bits::nwr);
    auto input_range = osmium::io::make_input_iterator_range<osmium::OSMObject>(reader);

    header = reader.header();

    for (const auto &obj : input_range) {
      ++count;
      if ((count & 0xffff) == 0) {
        std::cout << obj.type() << " " << obj.id() << std::endl;
      }
      std::set<hsplitter::tile_t> tiles;

      if (obj.type() == osmium::item_type::node) {
        insert_tiles_from<osmium::Node>(obj, node_tiles, tiles);
        insert_tiles_from<osmium::Node>(obj, extra_node_tiles, tiles);

        const auto &node = reinterpret_cast<const osmium::Node &>(obj);
        for (auto tile_id : tiles) {
          auto &file = grid.tile(tile_id);
          file.write(node);
        }

      } else if (obj.type() == osmium::item_type::way) {
        insert_tiles_from<osmium::Way>(obj, way_tiles, tiles);

        const auto &way = reinterpret_cast<const osmium::Way &>(obj);
        for (auto tile_id : tiles) {
          auto &file = grid.tile(tile_id);
          file.write(way);
        }

      } else if (obj.type() == osmium::item_type::relation) {
        insert_tiles_from<osmium::Relation>(obj, rel_tiles, tiles);
        insert_tiles_from<osmium::Relation>(obj, extra_rel_tiles, tiles);

        const auto &rel = reinterpret_cast<const osmium::Relation &>(obj);
        for (auto tile_id : tiles) {
          auto &file = grid.tile(tile_id);
          file.write(rel);
        }
      }

      for (auto tile_id : tiles) {
        all_tiles.insert(tile_id);
      }
    }

    reader.close();
  }
  grid.close();

  std::cerr << hsplitter::g_evictions << " evictions, " << hsplitter::g_flushes << " flushes." << std::endl;

  std::cerr << "Converting buffers..." << std::endl;
  count = 0;

  //std::vector<std::future<void> > threads;
  for (auto tile_id : all_tiles) {
    ++count;
    if ((count & 0xff) == 0) {
      std::cerr << count << " / " << all_tiles.size() << std::endl;
    }

    //auto thread = osmium::thread::Pool::instance().submit([tile_id, &header] () {
    // TODO: configurable, get location for tmp from grid
    std::string input_filename = (boost::format("tiles/%1%.buf") % tile_id).str();
    std::string output_filename = (boost::format("tiles/%1%.osh.pbf") % tile_id).str();

    int fd = ::open(input_filename.c_str(), O_RDONLY);
    if (fd == -1) {
      throw std::runtime_error("Unable to open file.");
    }

    struct ::stat status;
    if (::fstat(fd, &status) == -1) {
      throw std::runtime_error("Unable to get file status.");
    }

    osmium::util::MemoryMapping mapping(status.st_size, osmium::util::MemoryMapping::mapping_mode::readonly, fd, 0);
    osmium::memory::Buffer buffer(mapping.get_addr<unsigned char>(), status.st_size);
    osmium::io::File outfile(output_filename, "osh.pbf");
    osmium::io::Writer writer(outfile, header, osmium::io::overwrite::allow);
    writer(std::move(buffer));
    writer.close();
    mapping.unmap();
    ::close(fd);
    // TODO: rm buffer file
    //  });
    //threads.push_back(std::move(thread));
  }

  // for (auto &thread : threads) {
  //   osmium::thread::wait_until_done(thread);
  // }

  return 0;
}
