#include <osmium/io/any_input.hpp>
#include <osmium/io/any_output.hpp>
#include <osmium/io/any_compression.hpp>
#include <osmium/io/input_iterator.hpp>
#include <osmium/thread/pool.hpp>
#include <osmium/thread/util.hpp>
#include <osmium/util/memory_mapping.hpp>

#include "splitter.hpp"

#include <boost/format.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>

#include <iostream>
#include <map>
#include <set>
#include <cassert>
#include <fstream>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

namespace bmi = boost::multi_index;

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

size_t g_evictions = 0, g_flushes = 0;

template <typename T>
struct builder_trait {};

template <>
struct builder_trait<osmium::Node> {
  typedef osmium::builder::NodeBuilder type;
};

template <>
struct builder_trait<osmium::Way> {
  typedef osmium::builder::WayBuilder type;
};

template <>
struct builder_trait<osmium::Relation> {
  typedef osmium::builder::RelationBuilder type;
};

struct tile_file {
  typedef osmium::memory::Buffer buffer;
  typedef std::unique_ptr<buffer> buffer_ptr;
  typedef hsplitter::tile_t tile_t;

  tile_file()
    : m_id(0), m_buffer() {
  }

  tile_file(tile_t id, size_t capacity)
    : m_id(id),
      m_buffer(new buffer(capacity, osmium::memory::Buffer::auto_grow::yes)) {
    set_callback();
  }

  // tile_file(tile_t id, buffer_ptr &&buffer)
  //   : m_id(id), m_buffer(std::move(buffer)) {
  //   set_callback();
  // }

  tile_file(tile_file &&tf) : m_id(tf.m_id), m_buffer() {
    swap_buffer(tf.m_buffer);
  }
  //tile_file(tile_file &&) = delete;
  tile_file(const tile_file &) = delete;
  const tile_file &operator=(tile_file &&tf) {
    m_id = tf.m_id;
    swap_buffer(tf.m_buffer);
    return *this;
  }
  //const tile_file &operator=(tile_file &&tf) = delete;
  const tile_file &operator=(const tile_file &) = delete;

  void write(const osmium::OSMObject &obj) {
    assert(!empty());
    m_buffer->add_item(obj);
    m_buffer->commit();

    // auto itr = m_buffer->get_iterator<osmium::OSMObject>(offset);
    // assert(itr->id() == obj.id());
    // assert(itr->version() == obj.version());
    // assert(itr->next() == (m_buffer->data() + m_buffer->committed()));
  }

  void flush() {
    assert(!empty());
    static_flush(m_id, *m_buffer);
  }

  static void static_flush(hsplitter::tile_t id, osmium::memory::Buffer &buffer) {
    assert(buffer);
    if (buffer.begin() != buffer.end()) {
      // TODO: configurable, default to $TMPDIR
      std::string filename = (boost::format("tiles/%1%.buf") % id).str();
      std::ofstream file(filename, std::ios::binary | std::ios::ate | std::ios::app);
      file.write(reinterpret_cast<const char *>(buffer.data()), buffer.committed());
      buffer.clear();
      //assert(buffer.capacity() == 100*1024);
      std::fill(buffer.data(), buffer.data() + buffer.capacity(), 0);

      ++g_flushes;
    }
  }

  void swap(tile_file &tf) {
    std::swap(m_id, tf.m_id);
    swap_buffer(tf.m_buffer);
  }

  bool empty() const {
    return !m_buffer;
  }

  tile_t m_id;

private:
  buffer_ptr m_buffer;

  void swap_buffer(buffer_ptr &other) {
    if (m_buffer) {
      m_buffer->set_full_callback([](osmium::memory::Buffer&){});
    }
    if (other) {
      other->set_full_callback([](osmium::memory::Buffer&){});
    }
    std::swap(m_buffer, other);
    set_callback();
  }

  void set_callback() {
    const hsplitter::tile_t id = m_id;
    if (m_buffer) {
      m_buffer->set_full_callback([id](osmium::memory::Buffer &b) {
          static_flush(id, b);
        });
    }
  }
};

struct tile_grid {
  typedef bmi::multi_index_container<
    tile_file,
    bmi::indexed_by<
      bmi::sequenced<>,
      bmi::hashed_unique<bmi::member<tile_file, hsplitter::tile_t, &tile_file::m_id> >
      >
    > mru_tile_set;

  tile_grid(size_t n_tiles, size_t capacity) {
    for (size_t i = 0; i < n_tiles; ++i) {
      tile_file tf(0, capacity);
      m_free_tiles.emplace_back(std::move(tf));
    }
  }

  tile_file &tile(hsplitter::tile_t id) {
    auto &tile_idx = m_open_tiles.get<1>();
    auto itr = tile_idx.find(id);

    //std::cerr << m_open_tiles.size() << " open, " << m_free_tiles.size() << " free, " << g_evictions << " evictions, " << g_flushes << " flushes." << std::endl;
    if (itr == tile_idx.end()) {
      if (m_free_tiles.empty()) {
        evict_lru_tile();
      }

      tile_file tf;
      tf.swap(m_free_tiles.front());
      m_free_tiles.pop_front();
      tf.m_id = id;
      auto pair = tile_idx.emplace(std::move(tf));
      itr = pair.first;
    }

    m_open_tiles.relocate(m_open_tiles.begin(), m_open_tiles.project<0>(itr));
    assert(itr->m_id == id);
    return const_cast<tile_file &>(*itr);
  }

  void evict_lru_tile() {
    assert(m_open_tiles.size() > 0);

    auto &tile = const_cast<tile_file &>(m_open_tiles.back());
    tile.flush();
    m_free_tiles.emplace_front();
    m_free_tiles.front().swap(tile);
    assert(m_open_tiles.back().empty());
    m_open_tiles.pop_back();

    ++g_evictions;
  }

  void close() {
    for (auto &tile : m_open_tiles) {
      const_cast<tile_file &>(tile).flush();
    }
  }

private:
  std::list<tile_file> m_free_tiles;
  mru_tile_set m_open_tiles;
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

  std::cerr << g_evictions << " evictions, " << g_flushes << " flushes." << std::endl;

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
