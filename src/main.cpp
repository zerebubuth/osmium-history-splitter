#include <osmium/io/any_input.hpp>
#include <osmium/io/any_output.hpp>
#include <osmium/io/any_compression.hpp>
#include <osmium/io/input_iterator.hpp>
#include <osmium/thread/pool.hpp>
#include <osmium/thread/util.hpp>
#include <osmium/util/memory_mapping.hpp>

#include "splitter.hpp"
#include "tile_grid.hpp"
#include "tile_map_array.hpp"

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
const tile_map_subarray::cont_t tile_map_array::m_empty;
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

template <typename T>
struct histogram {
  histogram() : m_counts(sizeof(T) * CHAR_BIT + 1, 0) {}
  void insert(T v) {
    static const uint64_t b[] = {0x2ull, 0xCull, 0xF0ull, 0xFF00ull, 0xFFFF0000ull, 0xFFFFFFFF00000000ull};
    static const uint64_t S[] = {1, 2, 4, 8, 16, 32};
    if (v == 0) {
      ++m_counts[0];
    } else {
      T r = 0;
      for (int i = 5; i >= 0; --i) {
        if (v & b[i]) {
          v >>= S[i];
          r |= S[i];
        }
      }
      ++r;
      assert(r < m_counts.size());
      ++m_counts[r];
    }
  }
  std::vector<size_t> m_counts;
};

template <typename T>
std::ostream &operator<<(std::ostream &out, const histogram<T> &h) {
  //out << "{";
  for (size_t i = 0; i < h.m_counts.size(); ++i) {
    if (i > 0) {
      out << ", ";
    }
    out << h.m_counts[i];
  }
  //out << "}";
  return out;
}

void print_histograms(const hsplitter::tile_map_array &array, std::ostream &out) {
  uint64_t last_key = 0;
  uint32_t last_val = 0;
  size_t inner_count = 0;
  histogram<uint64_t> keys_histogram;
  histogram<uint32_t> values_histogram;
  histogram<size_t> count_histogram;
  const size_t array_size = array.m_array.size();

  for (size_t i = 0; i < array_size; ++i) {
    const hsplitter::tile_map_subarray &sub = array.m_array[i];

    if (sub.m_unsorted_count > 0) {
      sub.sort_array();
    }
    std::cout << "array size: " << sub.m_array.size() << " " << sub.m_array.capacity() << std::endl;

    for (std::pair<uint32_t, uint32_t> p : sub.m_array) {
      const uint64_t key = uint64_t(p.first) | (uint64_t(i) << 32);
      if (key == last_key) {
        assert(p.second > last_val);
        ++inner_count;
        uint32_t dval = p.second - last_val;
        values_histogram.insert(dval);

      } else {
        assert(key > last_key);
        uint64_t dkey = key - last_key;
        keys_histogram.insert(dkey);
        count_histogram.insert(inner_count);
        inner_count = 1;
        last_key = key;
        values_histogram.insert(p.second);
      }
      last_val = p.second;
    }
  }

  out << " keys,   " << keys_histogram << std::endl;
  out << " values, " << values_histogram << std::endl;
  out << " counts, " << count_histogram << std::endl;
}
} // anonymous namespace

int main(int argc, char *argv[]) {
  using tileset_t = hsplitter::tile_map_array;

  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " INFILE\n";
    return 1;
  }

  tileset_t node_tiles, way_tiles, extra_node_tiles, rel_tiles, extra_rel_tiles;

  std::cerr << "First pass..." << std::endl;
  { // first pass through file
    osmium::io::File infile(argv[1]);
    osmium::io::Reader reader(infile, osmium::osm_entity_bits::nwr);
    auto itr = osmium::io::make_input_iterator_range<osmium::OSMObject>(reader).begin();

    const auto end = osmium::io::InputIterator<osmium::io::Reader, osmium::OSMObject>();
    std::cerr << " >> nodes" << std::endl;
    node_tiles = hsplitter::tiles_for_nodes<tileset_t, 14>(itr, end);
    node_tiles.freeze();

    std::cerr << " >> ways" << std::endl;
    auto pair_way = hsplitter::tiles_for_ways(itr, end, node_tiles);
    way_tiles = std::move(pair_way.first);
    extra_node_tiles = std::move(pair_way.second);
    way_tiles.freeze();
    extra_node_tiles.freeze();

    std::cerr << " >> relations" << std::endl;
    auto pair_rel = hsplitter::tiles_for_relations(itr, end, node_tiles,
                                                   way_tiles, extra_node_tiles);
    rel_tiles = std::move(pair_rel.first);
    extra_rel_tiles = std::move(pair_rel.second);
    rel_tiles.freeze();
    extra_rel_tiles.freeze();

    reader.close();
  }

  std::cout << "node_tiles:\n"; print_histograms(node_tiles, std::cout);
  std::cout << "way_tiles:\n"; print_histograms(way_tiles, std::cout);
  std::cout << "extra_node_tiles:\n"; print_histograms(extra_node_tiles, std::cout);
  std::cout << "rel_tiles:\n"; print_histograms(rel_tiles, std::cout);
  std::cout << "extra_rel_tiles:\n"; print_histograms(extra_rel_tiles, std::cout);

  return 1;

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
