#include <osmium/io/any_input.hpp>
#include <osmium/io/any_compression.hpp>
#include <osmium/io/input_iterator.hpp>

#include "splitter.hpp"

#include <iostream>
#include <map>
#include <set>


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

// struct tile_map_subarray {
//   ??? &operator[](uint32_t k) {
//   }
// };

// struct tile_map_array {
//   typedef key_type int64_t;
//   typedef mapped_type ???;

//   mapped_type &operator[](key_type k) {
//     ASSERT(k >= 0);
//     uint32_t idx = k >> 32;
//     uint32_t prt = k & ((uint64_t(1) << 32) - 1);
//     if (idx >= m_array.size()) {
//       m_array.resize(idx + 1);
//     }
//     tile_map_subarray &sub = m_array[idx];
//     return sub[prt];
//   }

//   std::vector<tile_map_subarray> m_array;
// };

} // anonymous namespace

int main(int argc, char *argv[]) {
  using tileset_t = tileset;

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
