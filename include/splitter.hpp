#ifndef OSMIUM_HISTORY_SPLITTER_SPLITTER_HPP
#define OSMIUM_HISTORY_SPLITTER_SPLITTER_HPP

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

#include <osmium/geom/coordinates.hpp>
#include <osmium/geom/mercator_projection.hpp>
#include <osmium/osm/node.hpp>

namespace hsplitter {

using tile_t = uint32_t;

inline uint32_t interleave(uint32_t n) {
  n &= 0x0000ffffUL;
  n = (n | (n <<  8)) & 0x00ff00ffUL;
  n = (n | (n <<  4)) & 0x0f0f0f0fUL;
  n = (n | (n <<  2)) & 0x33333333UL;
  n = (n | (n <<  1)) & 0x55555555UL;
  return n;
}

inline tile_t morton_code(uint16_t x, uint16_t y) {
  return (interleave(x) << 1) | interleave(y);
}

template <typename Item>
struct osmium_item_type {};

template <>
struct osmium_item_type<osmium::Node> {
  static const osmium::item_type type = osmium::item_type::node;
};

template <typename Item, typename Iterator, typename Func>
inline void iterate(Iterator &it, const Iterator &end, Func func) {
  while ((it != end) && (it->type() == osmium_item_type<Item>::type)) {
    const Item &item = static_cast<const Item &>(*it);
    func(item);
    ++it;
  }
}

template <typename TileSet, typename Iterator>
inline TileSet tiles_for_nodes(Iterator &it, const Iterator &end) {
  static const double max_coord = osmium::geom::detail::max_coordinate_epsg3857;
  TileSet tiles;

  iterate<osmium::Node>(it, end, [&tiles] (const osmium::Node &node) {
      if (node.visible() && node.location().valid()) {
        // note: we've already checked for validity, no need to do it
        // again.
        double lat = node.location().lat_without_check();

        // nodes outside this range can't be projected to web mercator
        // tiles.
        if ((lat >= -osmium::geom::MERCATOR_MAX_LAT) &&
            (lat <= osmium::geom::MERCATOR_MAX_LAT)) {
          // project to mercator space - we'll use that for splitting
          // into tiles so that it matches the z/x/y (for z=16) that
          // people expect.
          osmium::geom::Coordinates coords(node.location().lon_without_check(), lat);
          coords = osmium::geom::lonlat_to_mercator(coords);

          // tiles are named by the morton code of the x,y coord pair.
          // clamp tile coords to 0, 65535 to avoid overflow.
          double x = double(1 << 16) * 0.5 * (coords.x + max_coord) / max_coord;
          if (x <     0.0) { x =     0.0; }
          if (x > 65535.0) { x = 65535.0; }
          double y = double(1 << 16) * 0.5 * (max_coord - coords.y) / max_coord;
          if (y <     0.0) { y =     0.0; }
          if (y > 65535.0) { y = 65535.0; }
          
          auto tile_id = morton_code(uint16_t(x), uint16_t(y));
          tiles[node.id()].insert(tile_id);
        }
      }
    });

  return std::move(tiles);
}

} // namespace hsplitter

#endif /* OSMIUM_HISTORY_SPLITTER_SPLITTER_HPP */
