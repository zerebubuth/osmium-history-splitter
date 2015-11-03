# CMake config for osmium-history-splitter
#
# by Matt Amos <zerebubuth@gmail.com>
# based on osmium-tool's CMakeLists.txt by Jochen Topf.

file(GLOB OSMIUM_HISTORY_SPLITTER_VX2FMT_SOURCE_FILES *.cpp)
add_library(osmium-history-splitter-vx2fmt ${OSMIUM_HISTORY_SPLITTER_VX2FMT_SOURCE_FILES})
target_link_libraries(osmium-history-splitter-vx2fmt ${Boost_LIBRARIES} ${OSMIUM_LIBRARIES})
