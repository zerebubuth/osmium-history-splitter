# CMake config for osmium-history-splitter
#
# by Matt Amos <zerebubuth@gmail.com>
# based on osmium-tool's CMakeLists.txt by Jochen Topf.

file(GLOB OSMIUM_HISTORY_SPLITTER_SOURCE_FILES *.cpp)
add_executable(osmium-history-splitter ${OSMIUM_HISTORY_SPLITTER_SOURCE_FILES})
target_link_libraries(osmium-history-splitter ${Boost_LIBRARIES} ${OSMIUM_LIBRARIES})
install(TARGETS osmium-history-splitter DESTINATION bin)
