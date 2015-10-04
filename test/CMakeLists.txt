# CMake config for osmium-history-splitter
#
# by Matt Amos <zerebubuth@gmail.com>
# based on osmium-tool's CMakeLists.txt by Jochen Topf.

include_directories(include)
include_directories(../src)
include_directories(../include)


file(GLOB UNIT_TEST_SOURCES *.cpp)
add_executable(unit_tests ${UNIT_TEST_SOURCES})
target_link_libraries(unit_tests ${Boost_LIBRARIES} ${OSMIUM_LIBRARIES})
add_test(NAME unit_tests COMMAND unit_tests)
set_tests_properties(unit_tests PROPERTIES
  LABELS "unit"
  ENVIRONMENT "OSMIUM_HISTORY_SPLITTER_TEST_DATA_DIR=${CMAKE_CURRENT_SOURCE_DIR}")

