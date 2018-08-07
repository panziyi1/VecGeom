# CMAKE for root_persistency test

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/../include ${VECGEOM_INCLUDE_DIR})

add_executable(RootPersistencyTest "RootPersistencyTest.cpp" "unplaced_test.cpp" "logical_test.cpp" "vector_test.cpp" "proxy.cpp")
install(TARGETS RootPersistencyTest DESTINATION ${CMAKE_INSTALL_PREFIX})
target_link_libraries(RootPersistencyTest RootPClassesDict ${VECGEOM_LIBRARIES} ${ROOT_LIBRARIES})