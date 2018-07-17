# CMAKE for root_persistency test

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/../include ${VECGEOM_INCLUDE_DIR})

add_executable(test_unplacedvol "test_unplacedvol.cpp")
install(TARGETS test_unplacedvol DESTINATION ${CMAKE_INSTALL_PREFIX})
target_link_libraries(test_unplacedvol RootPClassesDict ${VECGEOM_LIBRARIES} ${ROOT_LIBRARIES})







