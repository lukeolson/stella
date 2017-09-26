add_executable(atest test.cc)
target_link_libraries(atest stella cedar ${cedar_LIBRARIES} ${cedar-deps})
