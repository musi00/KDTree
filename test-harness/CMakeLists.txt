add_executable(test-harness test-harness.cpp)
#target_link_libraries(example_pcl_don don ${Boost_LIBRARIES} ${PCL_LIBRARIES})
# enable C++11 option for this target
set_property(TARGET test-harness PROPERTY CXX_STANDARD 11)
set_property(TARGET test-harness PROPERTY CXX_STANDARD_REQUIRED ON)
