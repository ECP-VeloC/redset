###############
# Build tests
###############
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR})

################
# Add tests to ctest
################

IF(REDSET_LINK_STATIC)
   SET(redset_lib redset::redset-static)
ELSE()
   SET(redset_lib redset::redset)
ENDIF()

ADD_EXECUTABLE(redset_test test_redset.c)
TARGET_LINK_LIBRARIES(redset_test ${redset_lib})
REDSET_ADD_TEST_PARALLEL(redset_test 256 "")

ADD_EXECUTABLE(test_config test_config.c)
TARGET_LINK_LIBRARIES(test_config ${redset_lib})
REDSET_ADD_TEST_SERIAL(test_config 256 "")

####################
# make a verbose "test" target named "check"
####################

ADD_CUSTOM_TARGET(check COMMAND ${CMAKE_CTEST_COMMAND} --verbose)
