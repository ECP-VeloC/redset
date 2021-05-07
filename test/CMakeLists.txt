###############
# Build tests
###############
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR})

################
# Add tests to ctest
################

IF(REDSET_LINK_STATIC)
   SET(redset_lib redset-static)
ELSE()
   SET(redset_lib redset)
ENDIF()

ADD_EXECUTABLE(redset_test test_redset.c)
TARGET_LINK_LIBRARIES(redset_test ${SPATH_EXTERNAL_LIBS} ${redset_lib})
ADD_TEST(NAME redset_test COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 -N 2 ./redset_test)

ADD_EXECUTABLE(test_config test_config.c)
TARGET_LINK_LIBRARIES(test_config ${SPATH_EXTERNAL_LIBS} ${redset_lib})
ADD_TEST(NAME test_config COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ./test_config)

####################
# make a verbose "test" target named "check"
####################

ADD_CUSTOM_TARGET(check COMMAND ${CMAKE_CTEST_COMMAND} --verbose)
