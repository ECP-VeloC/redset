###############
# Build tests
###############

################
# Add tests to ctest
################

ADD_EXECUTABLE(redset_test test_redset.c)
TARGET_LINK_LIBRARIES(redset_test ${SPATH_EXTERNAL_LIBS} redset)
REDSET_ADD_TEST(redset_test 256 "")
#ADD_TEST(NAME redset_test COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3 -N3 ./redset_test)

####################
# make a verbose "test" target named "check"
####################

ADD_CUSTOM_TARGET(check COMMAND ${CMAKE_CTEST_COMMAND} --verbose)
