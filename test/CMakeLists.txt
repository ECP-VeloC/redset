###############
# Build tests
###############

LIST(APPEND er_test_srcs
    test_redset.c
    test_er.c
)

ADD_EXECUTABLE(test_er ${er_test_srcs})
TARGET_LINK_LIBRARIES(test_er er kvtree)

################
# Add tests to ctest
################

# rountdrip tests
ADD_TEST(roundtrip_tests test_er)

# print test script
FIND_PROGRAM(PRINT_TEST ${CMAKE_CURRENT_SOURCE_DIR}/print_test)
# ADD_TEST(NAME print_test COMMAND ${PRINT_TEST} ${PROJECT_BINARY_DIR}/src/kvtree_print_file ${CMAKE_CURRENT_SOURCE_DIR}/files)

####################
# make a verbose "test" target named "check"
####################

ADD_CUSTOM_TARGET(check COMMAND ${CMAKE_CTEST_COMMAND} --verbose)
