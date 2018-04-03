# Header files need not be enumerated
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_BINARY_DIR})

###########
# LIB REDSET #
###########

# Install header files
LIST(APPEND libredset_install_headers
	redset.h
)
INSTALL(FILES ${libredset_install_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

LIST(APPEND libredset_srcs
    mpi_rank_str.c
    redset_util.c
    redset_io.c
    redset.c
    redset_single.c
    redset_partner.c
    redset_xor.c
)

# KVTREE Library
ADD_LIBRARY(redset_o OBJECT ${libredset_srcs})

ADD_LIBRARY(redset SHARED $<TARGET_OBJECTS:redset_o>)
TARGET_LINK_LIBRARIES(redset ${KVTREE_EXTERNAL_LIBS})
SET_TARGET_PROPERTIES(redset PROPERTIES OUTPUT_NAME redset CLEAN_DIRECT_OUTPUT 1)
INSTALL(TARGETS redset DESTINATION ${CMAKE_INSTALL_LIBDIR})

ADD_LIBRARY(redset-static STATIC $<TARGET_OBJECTS:redset_o>)
IF(KVTREE_LINK_STATIC)
        SET_TARGET_PROPERTIES(redset-static PROPERTIES LINK_SEARCH_START_STATIC 1)
        SET_TARGET_PROPERTIES(redset-static PROPERTIES LINK_SEARCH_END_STATIC 1)
ENDIF(KVTREE_LINK_STATIC)
TARGET_LINK_LIBRARIES(redset-static ${KVTREE_EXTERNAL_LIBS})

SET_TARGET_PROPERTIES(redset-static PROPERTIES OUTPUT_NAME redset CLEAN_DIRECT_OUTPUT 1)
INSTALL(TARGETS redset-static DESTINATION ${CMAKE_INSTALL_LIBDIR})
