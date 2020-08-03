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

LIST(APPEND cliredset_noMPI_srcs
    redset_err_serial.c
    redset_util.c
    redset_io.c
    redset_lofi.c
    redset_partner_serial.c
    redset_xor_serial.c
)

LIST(APPEND libredset_srcs
    redset_err_mpi.c
    redset_util.c
    redset_util_mpi.c
    redset_io.c
    redset_lofi.c
    redset.c
    redset_single.c
    redset_partner.c
    redset_xor.c
    redset_reedsolomon.c
)

# REDSET Library
ADD_LIBRARY(redset_o OBJECT ${libredset_srcs})

ADD_LIBRARY(redset SHARED $<TARGET_OBJECTS:redset_o>)
TARGET_LINK_LIBRARIES(redset ${REDSET_EXTERNAL_LIBS})
SET_TARGET_PROPERTIES(redset PROPERTIES OUTPUT_NAME redset CLEAN_DIRECT_OUTPUT 1)
INSTALL(TARGETS redset DESTINATION ${CMAKE_INSTALL_LIBDIR})

ADD_LIBRARY(redset-static STATIC $<TARGET_OBJECTS:redset_o>)
IF(REDSET_LINK_STATIC)
        SET_TARGET_PROPERTIES(redset-static PROPERTIES LINK_SEARCH_START_STATIC 1)
        SET_TARGET_PROPERTIES(redset-static PROPERTIES LINK_SEARCH_END_STATIC 1)
ENDIF(REDSET_LINK_STATIC)
TARGET_LINK_LIBRARIES(redset-static ${REDSET_EXTERNAL_LIBS})

SET_TARGET_PROPERTIES(redset-static PROPERTIES OUTPUT_NAME redset CLEAN_DIRECT_OUTPUT 1)
INSTALL(TARGETS redset-static DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Non-MPI library for CLI
ADD_LIBRARY(redset_base SHARED ${cliredset_noMPI_srcs})
TARGET_LINK_LIBRARIES(redset_base ${REDSET_EXTERNAL_LIBS})
INSTALL(TARGETS redset_base DESTINATION ${CMAKE_INSTALL_LIBDIR})

ADD_LIBRARY(redset_base-static STATIC ${cliredset_noMPI_srcs})
TARGET_LINK_LIBRARIES(redset_base-static ${REDSET_EXTERNAL_LIBS})
SET_TARGET_PROPERTIES(redset_base-static PROPERTIES OUTPUT_NAME redset_base CLEAN_DIRECT_OUTPUT 1)
INSTALL(TARGETS redset_base-static DESTINATION ${CMAKE_INSTALL_LIBDIR})

##############
# CLI REDSET #
##############

## Individual binaries generated from a single .c file
#LIST(APPEND cliredset_c_bins
#	redset_rebuild_xor
#)
#
## Build and install C binaries
#FOREACH(bin IN ITEMS ${cliredset_c_bins})
#	ADD_EXECUTABLE(${bin} ${bin}.c)
#	TARGET_LINK_LIBRARIES(${bin} redset_base)
#	INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${bin} DESTINATION ${CMAKE_INSTALL_BINDIR})
#ENDFOREACH(bin IN ITEMS ${cliredset_c_bins})
