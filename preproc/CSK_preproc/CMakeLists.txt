#
#	$Id: $
#

if (HAVE_M_LIBRARY)
	# link the math library
	set (LINK_LIBS m)
endif (HAVE_M_LIBRARY)

include_directories (include ../../gmtsar ../S1A_preproc/include ${HDF5_INCLUDE_DIR} ${GMT_INCLUDE_DIR})

set (LINK_LIBS ${LINK_LIBS} ${HDF5_LIBRARY} gmtsar xmlC)

if (WIN32)
add_definitions(/DH5_BUILT_AS_DYNAMIC_LIB)
endif (WIN32)

add_executable (make_raw_csk src_raw/make_raw_csk.c)
target_link_libraries (make_raw_csk ${LINK_LIBS})

add_executable (make_slc_csk src_slc/make_slc_csk.c)
target_link_libraries (make_slc_csk ${LINK_LIBS})

# add the install targets
install (TARGETS make_raw_csk make_slc_csk
	ARCHIVE DESTINATION lib
	COMPONENT Runtime
	LIBRARY DESTINATION lib
	COMPONENT Runtime
	RUNTIME DESTINATION bin
	COMPONENT Runtime)

# vim: textwidth=78 noexpandtab tabstop=2 softtabstop=2 shiftwidth=2
