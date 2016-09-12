#
#	$Id: $
#

if (HAVE_M_LIBRARY)
	# link the math library
	set (LINK_LIBS m)
endif (HAVE_M_LIBRARY)

include_directories (include ../../gmtsar ../ALOS_preproc/include ${GMT_INCLUDE_DIR})

set (LINK_LIBS ${LINK_LIBS} gmtsar ALOS)

add_executable (ers_line_fixer ers_line_fixer/ers_line_fixer.c)
target_link_libraries (ers_line_fixer ${LINK_LIBS})

add_executable (read_data_file_ccrs read_data_file_ccrs/read_data_file_ccrs.c)
target_link_libraries (read_data_file_ccrs ${LINK_LIBS})

add_executable (read_data_file_dpaf read_data_file_dpaf/read_data_file_dpaf.c)
target_link_libraries (read_data_file_dpaf ${LINK_LIBS})

add_executable (read_sarleader_ccrs read_sarleader_ccrs/read_sarleader_ccrs.c read_sarleader_ccrs/make_prm_ccrs.c)
target_link_libraries (read_sarleader_ccrs ${LINK_LIBS})

add_executable (read_sarleader_dpaf read_sarleader_dpaf/read_sarleader_dpaf.c read_sarleader_dpaf/make_prm_dpaf.c)
target_link_libraries (read_sarleader_dpaf ${LINK_LIBS})

add_executable (ccrs2dpaf utility/ccrs2dpaf.c)
target_link_libraries (ccrs2dpaf ${LINK_LIBS})


# add the install targets
install (TARGETS ers_line_fixer read_data_file_ccrs read_data_file_dpaf read_sarleader_ccrs read_sarleader_dpaf ccrs2dpaf
	ARCHIVE DESTINATION lib
	COMPONENT Runtime
	LIBRARY DESTINATION lib
	COMPONENT Runtime
	RUNTIME DESTINATION bin
	COMPONENT Runtime)

# vim: textwidth=78 noexpandtab tabstop=2 softtabstop=2 shiftwidth=2
