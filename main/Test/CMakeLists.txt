set(CompactStar_main_Test_SRC_Files
	tov_debug_main
	psr_debug_main
	spin_evol_main
	spin_evol_2_main
	spin_therm_evol_main
	ns_build_main
	spin_therm_evol_2_main
	compose_thermo_main
)

foreach(src ${CompactStar_main_Test_SRC_Files})
	add_executable(${src} "${src}.cpp")
	target_link_libraries(${src} PUBLIC CompactStar)
	target_include_directories(${src} PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}/include
		${GSL_INCLUDE_DIRS}
		${DEP_INC_DIR}
	)

	install(FILES "${src}.cpp" DESTINATION main/Test)
	install(TARGETS ${src} DESTINATION main/Test/bin)
endforeach()