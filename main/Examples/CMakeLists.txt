set(CompactStar_main_Example_SRC_Files

    sig_omega_rho_nstar
    Fermi_gas_many
)

foreach(src ${CompactStar_main_Example_SRC_Files})
    add_executable(${src} "${src}.cpp")  
    target_link_libraries(${src} PUBLIC CompactStar)
    target_include_directories(${src} PUBLIC  
                                ${CMAKE_CURRENT_SOURCE_DIR}/include
                                ${GSL_INCLUDE_DIRS}
                                ${DEP_INC_DIR}
                                )

    install(FILES "${src}.cpp" DESTINATION main/Example)
    install(TARGETS ${src} DESTINATION main/Example/bin)
endforeach()