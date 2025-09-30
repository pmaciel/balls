foreach(test
        cities
        )
        add_executable(${test} ${test}.cc)
        target_link_libraries(${test} PUBLIC balls)
        add_test(NAME test_${test} COMMAND ${test})
endforeach()
