find_package(Git)
if(Git_FOUND)
    execute_process(
        COMMAND
            sh -c
            "${GIT_EXECUTABLE} describe --abbrev=4 --dirty=-dirty --always --tags  | cut -c 2- | tr -d '\n' | sed --regexp-extended s/-/./"
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE PROJECT_VERSION_GIT)
endif()

# get project root directory
get_filename_component(CMAKE_PARENT_LIST_DIR ${CMAKE_PARENT_LIST_FILE}
                       DIRECTORY)
get_filename_component(CMAKE_PARENT_LIST_DIR ${CMAKE_PARENT_LIST_DIR} DIRECTORY)

# execute_process(
#     COMMAND sh -c "cat version.txt | tr -d '\n'"
#     WORKING_DIRECTORY "${CMAKE_PARENT_LIST_DIR}"
#     OUTPUT_VARIABLE PROJECT_VERSION)