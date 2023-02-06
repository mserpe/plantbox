
# Projectname: ${projectname}
# PROJECTNAME: ${PROJECTNAME_UPPER}
# path: ${librarypath}

get_filename_component(Folder ${CMAKE_CURRENT_LIST_DIR} NAME)
string(REPLACE " " "_" Folder ${Folder})

message("This is ${Folder}")

file(GLOB TESTSOURCES ./*.cpp)
file(GLOB TESTHEADERS ./*.h)


add_executable(${Folder}
  ${TESTSOURCES}
  ${TESTHEADERS}
)

target_include_directories(${Folder}
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
  PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_LIST_DIR}/../src/

)

target_link_libraries(${Folder} PRIVATE CPlantBox)

