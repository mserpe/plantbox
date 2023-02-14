#
# 1. Make CPlantBox library
#
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
include_directories(${PROJECT_BINARY_DIR}/src/external/suitsparse/include)
include_directories(${PROJECT_BINARY_DIR}/src/external/sundials/include)
add_library(CPlantBox SHARED                        
            organparameter.cpp
            rootparameter.cpp
			seedparameter.cpp
			leafparameter.cpp
			stemparameter.cpp
            Organ.cpp
            Root.cpp
            RootDelay.cpp            
            Seed.cpp
			Stem.cpp
			Leaf.cpp
			Plant.cpp
            Organism.cpp
            Plant.cpp            
            RootSystem.cpp
            MappedOrganism.cpp
            XylemFlux.cpp
     		sdf.cpp
            SegmentAnalyser.cpp            
            tropism.cpp            
			external/tinyxml2/tinyxml2.cpp
            external/aabbcc/AABB.cc
            external/gauss_legendre/gauss_legendre.cpp
			)

set_target_properties(CPlantBox PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/)

# makign cmake copy the files to the debug folder
set_target_properties(CPlantBox PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/debug_test/Debug)
# checking if the destination is correct by a message
message(STATUS "CMAKE_BINARY_DIR: " ${CMAKE_BINARY_DIR}/debug_test/Debug)
# making sure with a custom command
add_custom_command(TARGET CPlantBox POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy
                       $<TARGET_FILE:CPlantBox> ${CMAKE_BINARY_DIR}/debug_test/Debug)
#adding the respective target
add_custom_target(CPlantBox_copy ALL DEPENDS CPlantBox)

target_link_libraries(CPlantBox
						)

 install(TARGETS CPlantBox CONFIGURATIONS Debug RUNTIME_DESTINATION ${CMAKE_BINARY_DIR}/debug_test/${CMAKE_BUILD_TYPE})

#
# 2. Make CPlantBox Pyhthon binding
#
# find_package(pybind11 REQUIRED) 
# set(PYBIND11_PYTHON_VERSION 3.7)
find_package( PythonInterp 3.7 REQUIRED )
find_package( PythonLibs 3.7 REQUIRED )
add_subdirectory(external/pybind11)

pybind11_add_module(plantbox SHARED 
            sdf.cpp
            organparameter.cpp
            Organ.cpp
            Organism.cpp
            
            rootparameter.cpp
			seedparameter.cpp
			leafparameter.cpp
			stemparameter.cpp			
            
            Root.cpp
            RootDelay.cpp
            Seed.cpp
			Stem.cpp
			Leaf.cpp
			Plant.cpp
			            
            Plant.cpp            
            RootSystem.cpp
            MappedOrganism.cpp
			XylemFlux.cpp
           
            SegmentAnalyser.cpp            
            tropism.cpp
            
			external/tinyxml2/tinyxml2.cpp
            external/aabbcc/AABB.cc
            external/gauss_legendre/gauss_legendre.cpp
	    	
			PyPlantBox.cpp
								
			)
target_link_libraries(plantbox PUBLIC)
			
set_target_properties(plantbox PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/)
