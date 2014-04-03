

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -std=c++11")
endif()

set(LIBMAW_INCLUDE_PATH /home/chris/Documents/PhD/code/libmaw/src/)
set(LIBMAW_LIBRARY /home/chris/Documents/PhD/code/libmaw/lib/libmaw.a)

#sym link /opt/matlab-7sp1 to the real matlab root will make this work properly
find_package(Matlab REQUIRED)
find_package(Armadillo 3.0 REQUIRED)


set( HEADER_FILES 
	${AMRADILLO_INCLUDE_DIR}
	${LIBMAW_INCLUDE_PATH} 
	${MATLAB_INCLUDE_DIR})

set( LIB_FILES 
	${ARMADILLO_LIBRARIES} 
	${MATLAB_ENG_LIBRARY} 
	${MATLAB_MX_LIBRARY}  
	${LIBMAW_LIBRARY} )
 
set( SOURCES 
	main.cpp
	triangle.cpp
	triangulation.cpp
	bounding_rect.cpp
)

include_directories(${HEADER_FILES})
add_executable(umbra ${SOURCES}  ${HEADER_FILES})

target_link_libraries (umbra ${LIB_FILES})