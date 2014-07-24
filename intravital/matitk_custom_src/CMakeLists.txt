PROJECT(MATITK)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#-----------------------------------------------------------------------
# ADD CODE TO FIND CMAKE PACKAGES (Libraries) USED BY THE PROJECT HERE
#------------------------------------------------------------------------

	#---------------------
	#  Find Matlab
	#--------------------- 	
	#FIND_PACKAGE(Matlab)
	
	#IF(MATLAB_FOUND)   
	#ELSE(MATLAB_FOUND)   
	#	MESSAGE(FATAL_ERROR "Please Set MATLAB LIBRARIES.") 
	#ENDIF(MATLAB_FOUND)  
	
	SET(MATLAB_FOUND 0)

	#MESSAGE("CMake Generator chosen - ${CMAKE_GENERATOR}" )

	IF(WIN32)

	  #MESSAGE( "This is a windows system ..." )
	  
	  IF(${CMAKE_GENERATOR} MATCHES "Visual Studio 6")
	  
		SET(MATLAB_ROOT "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.0;MATLABROOT]/extern/lib/win32/microsoft")
		SET(MATLAB_DIR_ROOT "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.0;MATLABROOT]")
		
	  ELSE(${CMAKE_GENERATOR} MATCHES "Visual Studio 6")
	  
		IF(${CMAKE_GENERATOR} MATCHES "Visual Studio 7")
		  # Assume people are generally using 7.1,
		  # if using 7.0 need to link to: ../extern/lib/win32/microsoft/msvc70
		  SET(MATLAB_ROOT "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.0;MATLABROOT]/extern/lib/win32/microsoft")
		  SET(MATLAB_DIR_ROOT "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.0;MATLABROOT]")
		  
		ELSE(${CMAKE_GENERATOR} MATCHES "Visual Studio 7")
		
		  IF(${CMAKE_GENERATOR} MATCHES "Visual Studio 8 2005")

			  IF(${CMAKE_GENERATOR} MATCHES "Visual Studio 8 2005 Win64")
			   
				SET(MATLAB_ROOT "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.5;MATLABROOT]/extern/lib/win64/microsoft")	
				SET(MATLAB_DIR_ROOT "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.12;MATLABROOT]")			
			   
			  ELSE(${CMAKE_GENERATOR} MATCHES "Visual Studio 8 2005 Win64")
		  
				SET(MATLAB_ROOT "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.5;MATLABROOT]/extern/lib/win32/microsoft")
				SET(MATLAB_DIR_ROOT "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.5;MATLABROOT]")			
				
			  ENDIF(${CMAKE_GENERATOR} MATCHES "Visual Studio 8 2005 Win64")		  	
			
		  ELSE(${CMAKE_GENERATOR} MATCHES "Visual Studio 8 2005")	  

			IF(${CMAKE_GENERATOR} MATCHES "Visual Studio 9 2008")
			
				IF(${CMAKE_GENERATOR} MATCHES "Visual Studio 9 2008 Win64")
			   
					SET(MATLAB_ROOT "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.12;MATLABROOT]/extern/lib/win64/microsoft")	
					SET(MATLAB_DIR_ROOT "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.12;MATLABROOT]")			
			   
				ELSE(${CMAKE_GENERATOR} MATCHES "Visual Studio 9 2008 Win64")
		  
					SET(MATLAB_ROOT "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.12;MATLABROOT]/extern/lib/win32/microsoft")
					SET(MATLAB_DIR_ROOT "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.12;MATLABROOT]")			
				
				ENDIF(${CMAKE_GENERATOR} MATCHES "Visual Studio 9 2008 Win64")		  	
			
			ELSE(${CMAKE_GENERATOR} MATCHES "Visual Studio 9 2008")
			
				  IF(${CMAKE_GENERATOR} MATCHES "Borland")
				  
					# Same here, there are also: bcc50 and bcc51 directories
					SET(MATLAB_ROOT "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.0;MATLABROOT]/extern/lib/win32/microsoft/bcc54")
					SET(MATLAB_DIR_ROOT "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.0;MATLABROOT]")			
					
				  ELSE(${CMAKE_GENERATOR} MATCHES "Borland")
				  
					IF(MATLAB_FIND_REQUIRED)
					  MESSAGE(FATAL_ERROR "Generator not compatible: ${CMAKE_GENERATOR}")
					ENDIF(MATLAB_FIND_REQUIRED)
					
				  ENDIF(${CMAKE_GENERATOR} MATCHES "Borland")	  	

			ENDIF(${CMAKE_GENERATOR} MATCHES "Visual Studio 9 2008")
		  
		  ENDIF(${CMAKE_GENERATOR} MATCHES "Visual Studio 8 2005")
		  
		ENDIF(${CMAKE_GENERATOR} MATCHES "Visual Studio 7")
		
	  ENDIF(${CMAKE_GENERATOR} MATCHES "Visual Studio 6")
	  
	  MESSAGE( "Matlab Library Root - ${MATLAB_ROOT}" )
	  MESSAGE( "Matlab Dir Root - ${MATLAB_DIR_ROOT}" )
	  
	  FIND_LIBRARY(MATLAB_MEX_LIBRARY
		libmex
		${MATLAB_ROOT}
		)
	  FIND_LIBRARY(MATLAB_MX_LIBRARY
		libmx
		${MATLAB_ROOT}
		)
	  FIND_LIBRARY(MATLAB_ENG_LIBRARY
		libeng
		${MATLAB_ROOT}
		)

		
	  FIND_PATH(MATLAB_INCLUDE_DIR
		"mex.h"
		${MATLAB_DIR_ROOT}
		"${MATLAB_DIR_ROOT}/extern/include"
		)
		
	  FIND_PROGRAM(MATLAB_MEX_EXECUTABLE 
		  mex.bat
		  PATHS ${MATLAB_DIR_ROOT} ${MATLAB_DIR_ROOT}/bin
		 ) 
		 
	ELSE( WIN32 )
	
	  IF(CMAKE_SIZEOF_VOID_P EQUAL 4)
		# Regular x86
		SET(MATLAB_ROOT
		  /usr/local/matlab-7sp1/bin/glnx86/
		  /opt/matlab-7sp1/bin/glnx86/
		  $ENV{HOME}/matlab-7sp1/bin/glnx86/
		  $ENV{HOME}/redhat-matlab/bin/glnx86/
		  )
	  ELSE(CMAKE_SIZEOF_VOID_P EQUAL 4)
		# AMD64:
		SET(MATLAB_ROOT
		  /usr/local/matlab-7sp1/bin/glnxa64/
		  /opt/matlab-7sp1/bin/glnxa64/
		  $ENV{HOME}/matlab7_64/bin/glnxa64/
		  $ENV{HOME}/matlab-7sp1/bin/glnxa64/
		  $ENV{HOME}/redhat-matlab/bin/glnxa64/
		  )
	  ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 4)
	  
	  FIND_LIBRARY(MATLAB_MEX_LIBRARY
		mex
		${MATLAB_ROOT}
		)
	  FIND_LIBRARY(MATLAB_MX_LIBRARY
		mx
		${MATLAB_ROOT}
		)
	  FIND_LIBRARY(MATLAB_ENG_LIBRARY
		eng
		${MATLAB_ROOT}
		)
	  FIND_PATH(MATLAB_INCLUDE_DIR
		"mex.h"
		"/usr/local/matlab-7sp1/extern/include/"
		"/opt/matlab-7sp1/extern/include/"
		"$ENV{HOME}/matlab-7sp1/extern/include/"
		"$ENV{HOME}/redhat-matlab/extern/include/"
		)

	  FIND_PROGRAM(MATLAB_MEX_EXECUTABLE 
		  mex.bat
		  PATHS ${MATLAB_DIR_ROOT} ${MATLAB_DIR_ROOT}/bin
		 ) 
		
	ENDIF(WIN32)

	#MESSAGE( "Matlab Mex exe - ${MATLAB_MEX_EXECUTABLE}" )

	# This is common to UNIX and Win32:
	SET(MATLAB_LIBRARIES
	  ${MATLAB_MEX_LIBRARY}
	  ${MATLAB_MX_LIBRARY}
	  ${MATLAB_ENG_LIBRARY}
	)

	IF(MATLAB_INCLUDE_DIR AND MATLAB_LIBRARIES)
	  SET(MATLAB_FOUND 1)
	ENDIF(MATLAB_INCLUDE_DIR AND MATLAB_LIBRARIES)

	MARK_AS_ADVANCED(
	  MATLAB_LIBRARIES
	  MATLAB_MEX_LIBRARY
	  MATLAB_MX_LIBRARY
	  MATLAB_ENG_LIBRARY
	  MATLAB_INCLUDE_DIR
	  MATLAB_FOUND
	  MATLAB_ROOT
	  MATLAB_DIR_ROOT
	  MATLAB_MEX_EXECUTABLE
	)
	
	#---------------------
	#  Find ITK 
	#--------------------- 

	FIND_PACKAGE(ITK) 

	IF(ITK_FOUND)   
		INCLUDE(${ITK_USE_FILE}) 
	ELSE(ITK_FOUND)   
		MESSAGE(FATAL_ERROR "Please Set ITK_DIR.") 
	ENDIF(ITK_FOUND)  

	#MESSAGE(  "ITK Library Dir List --- ${ITK_LIBRARY_DIRS}" )	
	MESSAGE("ITK Version --- ${ITK_VERSION}")
	
	#---------------------
	#  ITK Libraries Needed
	#--------------------- 		
	SET(ITK_LIBRARIES_NEEDED

	${ITK_LIBRARIES}
	)

	STRING(REPLACE ";"
	 ";-I" NEW_ITK_INCLUDE_DIRS
	 "${ITK_INCLUDE_DIRS}"
	)	
	
	STRING(REPLACE ";"
	 ";-l" ITK_LIBRARIES_NEEDED_MEX
	 "${ITK_LIBRARIES_NEEDED}"
	)	
	
#---------------------------------------- 
#  INCLUDE DIRECTORIES 
#---------------------------------------- 

INCLUDE_DIRECTORIES(
${MATITK_BINARY_DIR} 
${MATITK_SOURCE_DIR}  

${MATLAB_INCLUDE_DIR}
${ITK_INCLUDE_DIRS}
)

#------------------------------- 
# Project Header(*.h) Files 
#-------------------------------  
SET(MATITK_SRCS_H   
	
MATITKTemplatedVariables.h

ParameterContainer.h
ParameterContainer.inl

itkcore.h
itkcore.inl

itkfiltercore.h
itkregistrationcore.h
itksegmentationcore.h

matitk_custom.h

seedcontainer.h
seedcontainer.inl

typedefs.inl
functionDescr.inl

)

#------------------------------- 
# Project Source(*.CPP) Files 
#-------------------------------  
SET(MATITK_SRCS_CXX   	

itkfiltercore.inl
itkregistrationcore.inl
itksegmentationcore.inl

matitk_custom.cxx

)  

ADD_LIBRARY(matitk MODULE ${MATITK_SRCS_H} ${MATITK_SRCS_CXX} )

ADD_DEFINITIONS(-DMATLAB_MEX_FILE)

SET_TARGET_PROPERTIES(matitk
PROPERTIES
LINK_FLAGS "/export:mexFunction"
)

TARGET_LINK_LIBRARIES(matitk ${MATLAB_LIBRARIES} ${ITK_LIBRARIES})

#MESSAGE( "ITK Library Default List -- ${ITK_LIBRARIES}" )
#MESSAGE( "MATLAB Library List -- ${MATLAB_LIBRARIES}" )	

#---------------------------------------- 
#  Mex File Generation
#---------------------------------------- 

IF (UNIX)
 
	ADD_CUSTOM_TARGET(matitk_custom ALL ${MATLAB_MEX_EXECUTABLE} -g ${MATITK_SOURCE_DIR}/matitk_custom.cxx -DUNIX_BUILD ${MATLAB_VERSION_SWITCH} -I${NEW_ITK_INCLUDE_DIRS} -I${ITK_SOURCE_DIR}/Code/Patented ${ITK_IO} ${ITK_ALGORITHMS} ${ITK_STATISTICS} ${ITK_SPATIALOBJECT} ${ITK_BASIC_FILTERS} ${ITK_COMMON} ${ITK_VNL_INST} ${ITK_VNL_ALGO} ${ITK_NRRDIO} ${ITK_GDCM} ${ITK_PNG} ${ITK_TIFF} ${ITK_NUMERICS} ${ITK_METAIO} ${ITK_NIFTIIO} ${ITK_ZNZ} ${ITK_SYS} ${ITK_ZLIB} ${ITK_EXPAT} ${ITK_DICOMPARSER} ${ITK_JPEG8} ${ITK_JPEG12} ${ITK_JPEG16} ${ITK_VNL} ${ITK_VCL} ${ITK_NETLIB} -lpthread -ldl -lm
       )
ELSE(UNIX)

	#MESSAGE( "MEX ITK Include Dirs -- ${NEW_ITK_INCLUDE_DIRS}" )		
	#MESSAGE( "MEX WINDOWS LIBS -- ${MS_WINDOWS_LIB_LIST}" )   
	  
	#SET(MEX_COMPILE_COMMAND "${MATLAB_MEX_EXECUTABLE} -largeArrayDims -g -c -O ${MATITK_SOURCE_DIR}/matitk.cxx -I${NEW_ITK_INCLUDE_DIRS} -I${ITK_SOURCE_DIR}/Code/Patented")
	#SET(MEX_LINK_COMMAND "${MATLAB_MEX_EXECUTABLE} -largeArrayDims -g matitk.obj -O -L${ITK_DIR}/bin/debug -l${ITK_LIBRARIES_NEEDED_MEX} ${MS_WINDOWS_LIB_LIST}")
	#ADD_CUSTOM_TARGET(matitklib ALL ${MEX_COMPILE_COMMAND}) 
	#ADD_CUSTOM_TARGET(matitklink ALL ${MEX_LINK_COMMAND} DEPENDS matitklib)
	#MESSAGE( "MEX COMPILE COMMAND ISSUED -- ${MEX_COMPILE_COMMAND}" )  
	#MESSAGE( "MEX LINK COMMAND ISSUED -- ${MEX_LINK_COMMAND}" ) 

	IF(${ITK_VERSION} MATCHES "4.-*")
		
		#MESSAGE("ITK Version is 4 and above")
		STRING(REPLACE ";"
			 "-4.1;-l" ITK_LIBRARIES_NEEDED_MEX
			 "${ITK_LIBRARIES_NEEDED}"
			)								
		SET(ITK_LIBRARIES_NEEDED_MEX ${ITK_LIBRARIES_NEEDED_MEX}-4.1)
		MESSAGE( "MEX ITK LIBS -- ${ITK_LIBRARY_DIRS}" ) 
		SET(MEX_COMPILE_COMMAND "${MATLAB_MEX_EXECUTABLE} -largeArrayDims -g -v -I${NEW_ITK_INCLUDE_DIRS} -L${ITK_DIR}/lib/debug -l${ITK_LIBRARIES_NEEDED_MEX} ${MATITK_SOURCE_DIR}/matitk_custom.cxx")
		MESSAGE( "MEX COMPILE COMMAND ISSUED -- ${MEX_COMPILE_COMMAND}" )  		
		ADD_CUSTOM_TARGET(matitkmexgen ALL ${MATLAB_MEX_EXECUTABLE} -largeArrayDims -g -v -I${NEW_ITK_INCLUDE_DIRS} -L${ITK_DIR}/lib/debug -l${ITK_LIBRARIES_NEEDED_MEX} ${MATITK_SOURCE_DIR}/matitk_custom.cxx)
	
	ELSE(${ITK_VERSION} MATCHES "4.-*")
	
		#MESSAGE("ITK Version is not 4 and above")
		STRING(REPLACE ";"
			 ";-l" ITK_LIBRARIES_NEEDED_MEX
			 "${ITK_LIBRARIES_NEEDED}"
			)			
		SET(MEX_COMPILE_COMMAND "${MATLAB_MEX_EXECUTABLE} -largeArrayDims -g -v -I${NEW_ITK_INCLUDE_DIRS} -L${ITK_DIR}/lib/debug -l${ITK_LIBRARIES_NEEDED_MEX} ${MATITK_SOURCE_DIR}/matitk_custom.cxx")
		MESSAGE( "MEX COMPILE COMMAND ISSUED -- ${MEX_COMPILE_COMMAND}" )  
		ADD_CUSTOM_TARGET(matitkmexgen ALL ${MATLAB_MEX_EXECUTABLE} -largeArrayDims -g -v -I${NEW_ITK_INCLUDE_DIRS} -L${ITK_DIR}/lib/debug -l${ITK_LIBRARIES_NEEDED_MEX} ${MATITK_SOURCE_DIR}/matitk_custom.cxx)	
		
	ENDIF(${ITK_VERSION} MATCHES "4.-*")
	
	
ENDIF(UNIX)

SET_DIRECTORY_PROPERTIES(PROPERTIES
 ADDITIONAL_MAKE_CLEAN_FILES "*.mex*"
)