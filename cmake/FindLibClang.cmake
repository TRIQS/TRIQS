MESSAGE(STATUS "-------- LibClang detection -------------")

SET(LIBCLANG_DETECTION_COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libclang)
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
 SET(LIBCLANG_DETECTION_COMMAND ${LIBCLANG_DETECTION_COMMAND}.dylib)
ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
 SET(LIBCLANG_DETECTION_COMMAND ${LIBCLANG_DETECTION_COMMAND}.so)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

EXECUTE_PROCESS(COMMAND ${LIBCLANG_DETECTION_COMMAND}
                OUTPUT_VARIABLE TRIQS_LIBCLANG_LOCATION OUTPUT_STRIP_TRAILING_WHITESPACE)
GET_FILENAME_COMPONENT(TRIQS_LIBCLANG_LOCATION ${TRIQS_LIBCLANG_LOCATION} ABSOLUTE)

MESSAGE(STATUS "LibClang location: ${TRIQS_LIBCLANG_LOCATION}")

SET(LIBCLANG_FLAGS_DETECTION_COMMAND "${CMAKE_CXX_COMPILER}" -E -x c++ -v -)
EXECUTE_PROCESS(COMMAND ${LIBCLANG_FLAGS_DETECTION_COMMAND}
                INPUT_FILE /dev/null
                ERROR_VARIABLE _compiler_output OUTPUT_QUIET)

STRING(REGEX MATCH "#include <...> search starts here:\n(.*)End of search list." _matches "${_compiler_output}")
STRING(REPLACE "\n" ";" CMAKE_MATCH_1 ${CMAKE_MATCH_1})
FOREACH(include_path IN ITEMS ${CMAKE_MATCH_1})
 STRING(STRIP ${include_path} include_path)
 SET(TRIQS_LIBCLANG_CXX_ADDITIONAL_FLAGS "${TRIQS_LIBCLANG_CXX_ADDITIONAL_FLAGS} -I${include_path}")
ENDFOREACH(include_path IN ITEMS ${CMAKE_MATCH_1})

MESSAGE(STATUS "LibClang additional flags: ${TRIQS_LIBCLANG_CXX_ADDITIONAL_FLAGS}")