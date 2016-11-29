cmake_minimum_required(VERSION 3.3)

# Use the default gcc/g++ version
find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

# require at least gcc 5.0
if (CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0)
    message(WARNING "Require at least gcc-5.0")
endif()


project(kraken)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++14 -O3 -w")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/bin)

# include directory in find path where all dependency modules exist
include_directories(./)

# Common Code sets
set(SOURCE_FILES_BASIC ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc) 
set(SOURCE_FILES_FFT ryggrad/extern/RealFFT/DynArray.hpp ryggrad/extern/RealFFT/FFTReal.hpp ryggrad/extern/RealFFT/OscSinCos.hpp) 
set(SOURCE_FILES_ANNOTQ ryggrad/src/general/AlignmentBlock.cc ryggrad/src/general/Coordinate.cc src/annotationQuery/AnnotationQuery.cc) 
set(SOURCE_FILES_COLA cola/src/cola/AlignmentCola.cc cola/src/cola/Cola.cc cola/src/cola/EditGraph.cc cola/src/cola/NSaligner.cc cola/src/cola/NSGAaligner.cc cola/src/cola/SWGAaligner.cc ryggrad/src/general/Alignment.cc)  
set(SOURCE_FILES_KRAKEN ryggrad/src/general/CodonTranslate.cc ryggrad/src/general/CrossCorr.cc src/kraken/KrakenConfig.cc src/kraken/KrakenMap.cc) 


# AnnotationQuery binaries
set(SOURCE_FILES_COMPAREANNOTWITHREF        ${SOURCE_FILES_ANNOTQ} ${SOURCE_FILES_BASIC} src/annotationQuery/CompareAnnotWithRef.cc) 
set(SOURCE_FILES_GETIOSINGLEEXONTRANSCRIPTS ${SOURCE_FILES_ANNOTQ} ${SOURCE_FILES_BASIC} src/annotationQuery/GetIOSingleExonTranscripts.cc) 
set(SOURCE_FILES_MAF2SATSUMA                ${SOURCE_FILES_ANNOTQ} ${SOURCE_FILES_BASIC} src/annotationQuery/MultiAlignParser.cc src/annotationQuery/MAF2Satsuma.cc)
set(SOURCE_FILES_RUNANNOTATIONQUERY         ${SOURCE_FILES_ANNOTQ} ${SOURCE_FILES_BASIC} src/annotationQuery/RunAnnotationQuery.cc)
set(SOURCE_FILES_RUNGTFCOMPARE              ${SOURCE_FILES_ANNOTQ} ${SOURCE_FILES_BASIC} src/annotationQuery/RunGTFCompare.cc)
set(SOURCE_FILES_THREEWAYANNOTCOMPARE       ${SOURCE_FILES_ANNOTQ} ${SOURCE_FILES_BASIC} src/annotationQuery/ThreeWayAnnotCompare.cc)
set(SOURCE_FILES_TWOWAYANNOTCOMPARE         ${SOURCE_FILES_ANNOTQ} ${SOURCE_FILES_BASIC} src/annotationQuery/TwoWayAnnotCompare.cc)
set(SOURCE_FILES_TRANSCRIPTINFO             ${SOURCE_FILES_ANNOTQ} ${SOURCE_FILES_BASIC} src/annotationQuery/TranscriptInfo.cc) 

add_executable(CompareAnnotWithRef        ${SOURCE_FILES_COMPAREANNOTWITHREF})
add_executable(GetIOSingleExonTranscripts ${SOURCE_FILES_GETIOSINGLEEXONTRANSCRIPTS})
add_executable(Maf2Satsuma                ${SOURCE_FILES_MAF2SATSUMA})
add_executable(RunAnnotationQuery         ${SOURCE_FILES_RUNANNOTATIONQUERY})
add_executable(RunGTFCompare              ${SOURCE_FILES_RUNGTFCOMPARE})
add_executable(ThreeWayAnnotCompare       ${SOURCE_FILES_THREEWAYANNOTCOMPARE})
add_executable(TwoWayAnnotCompare         ${SOURCE_FILES_TWOWAYANNOTCOMPARE})
add_executable(TranscriptInfo             ${SOURCE_FILES_TRANSCRIPTINFO})

# kraken binaries
set(SOURCE_FILES_ASSIGNKRAKENIDS ${SOURCE_FILES_BASIC} src/kraken/AssignKrakenIDs.cc) 
set(SOURCE_FILES_CLEANKRAKENFILES ${SOURCE_FILES_BASIC} src/kraken/CleanKrakenFile.cc) 
set(SOURCE_FILES_KRAKENEVALUATOR  ${SOURCE_FILES_ANNOTQ} ${SOURCE_FILES_BASIC} ${SOURCE_FILES_COLA} ${SOURCE_FILES_KRAKEN} ${SOURCE_FILES_FFT} src/annotationQuery/MultiAlignParser.cc src/kraken/KrakenEvaluator.cc) 
set(SOURCE_FILES_RUNKRAKEN       ${SOURCE_FILES_ANNOTQ} ${SOURCE_FILES_BASIC} ${SOURCE_FILES_COLA} ${SOURCE_FILES_KRAKEN} ${SOURCE_FILES_FFT} src/annotationQuery/MultiAlignParser.cc src/kraken/GTFTransfer.cc src/kraken/RunKraken.cc) 

add_executable(AssignKrakenIDs         ${SOURCE_FILES_ASSIGNKRAKENIDS})
add_executable(CleanKrakenFiles        ${SOURCE_FILES_CLEANKRAKENFILES})
add_executable(KrakenEvaluator         ${SOURCE_FILES_KRAKENEVALUATOR})
add_executable(RunKraken               ${SOURCE_FILES_RUNKRAKEN})


