#!/bin/bash

CODE_PATH=code
BUILD_PATH=build
ASSETS_PATH=data

OPTIMIZED_SWITCHES="-fno-builtin -O2 -ffast-math -ftrapping-math"
DEBUG_SWITCHES=""

IGNORE_WARNING_FLAGS="-Wno-unused-function -Wno-unused-variable -Wno-missing-braces -Wno-c++11-compat-deprecated-writable-strings"
OSX_DEPENDENCIES="-framework Cocoa -framework IOKit -framework CoreAudio -framework AudioToolbox"

# clang -g $DEBUG_SWITCHES -Wall $IGNORE_WARNING_FLAGS -lstdc++ -DINTERNAL $CODE_PATH/swift_cmd_line_test.cc -o $BUILD_PATH/swift_cmd_line_test
# clang++ -S -mllvm --x86-asm-syntax=intel $CODE_PATH/swift_cmd_line_test.cc -o $(BUILD_PATH)/swift_cmd_line_test

swiftc -g $DEBUG_SWITCHES $CODE_PATH/swift_cmd_line_test.swift -o $BUILD_PATH/swift_cmd_line_test

rm temp/out_001.ppm
$BUILD_PATH/swift_cmd_line_test
open temp/out_001.ppm

# so will need to output files directly from swift for this to work...
# would like to avoid this but not sure if we can...
# ffmpeg -framerate 25 -i image-%05d.jpg -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4

exit 0
