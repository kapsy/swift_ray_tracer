#!/bin/bash

CODE_PATH=code
BUILD_PATH=build
ASSETS_PATH=data

# OPTIMIZE_SWITCHES=""
OPTIMIZE_SWITCHES="-O"

swiftc -g $OPTIMIZE_SWITCHES $CODE_PATH/Math.swift $CODE_PATH/main.swift -o $BUILD_PATH/main

rm ./temp/*
rm ./out.mp4
$BUILD_PATH/main
# open temp/out_001.ppm
ffmpeg -i ./temp/out_%03d.ppm -c:v libx264 -crf 18 -preset slow -pix_fmt yuv420p -c:a copy out.mp4

exit 0
