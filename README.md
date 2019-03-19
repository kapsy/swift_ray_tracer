# About

This program was made to accompany _Graphics like Pixar using Swift_, a lightning talk given at _try! Swift Tokyo_ in 2019.

This is a simple ray tracer, or more specifically, path tracer, built using Swift.

Much of the code is derived from Pete Shirley's Ray Tracing Minibooks, with a few added extras.

# Overview

For each pixel, `ns` rays are stochastically fired into a simple scene. 

The final pixel color is just the average color of each ray, which is calculated using the material properties of whatever sphere the ray happens to intersect with.

To make things interesting, a simple looping camera animation was adding using Rodrigues rotation. Single frames are output as `ppm` files, and converted to an `mp4` using `ffmpeg` after the program has run.

# Algorithms

- Line sphere intersection
- Path tracing
- Lamberian scattering
- Fresnel equations
- Schlick approximation
- Perlin noise
- Rodrigues rotation formula

# Requirements

- Xcode 10
- Swift 4.2
- ffmpeg

# Installation

Run `build.sh` and everything should just work!

# Performance

The default render settings are quite slow. The whole program is single threaded, and not optimized at all.

To reduce the rendering time, try the following:

- Reducing the output image resolution.
- Reduce the number of rays per pixel.
- Reduce the frames drawn, either the duration of the animation or framerate.
