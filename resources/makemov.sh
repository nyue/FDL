#!/bin/bash
path_to_files=../build
ffmpeg -i ${path_to_files}/density_export_%04d.png -b 128000 density.mov