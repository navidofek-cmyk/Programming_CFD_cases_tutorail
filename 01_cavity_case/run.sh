#!/usr/bin/env bash

set -euo pipefail

g++ -std=c++20 -O2 -Wall -Wextra -pedantic main.cpp -o cavity
./cavity
