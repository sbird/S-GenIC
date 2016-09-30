#!/bin/bash
#Checkout the submodules before building.
git submodule init
git submodule update
cd GadgetReader
git submodule init
git submodule update
