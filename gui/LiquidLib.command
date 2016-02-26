#!/usr/bin/env bash

DIRNAME=`dirname "$0"`
cd $DIRNAME

python LiquidLib.gui.py

&& exit