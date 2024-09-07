#!/bin/bash

script_dir=$1
for file in ${script_dir}/*; do time bash $file; done

